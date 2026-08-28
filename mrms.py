# -*- coding: utf-8 -*-
"""MRMS radar rainfall -> peak 15/30/60-minute intensity fields (and a scalar I30).

Ported from ``stormscape.mrms`` (S. W. McCoy, MIT; the i15 stacking estimator is
D. Cavagna's ``MRMS_stack.py``, used with permission) and adapted to this
package's **GDAL** stack -- no rasterio/xarray. It reads NOAA MRMS grib2 straight
from the public S3 bucket with GDAL's GRIB driver (present in the QGIS 4.x GDAL
build), so it needs **no account and no extra installs**.

Why this matters here
---------------------
Cavagnaro's Q* needs a 30-minute rainfall intensity ``I30``. The workflow's
Synoptic path requires a gauge near the basin *and* an account with data access;
for the 2024 Bear Fire that account currently returns HTTP 403 (see
``SESSION_SUMMARY.md``). MRMS gives a gridded, gauge-independent ``I30`` over the
basin from public radar -- :func:`mrms_i30` returns a single number the pipeline
drops straight into Q*.

i15 estimator
-------------
MRMS ``PrecipRate`` is a 2-min instantaneous rate (mm/h); ``a2 = rate * 2/60`` is
the 2-min accumulation (mm). Over a trailing 16-min window (8 steps):
``i16 = sum(8) * 60/16`` and ``i14 = sum(last 7) * 60/14``; ``i15 = mean(i16,
i14)``. The running maximum over the storm gives ``i15max``. The 30- and 60-min
peaks use plain trailing windows (15 / 30 steps) scaled to mm/h.

Storm-window detection
----------------------
Hourly ``RadarOnly_QPE_01H`` is scanned over the UTC window covering the local
calendar day; the wettest hours over the AOI (> ``qpe_thresh``, capped at
``max_wet_hours``) are kept and 2-min ``PrecipRate`` is stacked over each
contiguous run (with a 14-min lead so the rolling i15 is defined from the run's
first wet minute).
"""
from __future__ import annotations

import datetime as _dt
import gzip
import io
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from ._compat import gdal
from .aoi import load_aoi

S3 = "https://noaa-mrms-pds.s3.amazonaws.com/CONUS"

# CONUS MRMS grid (from grib metadata): UL (-130, 55), 0.01 deg, 7000 x 3500.
G_W, G_N, G_RES, G_NX, G_NY = -130.0, 55.0, 0.01, 7000, 3500

# (product subdirectory, file prefix); S3 key = <dir>/<YYYYMMDD>/<prefix>_<dt>
PRODUCTS = {
    "PrecipRate": ("PrecipRate_00.00", "MRMS_PrecipRate_00.00"),
    "RadarOnly":  ("RadarOnly_QPE_01H_00.00", "MRMS_RadarOnly_QPE_01H_00.00"),
    "RQI":        ("RadarQualityIndex_00.00", "MRMS_RadarQualityIndex_00.00"),
    "MultiSensor":  ("MultiSensor_QPE_01H_Pass2_00.00",
                     "MRMS_MultiSensor_QPE_01H_Pass2_00.00"),
    "MultiSensor1": ("MultiSensor_QPE_01H_Pass1_00.00",
                     "MRMS_MultiSensor_QPE_01H_Pass1_00.00"),
}

# Defaults (overridable per call).
QPE_THRESH = 2.5            # mm; hourly areal-max above this = a "wet hour"
MAX_WET_HRS = 8             # cap processed wet hours (cost + i15-peak capture)
SCAN_PAD_H = (4, 10)        # UTC scan = [day 04:00, next-day 10:00] ~ local day
WORKERS = 12                # parallel MRMS downloads


class Missing(Exception):
    """File genuinely absent on the server (HTTP 404) -- do not retry."""


# --------------------------------------------------------------------------- #
# dates / windows
# --------------------------------------------------------------------------- #
def parse_date(date):
    """Accept a date/datetime, 'YYYYMMDD', or 'YYYY-MM-DD' -> datetime.date."""
    if isinstance(date, _dt.datetime):
        return date.date()
    if isinstance(date, _dt.date):
        return date
    s = str(date).strip().replace("-", "")
    return _dt.datetime.strptime(s, "%Y%m%d").date()


def aoi_window(bounds):
    """Pixel window ``(c0, r0, nx, ny)`` + GDAL geotransform for a WGS84 AOI.

    ``bounds`` is ``(W, S, E, N)`` in lon/lat (stormscape convention).
    """
    w, s, e, n = bounds
    c0 = max(int(np.floor((w - G_W) / G_RES)), 0)
    c1 = min(int(np.ceil((e - G_W) / G_RES)), G_NX)
    r0 = max(int(np.floor((G_N - n) / G_RES)), 0)
    r1 = min(int(np.ceil((G_N - s) / G_RES)), G_NY)
    nx, ny = c1 - c0, r1 - r0
    # top-left corner of the window, GDAL geotransform (north-up).
    gt = (G_W + c0 * G_RES, G_RES, 0.0, G_N - r0 * G_RES, 0.0, -G_RES)
    return (c0, r0, nx, ny), gt


# --------------------------------------------------------------------------- #
# transport: fetch one grib2.gz from S3 -> AOI-windowed numpy array via GDAL
# --------------------------------------------------------------------------- #
def _read_grib_window(raw_bytes, win):
    """Read band 1 of an in-memory grib2 over a pixel window -> float32 array."""
    c0, r0, nx, ny = win
    vp = f"/vsimem/mrms_{os.getpid()}_{id(raw_bytes)}.grib2"
    gdal.FileFromMemBuffer(vp, raw_bytes)
    try:
        ds = gdal.Open(vp)
        if ds is None:
            raise RuntimeError("GDAL could not open MRMS grib2 (GRIB driver?)")
        a = ds.GetRasterBand(1).ReadAsArray(c0, r0, nx, ny).astype("float32")
        ds = None
    finally:
        gdal.Unlink(vp)
    a[a < 0] = np.nan          # MRMS no-coverage / missing flags
    return a


def fetch(product, t, win):
    """Download one MRMS grib2 -> AOI-windowed array (values < 0 -> NaN).

    Retries only transient failures. A 404 means the timestep does not exist;
    we raise :class:`Missing` at once so absent files never burn backoff time.
    """
    import requests
    date, hms = t.strftime("%Y%m%d"), t.strftime("%H%M%S")
    pdir, prefix = PRODUCTS[product]
    url = f"{S3}/{pdir}/{date}/{prefix}_{date}-{hms}.grib2.gz"
    err = "?"
    for k in range(3):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                raw = gzip.GzipFile(fileobj=io.BytesIO(r.content)).read()
                return _read_grib_window(raw, win)
            if r.status_code == 404:
                raise Missing(f"{product} {date}-{hms}")
            err = f"HTTP {r.status_code}"
        except Missing:
            raise
        except Exception as e:                     # noqa: BLE001
            err = repr(e)[:120]
        time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"{product} {date}-{hms}: {err}")


def fetch_many(product, times, win, workers=WORKERS):
    """Parallel fetch; returns ``{t: array}`` (absent/failed timesteps omitted)."""
    out = {}

    def one(t):
        try:
            return t, fetch(product, t, win)
        except Exception:                          # noqa: BLE001
            return t, None

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for t, a in ex.map(one, list(times)):
            if a is not None:
                out[t] = a
    return out


# --------------------------------------------------------------------------- #
# estimators + storm-window detection
# --------------------------------------------------------------------------- #
def compute_i15(stack):
    """i15 (mm/h) from a trailing list of >=8 2-minute accumulations (mm)."""
    s = np.dstack(stack[-8:])
    i16 = np.nansum(s, axis=2) * 60 / 16
    i14 = np.nansum(s[:, :, 1:], axis=2) * 60 / 14
    return (i16 + i14) / 2


def _time_range(t0, t1, step_min):
    """List of datetimes from t0..t1 inclusive at step_min-minute spacing."""
    out, t, step = [], t0, _dt.timedelta(minutes=step_min)
    while t <= t1:
        out.append(t)
        t += step
    return out


def find_wet_hours(date0, win, qpe_thresh=QPE_THRESH, max_wet_hours=MAX_WET_HRS,
                   scan_pad_h=SCAN_PAD_H, workers=WORKERS):
    """Scan hourly QPE; return ``(wet_hours, scan)`` as lists of ``(t, qmax)``.

    ``wet_hours`` is time-sorted and capped at ``max_wet_hours``; ``scan`` is the
    full record used to locate the single peak hour.
    """
    start = _dt.datetime(date0.year, date0.month, date0.day, scan_pad_h[0])
    n_hours = 24 + scan_pad_h[1] - scan_pad_h[0]
    hours = [start + _dt.timedelta(hours=h) for h in range(n_hours)]
    arrs = fetch_many("RadarOnly", hours, win, workers=workers)
    scan = [(t, float(np.nanmax(arrs[t])) if t in arrs
             and np.isfinite(arrs[t]).any() else 0.0) for t in hours]
    wet = [rc for rc in scan if rc[1] > qpe_thresh]
    wet.sort(key=lambda rc: rc[1], reverse=True)
    if not wet:                                    # fall back to the best hour
        wet = sorted(scan, key=lambda rc: rc[1], reverse=True)[:1]
    wet = sorted(wet[:max_wet_hours], key=lambda rc: rc[0])
    return wet, scan


def contiguous_runs(hour_list):
    """Group sorted wet hours into contiguous runs (a gap > 1 h splits)."""
    runs, cur = [], [hour_list[0]]
    for t in hour_list[1:]:
        if (t - cur[-1]) <= _dt.timedelta(hours=1):
            cur.append(t)
        else:
            runs.append(cur)
            cur = [t]
    runs.append(cur)
    return runs


# --------------------------------------------------------------------------- #
# storm-day peak-intensity fields
# --------------------------------------------------------------------------- #
def i_storm_day(aoi, date, pad_deg=0.05, qpe_thresh=QPE_THRESH,
                max_wet_hours=MAX_WET_HRS, scan_pad_h=SCAN_PAD_H,
                workers=WORKERS, verbose=True):
    """Peak i15/i30/i60 rainfall-intensity fields for one storm-day over an AOI.

    Returns
    -------
    dict
        ``fields`` ({name: 2-D ndarray} for i15max, i30max, i60max, i2max,
        total, tpki15), ``geotransform`` (GDAL 6-tuple), ``crs`` ('EPSG:4326'),
        ``win`` (c0,r0,nx,ny), and ``meta`` (scalar summary).
    """
    bounds = load_aoi(aoi, pad_deg=pad_deg)
    date0 = parse_date(date)
    win, gt = aoi_window(bounds)
    _, _, nx, ny = win
    if nx <= 0 or ny <= 0:
        raise ValueError(f"AOI {bounds} is empty or outside the CONUS grid.")
    shape = (ny, nx)

    wet, scan = find_wet_hours(date0, win, qpe_thresh, max_wet_hours,
                               scan_pad_h, workers)
    peak_t = max(scan, key=lambda rc: rc[1])[0]
    qmax = max(rc[1] for rc in scan)
    if verbose:
        print(f"  MRMS {date0}: peak {peak_t:%m-%d %H}Z qmax={qmax:.1f} mm, "
              f"{len(wet)} wet hr", flush=True)

    i15_max = np.zeros(shape, np.float32)
    i30_max = np.zeros(shape, np.float32)
    i60_max = np.zeros(shape, np.float32)
    i2_max = np.zeros(shape, np.float32)
    total = np.zeros(shape, np.float32)
    tpki15 = np.full(shape, np.nan, np.float32)

    for run in contiguous_runs([rc[0] for rc in wet]):
        t0 = run[0] - _dt.timedelta(minutes=14)    # lead so rolling i15 is valid
        t1 = run[-1] + _dt.timedelta(hours=1)
        steps = _time_range(t0, t1, 2)
        arrs = fetch_many("PrecipRate", steps, win, workers=workers)
        stack = []
        for t in steps:
            if t not in arrs:                      # missing timestep -> reset
                stack = []
                continue
            a = arrs[t]
            a2 = np.nan_to_num(np.clip(a, 0, None)) * 2 / 60
            i2_max = np.fmax(i2_max, a)
            total = total + a2
            stack.append(a2)
            stack = stack[-30:]                    # keep the longest (60-min) window
            if len(stack) >= 8:                    # i15 = mean(i16, i14)
                i15 = compute_i15(stack)
                newmax = i15 > i15_max
                i15_max = np.where(newmax, i15, i15_max)
                tpki15 = np.where(newmax, t.hour + t.minute / 60.0, tpki15)
            if len(stack) >= 15:                   # i30: trailing 30 min
                i30 = np.nansum(np.dstack(stack[-15:]), axis=2) * 60 / 30
                i30_max = np.fmax(i30_max, i30)
            if len(stack) >= 30:                   # i60: trailing 60 min
                i60 = np.nansum(np.dstack(stack[-30:]), axis=2) * 60 / 60
                i60_max = np.fmax(i60_max, i60)

    fields = {"i15max": i15_max, "i30max": i30_max, "i60max": i60_max,
              "i2max": i2_max, "total": total, "tpki15": tpki15}
    meta = dict(date=str(date0), peak_utc=f"{peak_t:%Y%m%d-%H%M}",
                qmax_mm=float(qmax), n_wet_hr=int(len(wet)),
                i15max_aoi=float(np.nanmax(i15_max)),
                i30max_aoi=float(np.nanmax(i30_max)),
                i60max_aoi=float(np.nanmax(i60_max)))
    return dict(fields=fields, geotransform=gt, crs="EPSG:4326",
                win=win, meta=meta)


def save_fields(result, out_dir, key, which=None):
    """Write storm-day fields to ``out_dir/<key>_<field>.tif`` (GDAL GTiff, LZW).

    Returns the list of written paths.
    """
    os.makedirs(out_dir, exist_ok=True)
    gt = result["geotransform"]
    drv = gdal.GetDriverByName("GTiff")
    from osgeo import osr
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    paths = []
    for label, arr in result["fields"].items():
        if which and label not in which:
            continue
        path = os.path.join(out_dir, f"{key}_{label}.tif")
        ny, nx = arr.shape
        ds = drv.Create(path, nx, ny, 1, gdal.GDT_Float32,
                        options=["COMPRESS=LZW"])
        ds.SetGeoTransform(gt)
        ds.SetProjection(srs.ExportToWkt())
        band = ds.GetRasterBand(1)
        band.WriteArray(arr.astype("float32"))
        band.SetNoDataValue(float("nan"))
        ds = None
        paths.append(path)
    return paths


# --------------------------------------------------------------------------- #
# scalar I30 for the Q* pipeline
# --------------------------------------------------------------------------- #
def _sample_point(field, gt, lon, lat):
    """Nearest-cell value of a field at a lon/lat (NaN if outside)."""
    col = int((lon - gt[0]) / gt[1])
    row = int((lat - gt[3]) / gt[5])
    ny, nx = field.shape
    if 0 <= row < ny and 0 <= col < nx:
        return float(field[row, col])
    return float("nan")


def mrms_i30(aoi, date, reduce="areal_max", point=None, duration_min=30,
             pad_deg=0.05, verbose=True, **kwargs):
    """Single peak rainfall intensity (mm/h) over an AOI for one storm-day.

    A gauge-free ``I30`` for Cavagnaro's Q*. Builds the MRMS peak-intensity
    fields (:func:`i_storm_day`) and reduces the chosen duration's field to one
    number:

    * ``reduce='areal_max'`` -- the maximum over the AOI (the storm core; use
      when the whole small basin is under one convective cell);
    * ``reduce='areal_mean'`` -- the AOI mean of the peak field;
    * ``reduce='point'`` -- the value at ``point=(lon, lat)`` (e.g. the basin
      outlet / thalweg head), nearest-cell.

    ``duration_min`` selects i15/i30/i60 (default 30 -> the I30 Q* wants).
    Returns ``(intensity_mm_hr, meta)``; ``meta`` carries the full field summary
    and the reduction used, so it can be logged alongside the Synoptic path.
    """
    field_key = {15: "i15max", 30: "i30max", 60: "i60max"}.get(duration_min)
    if field_key is None:
        raise ValueError("duration_min must be 15, 30, or 60")
    res = i_storm_day(aoi, date, pad_deg=pad_deg, verbose=verbose, **kwargs)
    field = res["fields"][field_key]
    if reduce == "point":
        if point is None:
            raise ValueError("reduce='point' needs point=(lon, lat)")
        val = _sample_point(field, res["geotransform"], point[0], point[1])
    elif reduce == "areal_mean":
        val = float(np.nanmean(field)) if np.isfinite(field).any() else float("nan")
    else:                                          # areal_max (default)
        val = float(np.nanmax(field)) if np.isfinite(field).any() else float("nan")
    meta = dict(source="mrms", duration_min=duration_min, reduce=reduce,
                i30_mm_hr=val, point=point, **res["meta"])
    if verbose:
        print(f"       MRMS I{duration_min} ({reduce}) = {val:.3g} mm/hr")
    return val, meta
