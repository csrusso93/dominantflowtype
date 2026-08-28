# -*- coding: utf-8 -*-
"""Raster sampling and DEM-of-Difference construction."""
from __future__ import annotations

import math
import os

import numpy as np

from ._compat import gdal, _HAS_QGIS, QgsCoordinateReferenceSystem


class RasterSampler:
    """Bilinear point sampler + profile extractor over a single-band raster.

    Loads the whole band into memory (post-fire survey DEMs of a single basin
    are small).  Handles nodata and out-of-bounds gracefully (returns NaN).
    """

    def __init__(self, path, unit="auto"):
        """`unit` ('m' | 'ft' | 'auto') is the DEM's VERTICAL unit; feet are
        converted to metres so every sampled elevation/depth is in metres.
        'auto' inspects the raster CRS and converts if it looks like feet."""
        self.path = path
        ds = gdal.Open(path)
        if ds is None:
            raise IOError(f"cannot open raster: {path}")
        self.ds = ds
        self.vscale, _vreason = _resolve_scale(ds, unit)
        self.gt = ds.GetGeoTransform()          # (x0, dx, 0, y0, 0, dy)
        self.inv_gt = gdal.InvGeoTransform(self.gt)
        self.nx = ds.RasterXSize
        self.ny = ds.RasterYSize
        band = ds.GetRasterBand(1)
        self.nodata = band.GetNoDataValue()
        arr = band.ReadAsArray().astype("float64")
        if self.nodata is not None:
            arr[arr == self.nodata] = np.nan
        # very large sentinel nodata (e.g. 3.4e38) also -> NaN
        arr[np.abs(arr) > 1e30] = np.nan
        # convert vertical units to metres AFTER nodata masking
        if self.vscale != 1.0:
            arr *= self.vscale
            print(f"  [raster] {os.path.basename(path)} vertical scaled to metres "
                  f"(x{self.vscale:.5f}; {_vreason})")
        self.arr = arr
        self.wkt = ds.GetProjection()
        self.crs = None
        if _HAS_QGIS and self.wkt:
            self.crs = QgsCoordinateReferenceSystem.fromWkt(self.wkt)
        self.xres = abs(self.gt[1])
        self.yres = abs(self.gt[5])

    def sample(self, x, y):
        """Bilinear-interpolated value at world (x, y); NaN if nodata/off-grid."""
        px, py = gdal.ApplyGeoTransform(self.inv_gt, x, y)
        px -= 0.5
        py -= 0.5                     # cell-centre convention
        x0 = int(math.floor(px))
        y0 = int(math.floor(py))
        if x0 < 0 or y0 < 0 or x0 + 1 >= self.nx or y0 + 1 >= self.ny:
            # fall back to nearest valid inside-grid pixel
            xi = min(max(int(round(px)), 0), self.nx - 1)
            yi = min(max(int(round(py)), 0), self.ny - 1)
            return float(self.arr[yi, xi])
        fx = px - x0
        fy = py - y0
        v00 = self.arr[y0, x0]
        v10 = self.arr[y0, x0 + 1]
        v01 = self.arr[y0 + 1, x0]
        v11 = self.arr[y0 + 1, x0 + 1]
        vals = np.array([v00, v10, v01, v11])
        if np.all(np.isnan(vals)):
            return float("nan")
        # if some corners nodata, use nearest available (nan-robust)
        if np.any(np.isnan(vals)):
            valid = ~np.isnan(vals)
            return float(vals[valid].mean())
        top = v00 * (1 - fx) + v10 * fx
        bot = v01 * (1 - fx) + v11 * fx
        return float(top * (1 - fy) + bot * fy)

    def profile(self, x0, y0, x1, y1, step):
        """Return (s, z) arrays sampled from (x0,y0) to (x1,y1) every `step` m."""
        length = math.hypot(x1 - x0, y1 - y0)
        n = max(int(math.ceil(length / step)) + 1, 2)
        s = np.linspace(0.0, length, n)
        xs = np.linspace(x0, x1, n)
        ys = np.linspace(y0, y1, n)
        z = np.array([self.sample(px, py) for px, py in zip(xs, ys)])
        return s, z, xs, ys


# US survey / international foot -> metre. The two differ by 2 ppm, negligible
# for post-event bed change; the international value is used.
_FT_TO_M = 0.3048


def _explicit_scale(unit):
    """Fixed metre-multiplier for an explicit unit, or None for 'auto'."""
    u = (unit or "auto").strip().lower()
    if u in ("ft", "feet", "foot", "us-ft", "usft", "us survey foot"):
        return _FT_TO_M
    if u in ("m", "meter", "metre", "meters", "metres"):
        return 1.0
    if u in ("auto", ""):
        return None
    raise ValueError(f"unrecognised vertical unit {unit!r}; use 'm', 'ft' or 'auto'")


def _auto_scale_from_srs(ds):
    """Best-effort metre-multiplier from a raster's CRS (for unit='auto').

    Prefers an explicit vertical-axis unit; if the raster has no vertical CRS,
    falls back to the horizontal linear unit as a heuristic (common for US State
    Plane feet DEMs where Z shares the horizontal foot unit). Returns
    (scale, reason). Defaults to metres when nothing indicates feet.
    """
    srs = None
    try:
        srs = ds.GetSpatialRef()
    except Exception:
        srs = None
    if srs is None:
        return 1.0, "no CRS metadata; assuming metres"
    # 1) explicit vertical axis unit (compound / 3D CRS)
    vname = srs.GetAttrValue("VERT_CS|UNIT") or srs.GetAttrValue("VERT_CS")
    if vname:
        if any(k in vname.lower() for k in ("foot", "feet", "ft")):
            return _FT_TO_M, f"vertical CRS unit '{vname}' -> feet"
        return 1.0, f"vertical CRS unit '{vname}' -> metres"
    # 2) no vertical CRS: use the horizontal linear unit as a hint
    try:
        hname = srs.GetLinearUnitsName() or ""
    except Exception:
        hname = ""
    if any(k in hname.lower() for k in ("foot", "feet", "ft")):
        return _FT_TO_M, (f"no vertical CRS; horizontal unit '{hname}' suggests "
                          f"feet (heuristic — override with *_vertical_unit='m')")
    return 1.0, "no vertical CRS unit found; assuming metres"


def _resolve_scale(ds, unit):
    """Return (metre_multiplier, reason) honouring explicit units or auto-detect."""
    s = _explicit_scale(unit)
    if s is not None:
        return s, f"unit='{unit}'"
    return _auto_scale_from_srs(ds)


def reproject_raster(src_path, out_path, dst_epsg, resample="bilinear"):
    """Warp a raster to ``EPSG:<dst_epsg>`` (GeoTIFF). Returns ``out_path``.

    Used to bring a full-coverage routing DEM (e.g. USGS 3DEP in EPSG:5070) into
    the working CRS so A_us flow routing and the IDFVA hydrology DEM share one CRS
    with the point-cloud-derived DTM and change raster.
    """
    gdal.Warp(out_path, src_path, dstSRS=f"EPSG:{dst_epsg}", resampleAlg=resample,
              format="GTiff", creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
    print(f"  [raster] reprojected {os.path.basename(src_path)} -> EPSG:{dst_epsg}")
    return out_path


def build_dod(pre_path, post_path, out_path, pre_unit="auto", post_unit="auto"):
    """DEM-of-Difference = post - pre, resampled onto the POST grid.

    Convention: positive = DEPOSITION (surface rose), negative = EROSION.
    `pre_unit` / `post_unit` ('m' | 'ft' | 'auto') give each DEM's vertical unit;
    a DEM in feet is scaled to metres before differencing so the DoD is always in
    metres. Returns the RasterSampler of the written DoD.
    """
    # Resolve vertical units from each DEM's OWN CRS before the warp reprojects
    # the pre-DEM (warping would overwrite pre's CRS with post's).
    pre_scale, pre_reason = _resolve_scale(gdal.Open(pre_path), pre_unit)
    post = gdal.Open(post_path)
    post_scale, post_reason = _resolve_scale(post, post_unit)
    gt = post.GetGeoTransform()
    nx, ny = post.RasterXSize, post.RasterYSize
    xmin = gt[0]
    ymax = gt[3]
    xmax = xmin + nx * gt[1]
    ymin = ymax + ny * gt[5]
    # warp pre onto exact post grid
    warped = gdal.Warp(
        "", pre_path, format="MEM",
        outputBounds=(xmin, ymin, xmax, ymax),
        width=nx, height=ny, resampleAlg="bilinear",
        dstSRS=post.GetProjection(),
    )
    pre_arr = warped.GetRasterBand(1).ReadAsArray().astype("float64")
    pre_nd = warped.GetRasterBand(1).GetNoDataValue()
    post_band = post.GetRasterBand(1)
    post_arr = post_band.ReadAsArray().astype("float64")
    post_nd = post_band.GetNoDataValue()
    for a, nd in ((pre_arr, pre_nd), (post_arr, post_nd)):
        if nd is not None:
            a[a == nd] = np.nan
        a[np.abs(a) > 1e30] = np.nan
    # convert each DEM to metres AFTER nodata masking (scaling nodata is meaningless)
    if pre_scale != 1.0:
        pre_arr *= pre_scale
        print(f"  [DoD] pre-DEM scaled to metres (x{pre_scale:.5f}; {pre_reason})")
    if post_scale != 1.0:
        post_arr *= post_scale
        print(f"  [DoD] post-DEM scaled to metres (x{post_scale:.5f}; {post_reason})")
    dod = post_arr - pre_arr
    _report_dod_sanity(dod, abs(gt[1]), abs(gt[5]))
    drv = gdal.GetDriverByName("GTiff")
    out = drv.Create(out_path, nx, ny, 1, gdal.GDT_Float32,
                     options=["COMPRESS=DEFLATE", "TILED=YES"])
    out.SetGeoTransform(gt)
    out.SetProjection(post.GetProjection())
    ob = out.GetRasterBand(1)
    ob.SetNoDataValue(-9999.0)
    dod_w = np.where(np.isnan(dod), -9999.0, dod).astype("float32")
    ob.WriteArray(dod_w)
    out.FlushCache()
    out = None
    # DoD values are already metres; force unit='m' so auto-detect does not
    # re-scale them when the inherited (e.g. State Plane feet) CRS says otherwise.
    return RasterSampler(out_path, unit="m")


# Bed change beyond this (metres) is physically implausible for a single
# post-fire event and usually signals a vertical units/datum mismatch.
DOD_PLAUSIBLE_ABS_M = 10.0
# A feet->metre confusion inflates every value by ~3.28x; a median |dz| this
# large over an aligned pair is the classic fingerprint.
DOD_FEET_SUSPECT_MEDIAN_M = 3.0


def _report_dod_sanity(dod, xres, yres, big_thresh=DOD_PLAUSIBLE_ABS_M):
    """Print a quick DoD sanity report and warn on likely datum/units mismatch.

    Convention (see build_dod): positive = deposition, negative = erosion.
    """
    v = dod[np.isfinite(dod)]
    n_all = int(dod.size)
    n_valid = int(v.size)
    print("  [DoD sanity] post - pre (positive = deposition, negative = erosion)")
    if n_valid == 0:
        print("    -> no overlapping valid cells (pre/post do not intersect?)")
        return {}
    cell_m2 = (xres * yres) or 1.0
    lo, hi = np.percentile(v, [1, 99])
    med_abs = float(np.median(np.abs(v)))
    n_big = int(np.count_nonzero(np.abs(v) > big_thresh))
    net_vol = float(np.sum(v) * cell_m2)                 # signed net volume [m^3]
    ero_vol = float(-np.sum(v[v < 0]) * cell_m2)         # erosion volume [m^3]
    dep_vol = float(np.sum(v[v > 0]) * cell_m2)          # deposition volume [m^3]
    print(f"    valid cells   = {n_valid:,} / {n_all:,} "
          f"({100.0 * n_valid / n_all:.1f}% overlap)")
    print(f"    min / max     = {float(v.min()):+.2f} / {float(v.max()):+.2f} m")
    print(f"    median |dz|   = {med_abs:.2f} m   "
          f"(1-99%: {lo:+.2f} .. {hi:+.2f} m)")
    print(f"    |dz| > {big_thresh:g} m   = {n_big:,} cells "
          f"({100.0 * n_big / n_valid:.2f}%)")
    print(f"    volume: erosion {ero_vol:,.0f} m3 | deposition {dep_vol:,.0f} m3 "
          f"| net {net_vol:+,.0f} m3")
    # heuristics ----------------------------------------------------------------
    if med_abs > DOD_FEET_SUSPECT_MEDIAN_M:
        print(f"    [WARN] median |dz| = {med_abs:.2f} m is large. If one DEM is "
              f"in FEET, dividing it by 3.28084 would give ~{med_abs / 3.28084:.2f} m.")
    if abs(float(np.median(v))) > 1.0:
        print(f"    [WARN] median dz = {float(np.median(v)):+.2f} m (not ~0): a "
              f"uniform bias suggests a vertical datum offset (e.g. geoid vs "
              f"ellipsoid) between the two DEMs.")
    if float(v.max()) > big_thresh or float(v.min()) < -big_thresh:
        print(f"    [WARN] extremes exceed +/-{big_thresh:g} m; check DEM vertical "
              f"units (m vs ft) and datum, and mask survey-edge artefacts.")
    return dict(valid_cells=n_valid, total_cells=n_all, min_m=float(v.min()),
                max_m=float(v.max()), median_dz_m=float(np.median(v)),
                median_abs_m=med_abs, p1_m=float(lo), p99_m=float(hi),
                n_beyond_thresh=n_big, erosion_vol_m3=ero_vol,
                deposition_vol_m3=dep_vol, net_vol_m3=net_vol)
