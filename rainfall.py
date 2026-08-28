# -*- coding: utf-8 -*-
"""Rainfall intensity (I30) for Q*, from Synoptic gauges, MRMS radar, or a constant.

The Synoptic transport and the irregular-obs -> 1-minute interpolation are ported
from ``stormscape.gauges`` (S. W. McCoy, MIT), which in turn adapted them from the
USGS **FlowAlert** package (King, Rengers, Wedell & Fee, 2024; CC0). Two changes
over the previous implementation, both from stormscape's method:

* request Synoptic's **derived precip** service (``precip=1``) and read the
  standardised ``precip_intervals_set_1d`` (mm since the previous ob) /
  ``precip_accumulated_set_1d`` (cumulative mm), instead of guessing among the
  many raw ``precip_*`` variables (which differ gauge to gauge); and
* interpolate the irregular gauge record onto a regular **1-minute** grid, then
  run rolling-window peak estimators on it -- the 15-minute peak via the same
  ``(i16 + i14) / 2`` estimator MRMS uses, so gauge and radar I15 are defined
  identically. QC is requested and bad values removed
  (``qc=on, qc_flags=on, qc_remove_data=on``).

Backwards compatibility: :func:`synoptic_i30(station_id, date, cfg)` keeps its old
signature and fallback behaviour. New: :func:`synoptic_i30_bbox` discovers gauges
in an AOI (no station id needed), and :func:`rainfall_i30` is a single dispatcher
the pipeline calls with ``cfg.rainfall_source`` in ``{synoptic, mrms, constant}``.
"""
from __future__ import annotations

import datetime as _dt

import numpy as np

SYNOPTIC_BASE = "https://api.synopticdata.com/v2"
META_URL = f"{SYNOPTIC_BASE}/stations/metadata"
TS_URL = f"{SYNOPTIC_BASE}/stations/timeseries"

# Synoptic derived-precip variables requested with precip=1 (units mm).
ACC_VAR = "precip_accumulated_set_1d"     # cumulative precip over the request
INT_VAR = "precip_intervals_set_1d"       # precip since the previous observation


# --------------------------------------------------------------------------- #
# token
# --------------------------------------------------------------------------- #
def synoptic_token(cfg):
    """Return a usable Synoptic token.

    Prefers exchanging the master API key for a fresh token via the /v2/auth
    endpoint (recommended by Synoptic); falls back to any static token given.
    """
    import requests
    if getattr(cfg, "synoptic_api_key", ""):
        try:
            a = requests.get(f"{SYNOPTIC_BASE}/auth",
                             params=dict(apikey=cfg.synoptic_api_key), timeout=30)
            tok = a.json().get("TOKEN")
            if tok:
                return tok
        except Exception:
            pass
    return getattr(cfg, "synoptic_token", "")


def _raise_for_synoptic(r):
    """Surface Synoptic HTTP errors with their human-readable message.

    A revoked token or an account without an active data plan returns HTTP 403
    ``"Unauthorized. Review your account access settings in the customer
    console."`` -- propagate that text rather than a downstream KeyError.
    """
    if r.status_code == 200:
        return
    try:
        body = r.json()
        msg = body.get("message") or body.get("ERROR") or str(body)
    except Exception:
        msg = r.text[:200]
    raise RuntimeError(f"Synoptic HTTP {r.status_code}: {msg}")


# --------------------------------------------------------------------------- #
# transport (ported from stormscape.gauges / FlowAlert)
# --------------------------------------------------------------------------- #
def get_stations(aoi, cfg, status="ACTIVE", pad_deg=0.0):
    """Synoptic stations within the AOI bbox -> list of dicts.

    Each dict has ``station_id, name, lon, lat, elevation, distance_deg`` (the
    last filled by :func:`_nearest`). ``aoi`` is anything
    :func:`dominantflowtype.aoi.load_aoi` accepts. Returns ``[]`` on zero
    results (RESPONSE_CODE 2).
    """
    import requests
    from .aoi import load_aoi
    bounds = load_aoi(aoi, pad_deg=pad_deg)
    bbox = ",".join(f"{v:.6f}" for v in bounds)      # W,S,E,N (Synoptic order)
    token = synoptic_token(cfg)
    if not token:
        raise RuntimeError("no Synoptic token/api-key supplied")
    r = requests.get(META_URL, timeout=60,
                     params=dict(token=token, bbox=bbox, complete="1"))
    _raise_for_synoptic(r)
    payload = r.json()
    if payload.get("SUMMARY", {}).get("RESPONSE_CODE") == 2 \
            or "STATION" not in payload:
        return []
    out = []
    for st in payload["STATION"]:
        try:
            lon = float(st["LONGITUDE"])
            lat = float(st["LATITUDE"])
        except (KeyError, TypeError, ValueError):
            continue
        stt = str(st.get("STATUS", "")).upper()
        if status and stt and stt != status.upper():
            continue
        out.append(dict(station_id=str(st["STID"]), name=st.get("NAME"),
                        lon=lon, lat=lat,
                        elevation=st.get("ELEVATION"), network=st.get("MNET_ID")))
    return out


def get_rainfall(station_ids, starttime, endtime, cfg, token=None):
    """Per-station precip time series from Synoptic over ``[start, end]`` (UTC).

    Returns ``{station_id: {'date_time': [...], INT_VAR: [...], ACC_VAR: [...]}}``.
    Over-large requests ("Querying too many station hours") are split in half and
    recombined. Ported from stormscape.gauges.get_rainfall.
    """
    import requests
    token = token or synoptic_token(cfg)
    fmt = "%Y%m%d%H%M"
    params = {
        "start": _utc(starttime).strftime(fmt),
        "end": _utc(endtime).strftime(fmt),
        "stid": ",".join(station_ids),
        "obtimezone": "UTC", "output": "json", "precip": "1",
        "qc": "on", "qc_flags": "on", "qc_remove_data": "on",
        "units": "metric,precip|mm", "token": token,
    }
    r = requests.get(TS_URL, params=params, timeout=120)
    _raise_for_synoptic(r)
    payload = r.json()
    if "STATION" in payload:
        out = {}
        for st in payload["STATION"]:
            obs = st.get("OBSERVATIONS", {}) or {}
            out[str(st["STID"])] = obs
        return out
    msg = payload.get("SUMMARY", {}).get("RESPONSE_MESSAGE", "")
    if "Querying too many station hours" in msg:     # split and recurse
        mid = starttime + (endtime - starttime) / 2
        first = get_rainfall(station_ids, starttime, mid, cfg, token=token)
        second = get_rainfall(station_ids, mid, endtime, cfg, token=token)
        for sid, d in second.items():
            first.setdefault(sid, {})
            for k, v in d.items():
                first[sid].setdefault(k, [])
                first[sid][k] = list(first[sid][k]) + list(v)
        return first
    raise RuntimeError(msg or "empty Synoptic response")


def _utc(t):
    """Coerce a datetime to UTC (naive datetimes are assumed already UTC)."""
    if t.tzinfo is None:
        return t.replace(tzinfo=_dt.timezone.utc)
    return t.astimezone(_dt.timezone.utc)


# --------------------------------------------------------------------------- #
# gauge -> per-minute -> peak intensity (ported from stormscape.gauges)
# --------------------------------------------------------------------------- #
def precipitation_per_minute(obs, precip_var=INT_VAR, time_var="date_time"):
    """Irregular gauge obs -> 1-minute per-minute precip (mm).

    Interpolates the cumulative gauge record onto a regular 1-minute grid, so
    tipping-bucket records reported at uneven intervals become a clean
    minute-resolution series the rolling-window estimators can run on. Ported
    from FlowAlert's ``rainfall.precipitation_per_minute``. Returns ``(times,
    cumulative, per_minute)``; empties for < 2 valid obs.
    """
    if not obs or precip_var not in obs or time_var not in obs:
        return [], np.array([]), np.array([])
    raw_t = obs[time_var]
    raw_v = obs[precip_var]
    tv = []
    for t, v in zip(raw_t, raw_v):
        if v is None:
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        tv.append((_parse_iso(t), max(fv, 0.0)))
    if len(tv) < 2:
        return [], np.array([]), np.array([])
    tv.sort(key=lambda p: p[0])
    times = np.array([p[0].timestamp() for p in tv])
    precip = np.array([p[1] for p in tv], dtype="float64")
    cumulative = np.cumsum(precip)
    grid = np.arange(times[0], times[-1] + 60, 60)
    interp = np.interp(grid, times, cumulative)
    per_minute = np.concatenate(([0.0], np.diff(interp)))
    out_times = [_dt.datetime.fromtimestamp(t, tz=_dt.timezone.utc) for t in grid]
    return out_times, interp, per_minute


def _parse_iso(t):
    return _dt.datetime.fromisoformat(str(t).replace("Z", "+00:00"))


def _peak_window(per_minute, minutes):
    """Peak trailing ``minutes``-window intensity (mm/h) over a 1-min series."""
    if len(per_minute) < minutes:
        return float("nan")
    sums = np.convolve(per_minute, np.ones(minutes), mode="valid")
    return float(np.nanmax(sums) * 60.0 / minutes)


def _peak_i15(per_minute):
    """Peak 15-min intensity via the ``(i16 + i14) / 2`` estimator (MRMS-matching)."""
    n = len(per_minute)
    if n < 16:
        return _peak_window(per_minute, 15)
    s16 = np.convolve(per_minute, np.ones(16), mode="valid")   # len n-15
    s14 = np.convolve(per_minute, np.ones(14), mode="valid")   # len n-13
    i16 = s16 * 60.0 / 16.0
    i14 = s14[2:2 + len(s16)] * 60.0 / 14.0
    return float(np.nanmax((i16 + i14) / 2.0))


def peak_intensity(obs, window_min=30):
    """Peak ``window_min``-minute intensity (mm/h) from one gauge's obs.

    Uses the derived ``precip_intervals`` variable interpolated to 1 minute; the
    15-minute window uses the radar-matching ``(i16+i14)/2`` estimator, other
    windows use plain trailing sums. NaN if the gauge has no usable precip.
    """
    _, _, per_minute = precipitation_per_minute(obs)
    if len(per_minute) < 2:
        return float("nan")
    return _peak_i15(per_minute) if window_min == 15 \
        else _peak_window(per_minute, window_min)


# --------------------------------------------------------------------------- #
# public: single-station I30 (backwards-compatible)
# --------------------------------------------------------------------------- #
def synoptic_i30(station_id, date_str, cfg, query_type="timeseries", tz="utc"):
    """Peak 30-min rainfall intensity I30 [mm/hr] at ONE Synoptic station.

    ``date_str`` = 'YYYY-MM-DD' (event day, UTC); queries that whole day.
    Returns ``(i30_mm_hr, meta)``. Falls back to ``cfg.default_i30_mm_hr`` on any
    error (bad token, 403 no-access, no precip variable, ...).
    """
    meta = {"station": station_id, "date": date_str, "source": "synoptic",
            "query_type": query_type}
    try:
        day = _dt.datetime.strptime(date_str, "%Y-%m-%d")
        end = day + _dt.timedelta(days=1)
        obs = get_rainfall([station_id], day, end, cfg)
        stn = obs.get(str(station_id)) or (list(obs.values())[0] if obs else None)
        if not stn or INT_VAR not in stn:
            raise RuntimeError(
                "no derived precip (precip_intervals_set_1d) for this station")
        i30 = peak_intensity(stn, cfg.rain_window_duration_min)
        if not np.isfinite(i30):
            raise RuntimeError("gauge reported no usable precip in the window")
        meta.update(variable=INT_VAR, n_obs=len(stn.get("date_time", [])),
                    i30_mm_hr=i30)
        return i30, meta
    except Exception as e:
        return _fallback(cfg, meta, e)


def synoptic_i30_bbox(aoi, date_str, cfg, reduce="max", point=None,
                      status="ACTIVE", pad_deg=0.05):
    """Peak I30 [mm/hr] from ALL Synoptic gauges in an AOI (no station id needed).

    Discovers gauges in the AOI, pulls each gauge's precip for the event day, and
    reduces the per-gauge peak I30 by ``reduce``:

    * ``'max'``     -- the wettest gauge (storm core; default);
    * ``'nearest'`` -- the gauge closest to ``point=(lon, lat)`` with usable
      precip (e.g. the basin outlet);
    * ``'mean'``    -- mean over gauges that reported precip.

    Returns ``(i30_mm_hr, meta)``; falls back to ``cfg.default_i30_mm_hr``.
    """
    meta = {"aoi": True, "date": date_str, "source": "synoptic_bbox",
            "reduce": reduce}
    try:
        day = _dt.datetime.strptime(date_str, "%Y-%m-%d")
        end = day + _dt.timedelta(days=1)
        stations = get_stations(aoi, cfg, status=status, pad_deg=pad_deg)
        if not stations:
            raise RuntimeError("no Synoptic stations in the AOI")
        ids = [s["station_id"] for s in stations]
        rain = {}
        for i in range(0, len(ids), 50):            # chunk to limit station-hours
            rain.update(get_rainfall(ids[i:i + 50], day, end, cfg))
        results = []
        for s in stations:
            i30 = peak_intensity(rain.get(s["station_id"], {}),
                                 cfg.rain_window_duration_min)
            if np.isfinite(i30):
                d = (float("nan") if point is None else
                     ((s["lon"] - point[0]) ** 2 + (s["lat"] - point[1]) ** 2) ** 0.5)
                results.append((s, i30, d))
        if not results:
            raise RuntimeError("no AOI gauge reported usable precip")
        if reduce == "nearest" and point is not None:
            s, i30, _ = min(results, key=lambda r: r[2])
        elif reduce == "mean":
            i30 = float(np.mean([r[1] for r in results]))
            s = max(results, key=lambda r: r[1])[0]
        else:                                        # max
            s, i30, _ = max(results, key=lambda r: r[1])
        meta.update(n_gauges=len(results), station=s["station_id"],
                    station_name=s.get("name"), i30_mm_hr=i30)
        return i30, meta
    except Exception as e:
        return _fallback(cfg, meta, e)


def _fallback(cfg, meta, exc):
    meta["error"] = str(exc)
    meta["i30_mm_hr"] = cfg.default_i30_mm_hr
    meta["fallback"] = True
    print(f"  [rainfall] Synoptic query failed ({exc}); using fallback "
          f"I30={cfg.default_i30_mm_hr} mm/hr (Cavagnaro optimised constant).")
    return cfg.default_i30_mm_hr, meta


# --------------------------------------------------------------------------- #
# dispatcher used by the pipeline
# --------------------------------------------------------------------------- #
def rainfall_i30(cfg, station_id=None, date_str=None, aoi=None, point=None):
    """Resolve I30 [mm/hr] according to ``cfg.rainfall_source``.

    * ``'mrms'``     -- gauge-free radar (needs ``aoi`` + ``date_str``); reduce
      and optional ``point`` come from ``cfg.mrms_reduce`` / ``point``. No
      account required.
    * ``'synoptic'`` -- ``station_id`` if given, else all gauges in ``aoi``
      (``synoptic_i30_bbox``). Needs a Synoptic account with data access.
    * ``'constant'`` -- ``cfg.default_i30_mm_hr`` (Cavagnaro's 9.7 mm/hr).

    Returns ``(i30_mm_hr, meta)``.
    """
    src = getattr(cfg, "rainfall_source", "synoptic")
    if src == "constant" or not date_str:
        return cfg.default_i30_mm_hr, {
            "source": "constant", "i30_mm_hr": cfg.default_i30_mm_hr,
            "note": "rainfall_source=constant or no event date supplied"}
    if src == "mrms":
        if aoi is None:
            return _fallback(cfg, {"source": "mrms"},
                             RuntimeError("MRMS needs an AOI (watershed)"))
        from . import mrms
        try:
            return mrms.mrms_i30(
                aoi, date_str, reduce=getattr(cfg, "mrms_reduce", "areal_max"),
                point=point, duration_min=cfg.rain_window_duration_min)
        except Exception as e:                       # noqa: BLE001
            return _fallback(cfg, {"source": "mrms"}, e)
    # synoptic (default)
    if station_id:
        return synoptic_i30(station_id, date_str, cfg)
    if aoi is not None:
        return synoptic_i30_bbox(aoi, date_str, cfg,
                                 reduce="nearest" if point else "max",
                                 point=point)
    return cfg.default_i30_mm_hr, {
        "source": "fallback", "i30_mm_hr": cfg.default_i30_mm_hr,
        "note": "no station id and no AOI supplied"}
