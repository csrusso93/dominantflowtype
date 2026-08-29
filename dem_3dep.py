# -*- coding: utf-8 -*-
"""USGS **3DEP** / The National Map DEMs for an AOI (stormscape's DEM import).

Ported from ``stormscape.dem`` (S. W. McCoy, MIT). 3DEP is the public elevation
layer behind The National Map, so this fetches the same bare-earth data the
interactive downloader serves, but scripted and AOI-clipped -- handy for a
full-coverage routing DEM when the local lidar tiles do not cover the whole
watershed.

Two access paths, same USGS 3DEP data, both account-free:

* :func:`get_dem_py3dep` -- the faithful stormscape path via **py3dep** (the
  HyRiver 3DEP client). Most capable (``dem_sources`` / ``coverage_fraction``
  availability checks) but needs ``py3dep`` + ``rioxarray`` installed
  (``pip install py3dep rioxarray`` in the OSGeo4W shell). Imported lazily, so
  this module loads without them.
* :func:`get_dem_rest` -- a **py3dep-free** fallback that GETs the 3DEP dynamic
  **ImageServer** (``exportImage``) and writes the GeoTIFF with GDAL. Works in a
  stock QGIS Python with no extra installs.

:func:`get_dem` tries py3dep, falls back to the REST path, and always returns a
path to a GeoTIFF on disk (this package works with file paths + ``RasterSampler``
/ ``FlowRouter``, not xarray).
"""
from __future__ import annotations

import os

from ._compat import gdal, osr
from .aoi import load_aoi

# 3DEP dynamic elevation ImageServer (bare-earth, metres).
IMAGESERVER = ("https://elevation.nationalmap.gov/arcgis/rest/services/"
               "3DEPElevation/ImageServer/exportImage")

# Resolutions 3DEP publishes, coarse->fine (metres).
STD_RESOLUTIONS = (60, 30, 10, 5, 3, 1)


# --------------------------------------------------------------------------- #
# py3dep path (faithful stormscape port)
# --------------------------------------------------------------------------- #
def _py3dep():
    try:
        import py3dep
    except ImportError as e:                              # pragma: no cover
        raise ImportError(
            "py3dep is required for this path. Install it with\n"
            "  pip install py3dep rioxarray   (in the OSGeo4W shell)\n"
            "or use get_dem_rest(), which needs no extra packages.") from e
    return py3dep


def dem_sources(aoi, res=None):
    """Query 3DEP source footprints over an AOI (which lidar projects cover it).

    Returns the GeoDataFrame from ``py3dep.query_3dep_sources`` (columns include
    ``dem_res`` like ``'1m'``/``'10m'``). ``res`` optionally filters to one or
    more resolution strings. Needs py3dep.
    """
    bounds = load_aoi(aoi)
    return _py3dep().query_3dep_sources(bounds, crs=4326, res=res)


def coverage_fraction(aoi, res="1m"):
    """Fraction of the AOI covered by 3DEP sources at resolution ``res`` (needs py3dep).

    Mirrors the 1 m-availability check used to decide whether a basin can be
    fetched at lidar resolution without seam artefacts from 10 m fill.
    """
    import geopandas as gpd
    from shapely.geometry import box
    bounds = load_aoi(aoi)
    geom = box(*bounds)
    try:
        src = _py3dep().query_3dep_sources(bounds, crs=4326, res=res)
    except Exception:                                  # noqa: BLE001
        return 0.0
    if src is None or not len(src):
        return 0.0
    src = src.to_crs(5070)
    g = gpd.GeoSeries([geom], crs=4326).to_crs(5070).iloc[0]
    inter = src.geometry.union_all().intersection(g)
    return float(inter.area / g.area) if g.area else 0.0


def get_dem_py3dep(aoi, resolution=10, dst_crs="EPSG:5070", out_path=None,
                   pad_deg=0.02):
    """Download a 3DEP DEM via py3dep, reproject, and write a GeoTIFF (needs py3dep).

    Returns ``out_path``. Faithful port of ``stormscape.dem.get_dem`` (minus the
    HyRiver-specific retry loop), writing to disk instead of returning xarray.
    """
    import rioxarray  # noqa: F401  (registers .rio)  # pragma: no cover
    bounds = load_aoi(aoi, pad_deg=pad_deg)
    from shapely.geometry import box
    dem = _py3dep().get_dem(box(*bounds), resolution=resolution, crs=4326)
    if dst_crs is not None:
        dem = dem.rio.reproject(dst_crs, resolution=resolution)
    dem = dem.where(dem > -1e4)                       # drop 3DEP fill/no-data
    dem.name = "elevation"
    if out_path is None:
        out_path = f"3dep_{resolution}m.tif"
    dem.rio.to_raster(out_path, compress="LZW")
    return out_path


# --------------------------------------------------------------------------- #
# py3dep-free REST path (GDAL only)
# --------------------------------------------------------------------------- #
def _bbox_to_crs(bounds, dst_epsg):
    """Reproject a ``(W,S,E,N)`` 4326 bbox to ``dst_epsg`` -> (minx,miny,maxx,maxy)."""
    w, s, e, n = bounds
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(int(dst_epsg))
    for sr in (src, dst):
        try:
            sr.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:                          # pragma: no cover
            pass
    ct = osr.CoordinateTransformation(src, dst)
    xs, ys = [], []
    for x, y in ((w, s), (w, n), (e, s), (e, n)):
        px, py, _ = ct.TransformPoint(x, y)
        xs.append(px)
        ys.append(py)
    return (min(xs), min(ys), max(xs), max(ys))


def get_dem_rest(aoi, resolution=10, dst_epsg=5070, out_path=None, pad_deg=0.02,
                 max_pixels=4096, timeout=180, retries=3, retry_wait=5,
                 verbose=True):
    """Download a 3DEP DEM from the ImageServer (``exportImage``) with GDAL only.

    Requests bare-earth elevation over the AOI in ``dst_epsg`` (default 5070,
    CONUS Albers metres, so pixels are square metres) at ``resolution`` metres,
    as a float32 GeoTIFF. The service caps a request at ~4100x15000 px; if the
    AOI/resolution would exceed ``max_pixels`` on a side the effective resolution
    is coarsened (and reported) rather than failing. The dynamic ImageServer is
    occasionally slow/overloaded (HTTP 5xx), so the request is retried up to
    ``retries`` times with a linear backoff. Returns ``out_path`` or ``None``.
    """
    import time
    import requests
    bounds = load_aoi(aoi, pad_deg=pad_deg)
    minx, miny, maxx, maxy = _bbox_to_crs(bounds, dst_epsg)
    wm, hm = maxx - minx, maxy - miny
    nx, ny = max(int(round(wm / resolution)), 1), max(int(round(hm / resolution)), 1)
    eff = resolution
    if max(nx, ny) > max_pixels:                     # coarsen to fit the cap
        eff = resolution * (max(nx, ny) / max_pixels)
        nx, ny = max(int(round(wm / eff)), 1), max(int(round(hm / eff)), 1)
        if verbose:
            print(f"  [3DEP] AOI too large for {resolution} m at the service cap; "
                  f"using ~{eff:.1f} m ({nx}x{ny} px).")
    if out_path is None:
        out_path = f"3dep_{int(round(eff))}m.tif"
    params = {
        "bbox": f"{minx},{miny},{maxx},{maxy}",
        "bboxSR": int(dst_epsg), "imageSR": int(dst_epsg),
        "size": f"{nx},{ny}", "format": "tiff", "pixelType": "F32",
        "noData": -9999, "interpolation": "RSP_BilinearInterpolation",
        "adjustAspectRatio": "false", "f": "image",
    }
    if verbose:
        print(f"  [3DEP] exportImage {nx}x{ny}px @ ~{eff:.1f} m, EPSG:{dst_epsg}")
    last = "?"
    for attempt in range(retries + 1):
        try:
            r = requests.get(IMAGESERVER, params=params, timeout=timeout)
            if r.status_code == 200 and r.content[:2] in (b"II", b"MM"):
                with open(out_path, "wb") as fh:
                    fh.write(r.content)
                # some ImageServer builds omit georeferencing tags; stamp them.
                _ensure_georef(out_path, minx, maxy, eff, dst_epsg)
                if verbose:
                    print(f"  [3DEP] saved -> {out_path}")
                return out_path
            last = f"HTTP {r.status_code}: {r.text[:120]}"
        except Exception as ex:                        # noqa: BLE001
            last = repr(ex)[:120]
        if attempt < retries:
            if verbose:
                print(f"  [3DEP] attempt {attempt + 1} failed ({last}); retrying ...")
            time.sleep(retry_wait * (attempt + 1))
    print(f"  [3DEP] download failed after {retries + 1} tries: {last}")
    return None


def _ensure_georef(path, ulx, uly, res, epsg):
    """Stamp geotransform/CRS onto a GeoTIFF if the service returned it unset."""
    ds = gdal.Open(path, gdal.GA_Update)
    if ds is None:
        return
    gt = ds.GetGeoTransform()
    if gt == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0) or not ds.GetProjection():
        ds.SetGeoTransform((ulx, res, 0.0, uly, 0.0, -res))
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        ds.SetProjection(srs.ExportToWkt())
    ds = None


def get_dem(aoi, resolution=10, out_path=None, dst_epsg=5070, pad_deg=0.02,
            prefer_py3dep=True, verbose=True):
    """Download a USGS 3DEP DEM for an AOI -> path to a GeoTIFF (account-free).

    Tries the faithful py3dep path first (if installed and ``prefer_py3dep``),
    then falls back to the GDAL ImageServer path so it always works in a stock
    QGIS Python. ``resolution`` is metres (1/3/5/10/30/60); 1 m needs lidar
    coverage over the AOI (check :func:`coverage_fraction` if py3dep is present).
    """
    if prefer_py3dep:
        try:
            return get_dem_py3dep(aoi, resolution=resolution,
                                  dst_crs=f"EPSG:{dst_epsg}", out_path=out_path,
                                  pad_deg=pad_deg)
        except Exception as e:                         # noqa: BLE001
            if verbose:
                print(f"  [3DEP] py3dep unavailable ({str(e)[:80]}); "
                      f"using the REST ImageServer path.")
    return get_dem_rest(aoi, resolution=resolution, dst_epsg=dst_epsg,
                        out_path=out_path, pad_deg=pad_deg, verbose=verbose)
