# -*- coding: utf-8 -*-
"""Area-of-interest helpers: normalise anything to a WGS84 bounding box.

Ported from ``stormscape.aoi`` (S. W. McCoy, MIT) and adapted to this package's
GDAL/QGIS stack (no geopandas/shapely dependency). The MRMS and 3DEP engines
only need the AOI as a geographic (EPSG:4326) bounding box; this module accepts
the four things a user is likely to have on hand -- a bbox tuple, an OGR-readable
vector file, a QGIS vector layer, or a QGIS extent+CRS -- and returns bounds in a
single, consistent convention.

Bounds convention
-----------------
All bounds are ``(west, south, east, north)`` in lon/lat degrees, matching
stormscape. (Note: :mod:`dominantflowtype.opentopo` uses the *other* ordering
``(s, w, n, e)``; convert with :func:`to_swne` at that boundary.)
"""
from __future__ import annotations

import os

from ._compat import gdal, ogr, osr


def pad_bounds(bounds, pad_deg):
    """Expand a ``(W, S, E, N)`` bbox by ``pad_deg`` degrees on every side."""
    w, s, e, n = bounds
    return (w - pad_deg, s - pad_deg, e + pad_deg, n + pad_deg)


def to_swne(bounds):
    """``(W, S, E, N)`` -> ``(S, W, N, E)`` for the OpenTopography helpers."""
    w, s, e, n = bounds
    return (s, w, n, e)


def _transform_to_wgs84(minx, miny, maxx, maxy, src_srs):
    """Reproject a bbox from ``src_srs`` (osr.SpatialReference) to EPSG:4326."""
    wgs = osr.SpatialReference()
    wgs.ImportFromEPSG(4326)
    # GDAL >= 3 honours authority axis order; force lon/lat so x=lon, y=lat.
    try:
        wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except AttributeError:                                  # pragma: no cover
        pass
    ct = osr.CoordinateTransformation(src_srs, wgs)
    xs, ys = [], []
    for x, y in ((minx, miny), (minx, maxy), (maxx, miny), (maxx, maxy)):
        px, py, _ = ct.TransformPoint(x, y)
        xs.append(px)
        ys.append(py)
    return (min(xs), min(ys), max(xs), max(ys))


def load_aoi(spec, pad_deg=0.0, layer=None):
    """Normalise an AOI specification to a ``(W, S, E, N)`` WGS84 bbox.

    Parameters
    ----------
    spec
        One of:
          * a 4-tuple/list ``(W, S, E, N)`` already in lon/lat degrees;
          * a path to an OGR-readable vector (``.shp``/``.gpkg``/``.geojson``);
            all features' extent is used and reprojected to 4326;
          * a QGIS ``QgsVectorLayer`` (its ``extent()`` + ``crs()`` are used);
          * an object exposing ``.extent()`` and ``.crs()`` (duck-typed QGIS
            layer) -- handled the same way.
    pad_deg
        Degrees to pad the returned bounds on each side (a small pad ~0.02-0.05
        keeps rolling-window rainfall estimators valid to the border).
    layer
        Optional layer name/index for multi-layer OGR sources (GeoPackage).

    Returns
    -------
    tuple
        ``(W, S, E, N)`` padded bounds in EPSG:4326.
    """
    # 1) plain bbox tuple already in lon/lat
    if (isinstance(spec, (tuple, list)) and len(spec) == 4
            and all(isinstance(v, (int, float)) for v in spec)):
        return pad_bounds(tuple(float(v) for v in spec), pad_deg)

    # 2) QGIS vector layer (duck-typed: has extent() + crs())
    if hasattr(spec, "extent") and hasattr(spec, "crs"):
        ext = spec.extent()
        src = osr.SpatialReference()
        src.ImportFromWkt(spec.crs().toWkt())
        bounds = _transform_to_wgs84(ext.xMinimum(), ext.yMinimum(),
                                     ext.xMaximum(), ext.yMaximum(), src)
        return pad_bounds(bounds, pad_deg)

    # 3) path to an OGR vector file
    if isinstance(spec, str) and os.path.exists(spec):
        ds = ogr.Open(spec)
        if ds is None:
            raise ValueError(f"cannot open vector AOI: {spec}")
        lyr = ds.GetLayer(layer) if layer is not None else ds.GetLayer(0)
        minx, maxx, miny, maxy = lyr.GetExtent()       # OGR order: xmin,xmax,ymin,ymax
        src = lyr.GetSpatialRef()
        if src is None:
            raise ValueError(f"{spec} has no CRS; cannot place it on Earth.")
        bounds = _transform_to_wgs84(minx, miny, maxx, maxy, src)
        ds = None
        return pad_bounds(bounds, pad_deg)

    raise TypeError(
        "aoi must be a (W,S,E,N) lon/lat tuple, a QGIS vector layer, or a path "
        f"to a vector file; got {spec!r}")


def bounds_center(bounds):
    """Centre ``(lon, lat)`` of a ``(W, S, E, N)`` bbox."""
    w, s, e, n = bounds
    return ((w + e) / 2.0, (s + n) / 2.0)
