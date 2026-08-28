r"""Stage 0 acquisition — fetch USGS 3DEP LiDAR point clouds for an AOI.

USGS 3DEP LiDAR is published as **Entwine Point Tiles (EPT)** / COPC in the public
AWS bucket ``usgs-lidar-public`` (one ``ept.json`` per project). This module:

* ``list_3dep_projects(aoi_bounds_4326)`` — query the public boundary index and
  print the project(s) whose footprint covers the AOI, with their EPT URLs.
* ``acquire_3dep_cloud(ept_url, watershed_path, out_path, cfg)`` — clip an EPT to
  the watershed (QGIS ``pdal:clip``) to a local LAZ/COPC you can feed into the
  point-cloud pathway as the *pre-event* cloud.

Typical use: pick your project on the USGS 3DEP LidarExplorer
(https://apps.nationalmap.gov/lidar-explorer/), copy its EPT/COPC URL, then clip.
``list_3dep_projects`` helps find candidates programmatically.

References
----------
USGS 3DEP LiDAR (usgs-lidar-public EPT on AWS); Entwine/EPT (entwine.io);
QGIS PDAL provider ``pdal:clip``.
"""
from __future__ import annotations

import json
import os
import urllib.request

from ._compat import osr

# Public index of 3DEP EPT resources (project footprints + names).
_RESOURCES_URL = ("https://raw.githubusercontent.com/hobu/usgs-lidar/master/"
                  "boundaries/resources.geojson")
_EPT_BASE = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public"


def _aoi_bounds_4326(watershed_path):
    """Return (minx, miny, maxx, maxy) of a vector layer in EPSG:4326."""
    from ._compat import ogr
    ds = ogr.Open(watershed_path)
    if ds is None:
        raise IOError(f"cannot open AOI vector: {watershed_path}")
    lyr = ds.GetLayer()
    src_srs = lyr.GetSpatialRef()
    x0, x1, y0, y1 = lyr.GetExtent()          # (minx, maxx, miny, maxy)
    dst = osr.SpatialReference(); dst.ImportFromEPSG(4326)
    if src_srs is not None:
        try:
            src_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except Exception:
            pass
        ct = osr.CoordinateTransformation(src_srs, dst)
        corners = [ct.TransformPoint(x, y)[:2]
                   for x in (x0, x1) for y in (y0, y1)]
        xs = [c[0] for c in corners]; ys = [c[1] for c in corners]
        return min(xs), min(ys), max(xs), max(ys)
    return x0, y0, x1, y1


def _bbox_overlaps(a, b):
    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])


def list_3dep_projects(watershed_path, max_show=25, timeout=60):
    """Print 3DEP EPT projects whose footprint bbox overlaps the AOI.

    Returns a list of (name, ept_url). Needs internet. Footprint index is served
    in EPSG:4326; the AOI is transformed to 4326 for the test.
    """
    aoi = _aoi_bounds_4326(watershed_path)
    print(f"  [acquire] AOI bbox (EPSG:4326): {tuple(round(v, 5) for v in aoi)}")
    print(f"  [acquire] fetching 3DEP resource index …")
    with urllib.request.urlopen(_RESOURCES_URL, timeout=timeout) as resp:
        gj = json.load(resp)
    hits = []
    for feat in gj.get("features", []):
        geom = feat.get("geometry") or {}
        coords = geom.get("coordinates")
        if not coords:
            continue
        # flatten polygon/multipolygon rings to a bbox
        pts = _flatten_coords(coords)
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        fbbox = (min(xs), min(ys), max(xs), max(ys))
        if _bbox_overlaps(aoi, fbbox):
            props = feat.get("properties", {})
            name = props.get("name") or props.get("Name") or "?"
            url = props.get("url") or f"{_EPT_BASE}/{name}/ept.json"
            hits.append((name, url))
    if not hits:
        print("  [acquire] no 3DEP project footprints overlap the AOI (or the "
              "index moved). Use the USGS 3DEP LidarExplorer to find a URL.")
    else:
        print(f"  [acquire] {len(hits)} candidate project(s):")
        for name, url in hits[:max_show]:
            print(f"        {name}\n            {url}")
    return hits


def _flatten_coords(coords):
    """Recursively flatten GeoJSON coordinate nesting to a list of (x, y)."""
    out = []
    if (isinstance(coords, (list, tuple)) and len(coords) >= 2
            and all(isinstance(v, (int, float)) for v in coords[:2])):
        return [(coords[0], coords[1])]
    for c in coords:
        out.extend(_flatten_coords(c))
    return out


def acquire_3dep_cloud(ept_url, watershed_path, out_path, cfg):
    """Clip a 3DEP EPT/COPC to the watershed via QGIS ``pdal:clip``.

    Returns the local output path, or ``None`` on dry-run. ``ept_url`` may be a
    remote ``ept.json``/COPC URL or a local file. Requires QGIS Processing.
    """
    from .pointcloud import _processing, _resolve, _run
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "clip_alg_id", "") or "pdal:clip",
                      [["clip"], ["pdal", "clip"]])
    params = {
        "INPUT": ept_url,
        "OVERLAY": watershed_path,
        "OUTPUT": out_path,
    }
    params.update(getattr(cfg, "acquire_params", {}) or {})
    res = _run(processing, alg_id, params,
               getattr(cfg, "pointcloud_dry_run", False),
               f"clip 3DEP EPT to watershed ({ept_url})")
    return out_path if res is not None else None
