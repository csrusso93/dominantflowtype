r"""Smoke-test presets — a fast, tiny end-to-end pass to validate the pipeline.

Coarse everything (small AOI clip + coarse subsample/resolution) so the full
cloud → M3C2 → Q\* → IDFVA-bridge chain runs in minutes on a laptop. Use it to
confirm the M3C2 sign convention, CRS handling and the bridge before committing to
a full-resolution run.

    import dominantflowtype as dft
    aoi = dft.make_test_aoi((738405.27, 4387968.70, 739533.81, 4388812.37),
                            "smoketest_aoi.gpkg", epsg=6339, size_m=200)
    dft.run(cfg=dft.smoke_config())      # at [PC-5], pick smoketest_aoi.gpkg
"""
from __future__ import annotations

import os

from ._compat import ogr, osr


def smoke_config(work_epsg=6339, **overrides):
    """Return a speed-tuned :class:`Config` for a first end-to-end pass.

    Coarse DTM/subsample and no denoise so the clipped subset flies. Pass keyword
    ``overrides`` to tweak any field (e.g. ``dtm_resolution=1.0``).
    """
    from .config import Config
    cfg = Config(
        work_epsg=work_epsg,       # normalise compound CRS to a clean code
        dtm_resolution=2.0,        # coarse DTM (fast)
        dtm_method="tin",          # gap-free
        denoise_sfm=False,         # skip on the small clip for speed
        classify_sfm_ground=True,  # still need bare earth for the DTM
        m3c2_subsample=1.0,        # coarse M3C2 (fast)
        m3c2_normal_scale=2.0,
        m3c2_cyl_radius=1.0,
        pointcloud_dry_run=False,  # actually run — it's meant to be quick
    )
    for k, v in overrides.items():
        if not hasattr(cfg, k):
            raise AttributeError(f"Config has no field {k!r}")
        setattr(cfg, k, v)
    return cfg


def make_test_aoi(bounds, out_path, epsg, size_m=200.0):
    """Write a small square AOI polygon centred in ``bounds`` (for `[PC-5]` clip).

    Parameters
    ----------
    bounds : (minx, miny, maxx, maxy)
        Extent of the cloud/DEM, in EPSG:``epsg`` map units (metres).
    out_path : str
        Output GeoPackage (single polygon feature).
    epsg : int
        CRS of ``bounds`` / the output (e.g. 6339).
    size_m : float
        Side length of the square AOI [m] (default 200).
    """
    minx, miny, maxx, maxy = bounds
    cx, cy = 0.5 * (minx + maxx), 0.5 * (miny + maxy)
    h = size_m / 2.0
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in [(cx - h, cy - h), (cx + h, cy - h), (cx + h, cy + h),
                 (cx - h, cy + h), (cx - h, cy - h)]:
        ring.AddPoint_2D(x, y)     # 2D — avoid spurious Z in the overlay polygon
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    drv = ogr.GetDriverByName("GPKG")
    if os.path.exists(out_path):
        drv.DeleteDataSource(out_path)
    ds = drv.CreateDataSource(out_path)
    lyr = ds.CreateLayer("aoi", srs, ogr.wkbPolygon)
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometry(poly)
    lyr.CreateFeature(feat)
    feat = None
    ds = None
    print(f"  [smoketest] wrote {size_m:.0f} m AOI -> {out_path} (EPSG:{epsg}, "
          f"centre {cx:.1f}, {cy:.1f})")
    return out_path
