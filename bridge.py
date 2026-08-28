# -*- coding: utf-8 -*-
r"""Bridge to Guido's IDFVA (Incremental Debris Flow Volume Analyzer).

`dominantflowtype` computes Q\* inside QGIS; **volume is delegated to IDFVA**,
which runs in its own Python venv (`pip install idfva`). We keep the packages
separate and hand IDFVA a self-contained *bundle* of files it can consume:

    bridge_bundle/
        flow_path.shp   the thalweg (IDFVA's "flow path of interest")
        post_dtm.tif    bare-earth DEM for hydrology (process_watersheds.py)
        change.tif      signed change raster (deposition +, erosion -)
        manifest.json   CRS/units/provenance for reproducibility

IDFVA hardcodes EPSG:26913 and reads ESRI Shapefiles (see INTEGRATION_PLAN.md);
the bundle therefore emits a `.shp` flow path and records the true EPSG in the
manifest so the IDFVA-side wrapper can reproject/parameterise correctly.
"""
from __future__ import annotations

import json
import os
import shutil

from ._compat import _HAS_QGIS, QgsVectorFileWriter


def export_bridge(bundle_dir, thalweg_layer, post_dtm_path, change_raster_path,
                  work_crs, meta=None, epsg_hint=None):
    r"""Write the IDFVA bridge bundle. Returns the bundle directory.

    Parameters
    ----------
    bundle_dir : str
        Output folder (created if needed); the bundle is written here.
    thalweg_layer : QgsVectorLayer
        Thalweg line layer, already in ``work_crs`` — becomes ``flow_path.shp``.
    post_dtm_path : str
        Post-event bare-earth DEM (from Stage 0 or a supplied DEM).
    change_raster_path : str or None
        Signed change raster; ``None`` for the post-only (Q\*-only) case, in which
        case no volume bundle is written and the function returns ``None``.
    work_crs : QgsCoordinateReferenceSystem
        CRS of all bundle layers (projected, metres).
    meta : dict, optional
        Extra provenance merged into ``manifest.json``.
    """
    if not _HAS_QGIS:
        print("  [bridge] QGIS unavailable; skipping IDFVA bundle export.")
        return None
    if not change_raster_path:
        print("  [bridge] no change raster (post-only run) → IDFVA volume bundle "
              "not applicable; Q* outputs only.")
        return None

    os.makedirs(bundle_dir, exist_ok=True)

    # 1) flow path -> ESRI Shapefile in the working CRS
    flow_shp = os.path.join(bundle_dir, "flow_path.shp")
    err = QgsVectorFileWriter.writeAsVectorFormat(
        thalweg_layer, flow_shp, "UTF-8", work_crs, "ESRI Shapefile")
    # writeAsVectorFormat returns (errCode, msg) on some builds; 0 == NoError
    if isinstance(err, (tuple, list)) and err and err[0] != 0:
        print(f"  [bridge] WARNING writing flow_path.shp: {err}")

    # 2) rasters -> copy into the bundle (self-contained)
    dtm_out = os.path.join(bundle_dir, "post_dtm.tif")
    change_out = os.path.join(bundle_dir, "change.tif")
    _copy_raster(post_dtm_path, dtm_out)
    _copy_raster(change_raster_path, change_out)

    # 3) manifest
    epsg = None
    try:
        epsg = work_crs.postgisSrid() or int(work_crs.authid().split(":")[1])
    except Exception:
        pass
    # Compound / User-Defined CRS (e.g. UTM10N+NAVD88 from Metashape) often has a
    # blank authid; fall back to the caller's work_epsg hint so the manifest still
    # carries a usable EPSG for the IDFVA side.
    if not epsg and epsg_hint:
        epsg = int(epsg_hint)
        print(f"  [bridge] CRS authid blank; using epsg hint {epsg} "
              f"(set Config.work_epsg to reproject clouds to a clean code).")
    manifest = {
        "producer": "dominantflowtype",
        "consumer": "idfva (separate venv)",
        "crs_authid": work_crs.authid() if work_crs else None,
        "epsg": epsg,
        "units": "metres",
        "change_convention": "post-minus-pre; positive=deposition, negative=erosion",
        "files": {
            "flow_path": os.path.basename(flow_shp) if os.path.exists(flow_shp) else None,
            "post_dtm": os.path.basename(dtm_out) if os.path.exists(dtm_out) else None,
            "change_raster": os.path.basename(change_out) if os.path.exists(change_out) else None,
        },
        "idfva_notes": [
            "IDFVA hardcodes EPSG:26913 in las2ras.py / generate_query_points.py / "
            "prepare_path.py — reproject to 'epsg' above or patch those constants.",
            "IDFVA reads ESRI Shapefiles via fiona/ogr; flow_path.shp is provided.",
        ],
    }
    if meta:
        manifest.update(meta)
    with open(os.path.join(bundle_dir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"  [bridge] IDFVA bundle -> {bundle_dir}")
    print(f"           flow_path.shp | post_dtm.tif | change.tif | manifest.json "
          f"(EPSG {epsg})")
    return bundle_dir


def _copy_raster(src, dst):
    """Copy a raster and its sidecars (.tfw/.aux.xml/.prj) into the bundle."""
    if not src or not os.path.exists(src):
        print(f"  [bridge] WARNING: raster not found, not bundled: {src}")
        return
    if os.path.abspath(src) == os.path.abspath(dst):
        return
    shutil.copy2(src, dst)
    base = os.path.splitext(src)[0]
    for ext in (".tfw", ".aux.xml", ".tif.aux.xml", ".prj", ".ovr"):
        side = base + ext if not ext.startswith(".tif") else src + ext[4:]
        if os.path.exists(side):
            try:
                shutil.copy2(side, os.path.join(os.path.dirname(dst),
                                                os.path.basename(side)))
            except Exception:
                pass
