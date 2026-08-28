# -*- coding: utf-8 -*-
"""Workbook (xlsx/CSV) output and styled QGIS layer construction."""
from __future__ import annotations

import os

import numpy as np

from ._compat import (
    _HAS_QGIS, QColor, QgsVectorLayer, QgsFeature, QgsField, QgsGeometry,
    QgsFillSymbol, QgsLineSymbol, QgsRendererRange, QgsGraduatedSymbolRenderer,
    QgsProject, QgsVectorFileWriter,
)


# Column layout of the output sheets.
PRE_COLUMNS = [
    "basin", "transect_id", "latitude", "longitude", "station_m",
    "slope", "depth_m", "width_m", "area_xs_m2", "wetted_perimeter_m",
    "hydraulic_radius_m", "area_upstream_m2", "velocity_ms", "Q_peak_m3s",
    "Q_fluv_m3s", "Qstar", "flow_type", "continuity", "usable", "flag_reason",
]
POST_EXTRA = ["dz_bed_m", "volume_m3"]     # appended for the post-event sheet


def write_workbook(rows_pre, rows_post, cfg, out_xlsx, rain_meta):
    """Write pre/post/usable-only sheets to xlsx (+ CSV fallbacks)."""
    import pandas as pd
    sheets = {}
    if rows_pre:
        sheets["pre_event"] = pd.DataFrame(rows_pre, columns=PRE_COLUMNS)
    if rows_post:
        cols = PRE_COLUMNS[:]
        # insert post-only columns before flags
        insert_at = cols.index("usable")
        cols = cols[:insert_at] + POST_EXTRA + cols[insert_at:]
        sheets["post_event"] = pd.DataFrame(rows_post, columns=cols)
        usable_df = sheets["post_event"][sheets["post_event"]["usable"] == True]
        sheets["usable_only"] = usable_df.reset_index(drop=True)
    elif rows_pre:
        usable_df = sheets["pre_event"][sheets["pre_event"]["usable"] == True]
        sheets["usable_only"] = usable_df.reset_index(drop=True)

    # rainfall metadata sheet
    sheets["rainfall_I30"] = pd.DataFrame([rain_meta])

    try:
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
            for name, df in sheets.items():
                df.to_excel(xw, sheet_name=name[:31], index=False)
        print(f"  [output] workbook -> {out_xlsx}")
    except Exception as ex:
        # CSV fallback only if the single-workbook write fails
        print(f"  [output] xlsx failed ({ex}); writing CSV fallbacks instead")
        for name, df in sheets.items():
            df.to_csv(os.path.splitext(out_xlsx)[0] + f"_{name}.csv", index=False)
    return sheets


# ---- QGIS layer construction & diverging styling ----------------------------
def _diverging_ramp_color(value, center, vmin, vmax, low, mid, high):
    """Interpolate a two-sided diverging color (QColor)."""
    def lerp(a, b, t):
        return QColor(int(a.red() + (b.red() - a.red()) * t),
                      int(a.green() + (b.green() - a.green()) * t),
                      int(a.blue() + (b.blue() - a.blue()) * t))
    if not np.isfinite(value):
        return QColor(180, 180, 180)
    if value <= center:
        denom = (center - vmin) or 1.0
        t = max(0.0, min(1.0, (value - vmin) / denom))
        return lerp(low, mid, t)
    denom = (vmax - center) or 1.0
    t = max(0.0, min(1.0, (value - center) / denom))
    return lerp(mid, high, t)


def _graduated_diverging(layer, field, center, low, mid, high,
                         n_classes=8, geom="polygon", unit=""):
    """Build a graduated renderer with a diverging ramp centred on `center`.

    The class edges are made **symmetric about `center`** so that the two sides
    always map to opposite ends of the ramp regardless of how lopsided the data
    are (e.g. max erosion -24 m vs max deposition +23 m). Without this, an
    asymmetric min/max lets every class land on one side of the ramp and the
    layer renders as a single-hue gradient instead of a true diverging one.
    """
    vals = [f[field] for f in layer.getFeatures()
            if f[field] is not None and np.isfinite(f[field])]
    if not vals:
        return
    vmin, vmax = min(vals), max(vals)
    # symmetric half-span around the center -> balanced, always-diverging classes
    span = max(abs(vmin - center), abs(vmax - center)) or 1.0
    below = np.linspace(center - span, center, n_classes // 2 + 1)
    above = np.linspace(center, center + span, n_classes // 2 + 1)
    edges = list(dict.fromkeys(list(below) + list(above)))
    suffix = f" {unit}" if unit else ""
    ranges = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        midv = 0.5 * (lo + hi)
        col = _diverging_ramp_color(midv, center, center - span, center + span,
                                    low, mid, high)
        if geom == "polygon":
            sym = QgsFillSymbol.createSimple(
                {"color": col.name(), "outline_color": "0,0,0,60",
                 "outline_width": "0.1"})
        else:
            sym = QgsLineSymbol.createSimple({"color": col.name(), "width": "0.8"})
        # signed labels so erosion (negative) vs deposition (positive) is explicit
        ranges.append(QgsRendererRange(lo, hi, sym,
                                       f"{lo:+.3g} – {hi:+.3g}{suffix}"))
    renderer = QgsGraduatedSymbolRenderer(field, ranges)
    layer.setRenderer(renderer)
    layer.triggerRepaint()


def _write_styled_gpkg(mem_layer, gpkg_path, layer_name):
    """Write a memory layer to GeoPackage AND persist its QGIS style.

    ``QgsVectorFileWriter`` copies only geometry + attributes, not the renderer,
    so a freshly opened .gpkg would lose the diverging palette and fall back to a
    default (single-hue) ramp. We store the style two ways so it survives:
      * a ``.qml`` sidecar (auto-loaded by QGIS for a single-layer file), and
      * embedded in the GeoPackage's ``layer_styles`` table as the default style.
    """
    QgsVectorFileWriter.writeAsVectorFormat(
        mem_layer, gpkg_path, "UTF-8", mem_layer.crs(), "GPKG")
    # .qml sidecar next to the file
    mem_layer.saveNamedStyle(os.path.splitext(gpkg_path)[0] + ".qml")
    # embed as the default style inside the GeoPackage itself
    try:
        uri = f"{gpkg_path}|layername={layer_name}"
        saved = QgsVectorLayer(uri, layer_name, "ogr")
        if saved.isValid():
            saved.setRenderer(mem_layer.renderer().clone())
            saved.saveStyleToDatabase(layer_name, "dominantflowtype diverging style",
                                      True, "")
    except Exception as ex:
        print(f"  [output] could not embed style in {os.path.basename(gpkg_path)}: {ex}")


def _make_memory_layer(name, geom_type, crs, fields_spec, features):
    """Create an in-memory QgsVectorLayer and populate it."""
    uri = f"{geom_type}?crs={crs.authid() or crs.toWkt()}"
    lyr = QgsVectorLayer(uri, name, "memory")
    pr = lyr.dataProvider()
    qfields = []
    for fname, ftype in fields_spec:
        qfields.append(QgsField(fname, ftype))
    pr.addAttributes(qfields)
    lyr.updateFields()
    feats = []
    for geom, attrs in features:
        feat = QgsFeature(lyr.fields())
        feat.setGeometry(geom)
        feat.setAttributes(attrs)
        feats.append(feat)
    pr.addFeatures(feats)
    lyr.updateExtents()
    return lyr


def add_qgis_layers(records, reaches, work_crs, cfg, out_dir, have_dod):
    """Add transect-lines, inundation polygons (Q*), and erosion-deposition
    polygons to the current QGIS project, styled with diverging palettes."""
    if not _HAS_QGIS:
        return
    from qgis.PyQt.QtCore import QVariant

    # ---- inundation polygons coloured by Q* (green<1, white=1, purple>1) ----
    poly_fields = [("transect_id", QVariant.Int), ("Qstar", QVariant.Double),
                   ("flow_type", QVariant.String), ("area_xs_m2", QVariant.Double),
                   ("depth_m", QVariant.Double), ("width_m", QVariant.Double),
                   ("usable", QVariant.String)]
    poly_feats = []
    for r in records:
        if r.get("polygon") is None:
            continue
        poly_feats.append((r["polygon"],
                           [r["transect_id"], _num(r.get("Qstar_post", r.get("Qstar"))),
                            r.get("flow_type_post", r.get("flow_type")),
                            _num(r.get("area_xs_m2_post", r.get("area_xs_m2"))),
                            _num(r.get("depth_m_post", r.get("depth_m"))),
                            _num(r.get("width_m_post", r.get("width_m"))),
                            "yes" if r.get("usable_post", r.get("usable")) else "no"]))
    if poly_feats:
        plyr = _make_memory_layer("DFT_inundation_Qstar", "Polygon", work_crs,
                                  poly_fields, poly_feats)
        _graduated_diverging(plyr, "Qstar", cfg.qstar_threshold,
                             QColor(0, 136, 55),    # green  (flood, Q*<1)
                             QColor(255, 255, 255),  # white  (Q*=1)
                             QColor(123, 50, 148),   # purple (debris flow, Q*>1)
                             geom="polygon", unit="")   # Q* is dimensionless
        QgsProject.instance().addMapLayer(plyr)
        _write_styled_gpkg(plyr, os.path.join(out_dir, "DFT_inundation_Qstar.gpkg"),
                           "DFT_inundation_Qstar")

    # ---- transect centrelines (reference) -----------------------------------
    line_fields = [("transect_id", QVariant.Int), ("Qstar", QVariant.Double)]
    line_feats = [(r["transect_geom"],
                   [r["transect_id"], _num(r.get("Qstar_post", r.get("Qstar")))])
                  for r in records if r.get("transect_geom") is not None]
    if line_feats:
        llyr = _make_memory_layer("DFT_transects", "LineString", work_crs,
                                  line_fields, line_feats)
        QgsProject.instance().addMapLayer(llyr)

    # ---- erosion / deposition polygons (only if DoD available) --------------
    if have_dod and reaches:
        ed_fields = [("reach_id", QVariant.Int), ("mean_dz_m", QVariant.Double),
                     ("area_m2", QVariant.Double), ("volume_m3", QVariant.Double)]
        ed_feats = [(rc["polygon"], [rc["reach_id"], _num(rc["mean_dz"]),
                                     _num(rc["area"]), _num(rc["volume"])])
                    for rc in reaches if rc.get("polygon") is not None]
        if ed_feats:
            elyr = _make_memory_layer("DFT_erosion_deposition", "Polygon",
                                      work_crs, ed_fields, ed_feats)
            _graduated_diverging(elyr, "mean_dz_m", 0.0,
                                 QColor(202, 0, 32),      # red  (erosion, dz<0)
                                 QColor(247, 247, 247),   # white (no change)
                                 QColor(5, 113, 176),     # blue (deposition, dz>0)
                                 geom="polygon", unit="m")   # bed change in metres
            QgsProject.instance().addMapLayer(elyr)
            _write_styled_gpkg(elyr,
                               os.path.join(out_dir, "DFT_erosion_deposition.gpkg"),
                               "DFT_erosion_deposition")


def _num(v):
    return float(v) if v is not None and np.isfinite(v) else None
