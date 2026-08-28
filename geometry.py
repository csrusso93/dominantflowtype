# -*- coding: utf-8 -*-
"""Vector/CRS helpers, transect casting and trimline-to-transect intersection."""
from __future__ import annotations

import os
import math

from ._compat import (
    QgsGeometry, QgsPointXY, QgsWkbTypes, QgsVectorLayer,
    QgsCoordinateTransform, QgsProject,
)


# ---------------------------------------------------------------------------
# vector / CRS helpers
# ---------------------------------------------------------------------------
def _load_vector(path):
    """Load a .shp or first layer of a .gpkg as a QgsVectorLayer."""
    if path.lower().endswith(".gpkg"):
        # take the first (or only) layer
        lyr = QgsVectorLayer(path, os.path.basename(path), "ogr")
        sub = lyr.dataProvider().subLayers()
        if sub:
            name = sub[0].split("!!::!!")[1]
            lyr = QgsVectorLayer(f"{path}|layername={name}", name, "ogr")
    else:
        lyr = QgsVectorLayer(path, os.path.basename(path), "ogr")
    if not lyr.isValid():
        raise IOError(f"invalid vector layer: {path}")
    return lyr


def _transform_geom(geom, src_crs, dst_crs):
    """Return a copy of geom transformed src->dst (no-op if equal)."""
    if src_crs is None or dst_crs is None or src_crs == dst_crs:
        return QgsGeometry(geom)
    xform = QgsCoordinateTransform(src_crs, dst_crs, QgsProject.instance())
    g = QgsGeometry(geom)
    g.transform(xform)
    return g


def _merged_lines(layer, dst_crs):
    """Return list of QgsGeometry (single LineStrings) in dst_crs."""
    src = layer.crs()
    out = []
    for f in layer.getFeatures():
        g = _transform_geom(f.geometry(), src, dst_crs)
        if g.isMultipart():
            for part in g.asMultiPolyline():
                if len(part) >= 2:
                    out.append(QgsGeometry.fromPolylineXY(part))
        else:
            pl = g.asPolyline()
            if len(pl) >= 2:
                out.append(QgsGeometry.fromPolylineXY(pl))
    return out


def _feature_flag(feat, cfg):
    """Return True if a feature is a MAPPED/certain trimline (or no flag field)."""
    if cfg.trimline_flag_field not in [fld.name() for fld in feat.fields()]:
        return True
    val = feat[cfg.trimline_flag_field]
    if val is None:
        return True
    return str(val).strip().lower() in cfg.trimline_certain_values


def _longest(geoms):
    if not geoms:
        raise ValueError("no line geometry found")
    return max(geoms, key=lambda g: g.length())


# ---------------------------------------------------------------------------
# transect casting & bank classification
# ---------------------------------------------------------------------------
def _densify_line(geom, spacing):
    """Return list of (x, y, s) stations every `spacing` m along a line geom."""
    length = geom.length()
    stations = []
    d = 0.0
    while d <= length + 1e-6:
        p = geom.interpolate(d)
        if p and not p.isEmpty():
            pt = p.asPoint()
            stations.append((pt.x(), pt.y(), d))
        d += spacing
    return stations, length


def cast_transects(thalweg_geom, cfg):
    """Cast transects perpendicular to the thalweg every `transect_spacing` m.

    Returns list of dicts with station point, unit tangent, and transect
    endpoints (left_end, right_end) at +/- transect_halfwidth.
    """
    stations, length = _densify_line(thalweg_geom, cfg.transect_spacing)
    out = []
    for (x, y, s) in stations:
        # tangent via finite difference along the line
        s0 = max(0.0, s - cfg.tangent_delta)
        s1 = min(length, s + cfg.tangent_delta)
        p0 = thalweg_geom.interpolate(s0).asPoint()
        p1 = thalweg_geom.interpolate(s1).asPoint()
        tx, ty = (p1.x() - p0.x()), (p1.y() - p0.y())
        tn = math.hypot(tx, ty)
        if tn < 1e-9:
            continue
        tx, ty = tx / tn, ty / tn
        # left-normal = (-ty, tx); right = (ty, -tx)
        lx, ly = -ty, tx
        L = cfg.transect_halfwidth
        left_end = (x + lx * L, y + ly * L)
        right_end = (x - lx * L, y - ly * L)
        out.append(dict(idx=len(out), s=s, x=x, y=y,
                        tx=tx, ty=ty, lx=lx, ly=ly,
                        left_end=left_end, right_end=right_end))
    return out


def _side_of(station, tx, ty, px, py):
    """Sign of cross(tangent, point-station): +1 left, -1 right."""
    dx, dy = px - station[0], py - station[1]
    cross = tx * dy - ty * dx
    return 1 if cross >= 0 else -1


def trimline_hits(transect, trimlines, cfg, post_sampler):
    """Intersect one transect with all MAPPED trimlines; classify L/R.

    Returns dict {'L': (station_dist, x, y, z) or None, 'R': ...} where
    station_dist is signed distance from the thalweg point along the transect
    (positive toward the left endpoint).
    """
    st = (transect["x"], transect["y"])
    line = QgsGeometry.fromPolylineXY(
        [QgsPointXY(*transect["left_end"]), QgsPointXY(*transect["right_end"])]
    )
    best = {"L": None, "R": None}
    for tg in trimlines:
        inter = line.intersection(tg)
        if inter is None or inter.isEmpty():
            continue
        pts = []
        if inter.wkbType() in (QgsWkbTypes.Point, QgsWkbTypes.PointZ):
            pts = [inter.asPoint()]
        elif inter.wkbType() in (QgsWkbTypes.MultiPoint, QgsWkbTypes.MultiPointZ):
            pts = inter.asMultiPoint()
        else:
            # line overlap: take midpoint
            c = inter.centroid()
            if c and not c.isEmpty():
                pts = [c.asPoint()]
        for p in pts:
            side = _side_of(st, transect["tx"], transect["ty"], p.x(), p.y())
            key = "L" if side > 0 else "R"
            dist = math.hypot(p.x() - st[0], p.y() - st[1]) * side
            z = post_sampler.sample(p.x(), p.y()) if post_sampler else float("nan")
            cur = best[key]
            # keep the closest crossing to the channel on each side
            if cur is None or abs(dist) < abs(cur[0]):
                best[key] = (dist, p.x(), p.y(), z)
    return best
