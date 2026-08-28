# -*- coding: utf-8 -*-
"""Shared third-party imports and small compatibility shims.

Centralises GDAL and (optional) QGIS imports so the rest of the package can do
``from ._compat import _HAS_QGIS, QgsGeometry, ...``.  When QGIS is unavailable
(e.g. running the numeric core outside a QGIS environment) the QGIS names are
defined as ``None`` so imports still succeed; only functions that actually touch
QGIS will fail, and only when they are called.
"""
from __future__ import annotations

import numpy as np

# numpy>=2.0 renamed trapz -> trapezoid; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# --- GDAL (always present in QGIS python) ------------------------------------
from osgeo import gdal, ogr, osr
gdal.UseExceptions()

# --- QGIS API (optional) -----------------------------------------------------
_QGIS_NAMES = [
    "QgsProject", "QgsVectorLayer", "QgsRasterLayer", "QgsFeature", "QgsFields",
    "QgsField", "QgsGeometry", "QgsPointXY", "QgsPoint", "QgsVectorFileWriter",
    "QgsCoordinateReferenceSystem", "QgsCoordinateTransform", "QgsWkbTypes",
    "QgsFeatureRequest", "QgsFillSymbol", "QgsLineSymbol", "QgsRendererRange",
    "QgsGraduatedSymbolRenderer", "QgsSymbol", "QgsMarkerSymbol",
]

try:
    from qgis.core import (
        QgsProject, QgsVectorLayer, QgsRasterLayer, QgsFeature, QgsFields,
        QgsField, QgsGeometry, QgsPointXY, QgsPoint, QgsVectorFileWriter,
        QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsWkbTypes,
        QgsFeatureRequest, QgsFillSymbol, QgsLineSymbol, QgsRendererRange,
        QgsGraduatedSymbolRenderer, QgsSymbol, QgsMarkerSymbol,
    )
    from qgis.PyQt.QtGui import QColor
    _HAS_QGIS = True
except Exception:  # pragma: no cover - lets the numeric core import alone
    _HAS_QGIS = False
    QColor = None
    for _n in _QGIS_NAMES:
        globals()[_n] = None
