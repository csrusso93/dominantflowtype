# -*- coding: utf-8 -*-
"""Upstream drainage area (A_us) via portable D8 flow accumulation.

A priority-flood depression fill (Barnes et al., 2014) followed by D8 steepest-
descent routing and topological accumulation, computed on a coarse, watershed-
clipped resample of the DEM.  numpy-only (no GRASS/SAGA dependency).
"""
from __future__ import annotations

import math

import numpy as np

from ._compat import gdal, ogr, osr, QgsCoordinateReferenceSystem
from .geometry import _transform_geom


def _priority_flood_fill(dem):
    r"""Barnes (2014) priority-flood depression filling, **with epsilon**.

    Cells are raised to just *above* the spill elevation, not exactly to it, so
    a filled depression drains monotonically back out toward its spill point.

    This matters. :func:`_d8_accumulation` assigns a receiver only where the
    steepest slope is ``> 0``. Filling a pit to exactly the spill level leaves a
    perfectly flat region in which no cell has a downslope neighbour, so every
    cell in it is a terminal sink and **all upstream flow stops there**.

    Measured on a 0.445 km2 burned catchment: filling to exactly the spill level
    left 546 terminal flat cells and A_us reached only 0.257 km2, i.e. 57.9 % of
    the basin. With the epsilon it is 0.448 km2 (100.7 %) -- a 1.74x correction.
    Since ``Q_fluv = A_us * I``, Q\* there had been over-predicted by 74 %.
    Catchments without significant depressions are unaffected (two others moved
    100.4 -> 100.8 % and 100.5 -> 100.6 %).

    NaN = nodata.
    """
    import heapq
    ny, nx = dem.shape
    filled = dem.copy()
    closed = ~np.isfinite(dem)
    pq = []
    seen = closed.copy()
    # seed with all valid border cells
    for i in range(ny):
        for j in (0, nx - 1):
            if not seen[i, j]:
                heapq.heappush(pq, (dem[i, j], i, j)); seen[i, j] = True
    for j in range(nx):
        for i in (0, ny - 1):
            if not seen[i, j]:
                heapq.heappush(pq, (dem[i, j], i, j)); seen[i, j] = True
    nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    while pq:
        elev, i, j = heapq.heappop(pq)
        for di, dj in nbrs:
            ni, nj = i + di, j + dj
            if 0 <= ni < ny and 0 <= nj < nx and not seen[ni, nj]:
                seen[ni, nj] = True
                ne = dem[ni, nj]
                if ne <= elev:
                    # Raise to just ABOVE the spill level. `<=` (not `<`) also
                    # catches cells already exactly at it, which would otherwise
                    # form a flat. One ULP is enough for slope > 0 in float64
                    # and is ~1e-13 m at these elevations, i.e. far below any
                    # DEM's real precision.
                    ne = np.nextafter(elev, np.inf)
                filled[ni, nj] = ne
                heapq.heappush(pq, (ne, ni, nj))
    return filled


def _d8_accumulation(filled, cellsize, valid):
    """D8 flow accumulation (cell counts) on a filled DEM. Returns accum array."""
    ny, nx = filled.shape
    nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    dist = np.array([1, 1, 1, 1, math.sqrt(2), math.sqrt(2),
                     math.sqrt(2), math.sqrt(2)]) * cellsize
    # steepest-descent receiver for every cell
    rec_i = np.full((ny, nx), -1, dtype=np.int64)
    rec_j = np.full((ny, nx), -1, dtype=np.int64)
    z = np.where(valid, filled, np.inf)
    for k, (di, dj) in enumerate(nbrs):
        zi0 = z[1:-1, 1:-1]
        zin = z[1 + di:ny - 1 + di, 1 + dj:nx - 1 + dj]
        with np.errstate(invalid="ignore"):
            slope = (zi0 - zin) / dist[k]
        if k == 0:
            best = np.full_like(zi0, -np.inf)
            bi = np.full(zi0.shape, -1, dtype=np.int64)
            bj = np.full(zi0.shape, -1, dtype=np.int64)
        upd = slope > best
        best = np.where(upd, slope, best)
        bi = np.where(upd, di, bi)
        bj = np.where(upd, dj, bj)
    ii, jj = np.mgrid[1:ny - 1, 1:nx - 1]
    has = best > 0
    rec_i[1:-1, 1:-1] = np.where(has, ii + bi, -1)
    rec_j[1:-1, 1:-1] = np.where(has, jj + bj, -1)

    # accumulate in order of descending elevation (valid DAG for D8)
    accum = np.where(valid, 1.0, 0.0)
    order = np.argsort(-np.where(valid, filled, -np.inf), axis=None)
    ri = rec_i.ravel(); rj = rec_j.ravel()
    acc = accum.ravel(); vf = valid.ravel()
    for idx in order:
        if not vf[idx]:
            continue
        di, dj = ri[idx], rj[idx]
        if di < 0:
            continue
        acc[di * nx + dj] += acc[idx]
    return acc.reshape(ny, nx)


class FlowRouter:
    """Compute A_us at arbitrary points by D8 accumulation on a resampled,
    watershed-clipped DEM. Portable (numpy only); optional GRASS backend."""

    def __init__(self, dem_path, watershed_layer, cfg, work_crs_wkt):
        self.cfg = cfg
        # 1. clip/resample DEM to coarse grid over the buffered watershed bbox
        ext = watershed_layer.extent()
        b = cfg.flow_accum_buffer
        bounds = (ext.xMinimum() - b, ext.yMinimum() - b,
                  ext.xMaximum() + b, ext.yMaximum() + b)
        res = cfg.flow_accum_res
        warped = gdal.Warp(
            "", dem_path, format="MEM",
            outputBounds=bounds, xRes=res, yRes=res, resampleAlg="bilinear",
            dstSRS=work_crs_wkt,
        )
        self.gt = warped.GetGeoTransform()
        self.inv_gt = gdal.InvGeoTransform(self.gt)
        band = warped.GetRasterBand(1)
        dem = band.ReadAsArray().astype("float64")
        nd = band.GetNoDataValue()
        if nd is not None:
            dem[dem == nd] = np.nan
        dem[np.abs(dem) > 1e30] = np.nan
        self.ny, self.nx = dem.shape

        # 2. rasterize watershed as an in-basin mask
        mask = self._rasterize_mask(watershed_layer, warped, work_crs_wkt)
        valid = np.isfinite(dem) & (mask > 0)
        dem = np.where(valid, dem, np.nan)

        # 3. fill + accumulate
        filled = _priority_flood_fill(dem)
        self.valid = np.isfinite(filled)
        self.accum = _d8_accumulation(filled, res, self.valid)
        self.cell_area = res * res

        # 4. warn if the DEM does not fully cover the watershed (A_us undercount)
        n_basin = int((mask > 0).sum())
        n_covered = int((self.valid & (mask > 0)).sum())
        if n_basin > 0:
            cov = n_covered / n_basin
            if cov < 0.98:
                print(f"  [A_us] WARNING: routing DEM covers only "
                      f"{cov*100:.1f}% of the watershed; upstream areas for "
                      f"cross-sections draining the uncovered part will be "
                      f"UNDER-estimated. Consider supplying a full-coverage "
                      f"DEM (e.g. USGS 10 m via OpenTopography) for routing.")

    def _rasterize_mask(self, layer, ref_ds, work_crs_wkt):
        drv = gdal.GetDriverByName("MEM")
        mem = drv.Create("", ref_ds.RasterXSize, ref_ds.RasterYSize, 1,
                         gdal.GDT_Byte)
        mem.SetGeoTransform(ref_ds.GetGeoTransform())
        mem.SetProjection(work_crs_wkt)
        # write features to an in-memory OGR layer
        ogr_drv = ogr.GetDriverByName("Memory")
        ods = ogr_drv.CreateDataSource("m")
        srs = osr.SpatialReference()
        srs.ImportFromWkt(work_crs_wkt)
        olyr = ods.CreateLayer("w", srs, ogr.wkbPolygon)
        src_crs = layer.crs()
        dst_crs = QgsCoordinateReferenceSystem.fromWkt(work_crs_wkt)
        for feat in layer.getFeatures():
            g = _transform_geom(feat.geometry(), src_crs, dst_crs)
            of = ogr.Feature(olyr.GetLayerDefn())
            of.SetGeometry(ogr.CreateGeometryFromWkt(g.asWkt()))
            olyr.CreateFeature(of)
        gdal.RasterizeLayer(mem, [1], olyr, burn_values=[1])
        return mem.GetRasterBand(1).ReadAsArray()

    def upstream_area(self, x, y):
        """Snap (x,y) to the max-accumulation cell within snap_radius; return m^2."""
        px, py = gdal.ApplyGeoTransform(self.inv_gt, x, y)
        pj, pi = int(px), int(py)
        r = int(math.ceil(self.cfg.snap_radius / abs(self.gt[1])))
        best = 0.0
        for i in range(max(0, pi - r), min(self.ny, pi + r + 1)):
            for j in range(max(0, pj - r), min(self.nx, pj + r + 1)):
                if self.valid[i, j] and self.accum[i, j] > best:
                    best = self.accum[i, j]
        return best * self.cell_area
