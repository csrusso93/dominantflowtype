# -*- coding: utf-8 -*-
"""End-to-end orchestration of the dominantflowtype workflow."""
from __future__ import annotations

import os

import numpy as np

from ._compat import (
    _HAS_QGIS, QgsCoordinateReferenceSystem, QgsCoordinateTransform,
    QgsProject, QgsPointXY, QgsGeometry,
)
from .config import Config
from .io_prompts import (
    Prompter, PromptCancelled, RASTER_FILTER, VECTOR_FILTER, CLOUD_FILTER,
)
from .raster import RasterSampler, build_dod, reproject_raster
from .geometry import (
    _load_vector, _merged_lines, _longest, cast_transects, trimline_hits,
)
from .cross_section import measure_cross_section
from .hydrology import FlowRouter
from .rainfall import rainfall_i30
from .opentopo import opentopo_search, opentopo_download
from . import aoi as _aoi
from . import dem_3dep
from .metrics import compute_qstar, usability
from .outputs import write_workbook, add_qgis_layers


class Pipeline:
    def __init__(self, cfg=None, prompter=None):
        self.cfg = cfg or Config()
        self.p = prompter or Prompter()

    # -- interactive input collection ----------------------------------------
    def collect_inputs(self):
        cfg = self.cfg
        print("=" * 72)
        print(" dominantflowtype — post-event debris-flow vs flood (Q*) workflow")
        print(" Cavagnaro et al. (2024), GRL, doi:10.1029/2024GL109768")
        print("=" * 72)

        # 0. OPTIONAL point-cloud pathway: pre/post clouds -> QGIS M3C2 -> a
        #    post DTM (for Q*) and a change raster (for IDFVA volume). If only a
        #    post cloud is given, no differencing happens and only Q* is produced.
        self._from_clouds = False
        self._change_raster = None
        self._hydrology_dem = None        # full-coverage DEM for A_us + IDFVA
        if self.p.ask_yes_no(
                "\n[0] Use point clouds (LiDAR/SfM) instead of ready DEMs? "
                "Runs QGIS M3C2 and rasterises to a DTM", default=False):
            self._collect_clouds()

        # 1. OPTIONAL DEM download (USGS 3DEP, account-free, or OpenTopography).
        #    Entirely bypassable: if you already have your pre-/post-event DEMs,
        #    answer "no" and go straight to [2]. When used, choose what it is FOR:
        #      * 'pre'     -> use it as the pre-event DEM (obtain pre-event
        #                     topography you don't otherwise have);
        #      * 'routing' -> a full-coverage DEM for A_us flow routing only
        #                     (leaves your pre/post analysis DEMs untouched);
        #      * 'save'    -> just download a copy, don't wire it into anything.
        self._refdem_source = "none"
        self._refdem_purpose = "none"
        if not self._from_clouds and self.p.ask_yes_no(
                "\n[1] (optional) Download a DEM (USGS 3DEP / OpenTopography)? "
                "Skip if you already have pre- and post-event DEMs",
                default=False):
            s = self.p.ask_text("    Source: '3dep' (no key) or 'opentopo'",
                                default="3dep").strip().lower()
            self._refdem_source = "opentopo" if s.startswith(("o", "ot")) else "3dep"
            self._refdem_purpose = self.p.ask_text(
                "    Use it as: 'pre' (pre-event DEM), 'routing' (A_us routing "
                "DEM), or 'save' (just download)", default="routing"
            ).strip().lower()
            if self._refdem_source == "3dep":
                cfg.dem3dep_resolution = int(self.p.ask_float(
                    "    3DEP resolution (m: 1/3/5/10/30)",
                    default=cfg.dem3dep_resolution))
            else:
                cfg.opentopo_api_key = self.p.ask_text(
                    "    OpenTopography API key", allow_blank=True)
                cfg.opentopo_demtype = self.p.ask_text(
                    "    DEM type (USGS1m/USGS10m/COP30/SRTMGL1)",
                    default=cfg.opentopo_demtype)

        # 2 & 3. DEMs (skippable). If the download at [1] is set to be the
        #        pre-event DEM, [2] is filled from that download instead of asked.
        pre_from_dl = self._refdem_purpose == "pre" and self._refdem_source != "none"
        if self._from_clouds:
            # DTM already produced from the post cloud by _collect_clouds()
            pre, post = self.pre_path, self.post_path
            print(f"[2/3] Using DTM rasterised from the post-event cloud: {post}")
        elif pre_from_dl:
            print("[2] Pre-event DEM -> will use the DEM downloaded at [1].")
            pre = None
            post = self.p.ask_path("[3] Post-event DEM (GeoTIFF)", skippable=True,
                                   file_filter=RASTER_FILTER)
        else:
            pre = self.p.ask_path("[2] Pre-event DEM (GeoTIFF)", skippable=True,
                                  file_filter=RASTER_FILTER)
            post = self.p.ask_path("[3] Post-event DEM (GeoTIFF)", skippable=True,
                                   file_filter=RASTER_FILTER)
        if pre is None and post is None and not pre_from_dl and not self._from_clouds:
            raise PromptCancelled("At least one DEM is required.")
        self.pre_path, self.post_path = pre, post
        self._pre_from_dl = pre_from_dl
        self.have_both = pre is not None and post is not None
        if not self.have_both and not pre_from_dl and not self._from_clouds:
            print("    -> Only one DEM provided: volume/DoD products will be "
                  "SKIPPED; Q* computed on the single DEM.")

        # 4-6. vectors
        self.trim_path = self.p.ask_path("[4] Trimlines (.shp/.gpkg)",
                                         file_filter=VECTOR_FILTER)
        self.thal_path = self.p.ask_path("[5] Thalweg (.shp/.gpkg)",
                                         file_filter=VECTOR_FILTER)
        self.ws_path = self.p.ask_path("[6] Watershed outline (.shp/.gpkg)",
                                       file_filter=VECTOR_FILTER)

        # output dir
        base = os.path.dirname(post or pre)
        cfg.output_dir = self.p.ask_dir(
            "Output directory", default=os.path.join(base, "dft_outputs"))
        os.makedirs(cfg.output_dir, exist_ok=True)

        # 7. rainfall (for the 30-min intensity I30)
        print("\n[7] Rainfall — 30-min intensity I30 for Q*")
        src = self.p.ask_text(
            "    Source: 'mrms' (radar, no account), 'synoptic' (gauges), "
            "or 'constant'", default=cfg.rainfall_source).strip().lower()
        cfg.rainfall_source = src if src in ("mrms", "synoptic", "constant") \
            else "synoptic"
        self.rain_station = ""
        self.rain_date = ""
        self.rain_query = "timeseries"
        if cfg.rainfall_source == "mrms":
            # gauge-free NOAA MRMS radar: only an event date is needed; the AOI
            # comes from the watershed. Solves the "no Synoptic access" blocker.
            self.rain_date = self.p.ask_text(
                "    Event/storm date (YYYY-MM-DD, UTC)", allow_blank=True)
            cfg.mrms_reduce = self.p.ask_text(
                "    Reduce field to I30 by (areal_max/areal_mean/point)",
                default=cfg.mrms_reduce)
        elif cfg.rainfall_source == "synoptic":
            # API key only: a request token is auto-minted from it via /v2/auth.
            cfg.synoptic_api_key = self.p.ask_text(
                "    Synoptic API key", default=cfg.synoptic_api_key,
                allow_blank=True)
            self.rain_station = self.p.ask_text(
                "    Station ID (blank = auto-discover gauges in the watershed)",
                allow_blank=True)
            self.rain_date = self.p.ask_text(
                "    Event date (YYYY-MM-DD)", allow_blank=True)
        cfg.rain_event_date = self.rain_date

    # -- point-cloud front-end (Stage 0) -------------------------------------
    def _collect_clouds(self):
        r"""Prompt for pre/post point clouds, run Stage 0 (QGIS M3C2 -> DTM +
        change raster), and wire the results in as if a post DEM was supplied.

        Post-only (no pre cloud) ⇒ no differencing: a DTM is produced for Q\*,
        volume is skipped, and no IDFVA bundle is written.
        """
        from . import pointcloud as _pc
        from . import acquire as _acq
        cfg = self.cfg
        print("\n[PC] Point-cloud pathway — QGIS 4.0 M3C2 (Lague et al., 2013)")
        pre_cloud = self.p.ask_path(
            "[PC-1] Pre-event point cloud (LAS/LAZ/COPC) — needed for VOLUME; "
            "skip for Q*-only", skippable=True, file_filter=CLOUD_FILTER)
        # optional: acquire a pre-event cloud from USGS 3DEP LiDAR if none given
        if pre_cloud is None and self.p.ask_yes_no(
                "[PC-1b] No pre-event cloud — fetch USGS 3DEP LiDAR for an AOI?",
                default=False):
            ws = self.p.ask_path("    Watershed/AOI (.shp/.gpkg) for the 3DEP clip",
                                 file_filter=VECTOR_FILTER)
            try:
                _acq.list_3dep_projects(ws)
            except Exception as ex:
                print(f"    [acquire] project lookup failed ({ex}); enter a URL "
                      f"from the USGS 3DEP LidarExplorer.")
            ept = self.p.ask_text("    EPT/COPC URL (ept.json or .copc.laz)",
                                  allow_blank=True)
            if ept:
                dest = os.path.join(os.path.dirname(ws), "pre_3dep.copc.laz")
                pre_cloud = _acq.acquire_3dep_cloud(ept, ws, dest, cfg)
        post_cloud = self.p.ask_path(
            "[PC-2] Post-event point cloud (LAS/LAZ/COPC)", file_filter=CLOUD_FILTER)
        ptype = self.p.ask_text(
            "[PC-3] Post-event cloud type: 'sfm' or 'lidar'",
            default="sfm").strip().lower()
        post_is_sfm = ptype.startswith("s")

        stage0_dir = os.path.join(os.path.dirname(post_cloud), "stage0_pointcloud")
        os.makedirs(stage0_dir, exist_ok=True)

        # Full-coverage DEM for A_us routing + IDFVA hydrology. The SfM DTM covers
        # only the channel, so a watershed-wide DEM (e.g. USGS 3DEP 10 m) is needed;
        # reproject it to work_epsg so everything shares one CRS.
        routing = self.p.ask_path(
            "[PC-4] Full-coverage watershed DEM for A_us routing + IDFVA hydrology "
            "(e.g. 3DEP 10 m) — recommended", skippable=True,
            file_filter=RASTER_FILTER)
        if routing:
            if cfg.work_epsg:
                routing = reproject_raster(
                    routing, os.path.join(stage0_dir, "routing_dem.tif"),
                    cfg.work_epsg)
            cfg.routing_dem_path = routing
            self._hydrology_dem = routing
        else:
            print("    -> no watershed DEM given: A_us and IDFVA hydrology will use "
                  "the channel-only DTM (catchment will be TRUNCATED).")

        # optional: clip both clouds to an AOI polygon first (fast first pass)
        clip_overlay = self.p.ask_path(
            "[PC-5] Clip clouds to an AOI polygon (watershed) before processing? "
            "Recommended for a fast first run — skip to use the full cloud",
            skippable=True, file_filter=VECTOR_FILTER)

        out = _pc.prepare_from_clouds(pre_cloud, post_cloud, post_is_sfm,
                                      stage0_dir, cfg, clip_overlay=clip_overlay)
        if not out.get("post_dtm"):
            raise PromptCancelled(
                "Stage 0 produced no post DTM — run "
                "dominantflowtype.pointcloud.diagnose() to check algorithm ids, "
                "or set cfg.pointcloud_dry_run=False.")

        self._from_clouds = True
        self.pre_path = None                 # no pre DEM: volume via IDFVA change
        self.post_path = out["post_dtm"]
        self._change_raster = out.get("change_raster")
        self._pre_from_dl = False
        self._refdem_source = self._refdem_purpose = "none"
        self.have_both = False
        if self._change_raster is None:
            print("    -> post-only clouds: differencing/volume SKIPPED; Q* only.")

    # -- main run ------------------------------------------------------------
    def run(self):
        if not _HAS_QGIS:
            raise RuntimeError("This pipeline must run inside QGIS (qgis.core).")
        cfg = self.cfg
        self.collect_inputs()

        # ---- optional DEM download FIRST (may supply the pre-event DEM) -----
        ws_layer = _load_vector(self.ws_path)
        self._download_optional_dem(ws_layer)

        # ---- load reference bed DEM & working CRS --------------------------
        if self.post_path is None and self.pre_path is None:
            raise PromptCancelled("No usable DEM (the [1] download failed and no "
                                  "pre/post DEM was provided).")
        bed_path = self.post_path or self.pre_path      # trimline WSE from post
        bed_unit = cfg.post_vertical_unit if self.post_path else cfg.pre_vertical_unit
        post_sampler = RasterSampler(bed_path, unit=bed_unit)
        pre_sampler = (RasterSampler(self.pre_path, unit=cfg.pre_vertical_unit)
                       if self.pre_path else None)
        work_crs = post_sampler.crs
        work_wkt = post_sampler.wkt

        # ---- DoD -----------------------------------------------------------
        dod_sampler = None
        if self.have_both:
            dod_path = os.path.join(cfg.output_dir, "DoD_post_minus_pre.tif")
            print("\n[proc] building DEM-of-Difference (post - pre) ...")
            dod_sampler = build_dod(self.pre_path, self.post_path, dod_path,
                                    pre_unit=cfg.pre_vertical_unit,
                                    post_unit=cfg.post_vertical_unit)
            # the DoD doubles as the change raster for the IDFVA bridge
            self._change_raster = dod_path

        # ---- load thalweg + trimlines in work CRS --------------------------
        thal_layer = _load_vector(self.thal_path)
        trim_layer = _load_vector(self.trim_path)
        basin = cfg.basin_id or self._read_basin(thal_layer, trim_layer)
        cfg.basin_id = basin

        thal_lines = _merged_lines(thal_layer, work_crs)
        thal_geom = _longest(thal_lines)
        trimlines = _merged_lines(trim_layer, work_crs)

        # ---- transects -----------------------------------------------------
        print("[proc] casting transects every %.1f m ..." % cfg.transect_spacing)
        transects = cast_transects(thal_geom, cfg)
        print(f"       {len(transects)} transects")

        # ---- rainfall I30 --------------------------------------------------
        # One dispatcher, source chosen by cfg.rainfall_source: 'mrms' (gauge-free
        # radar, AOI = watershed), 'synoptic' (a named gauge or auto-discovered
        # gauges in the watershed), or 'constant'. AOI + a representative basin
        # point (watershed centre, WGS84) drive the AOI-based paths.
        print(f"[proc] retrieving rainfall intensity "
              f"(I30, source={cfg.rainfall_source}) ...")
        try:
            ws_center = _aoi.bounds_center(_aoi.load_aoi(ws_layer))
        except Exception:                              # noqa: BLE001
            ws_center = None
        i30, rain_meta = rainfall_i30(
            cfg, station_id=self.rain_station or None,
            date_str=self.rain_date or None, aoi=ws_layer, point=ws_center)
        print(f"       I30 = {i30:.3g} mm/hr")

        # ---- flow accumulation (A_us) --------------------------------------
        # Prefer an explicit routing DEM (e.g. the just-downloaded full-coverage
        # 3DEP/OpenTopography DEM) over the analysis DEM, which may be clipped to
        # the survey footprint and undercount A_us near the watershed edge.
        routing_dem = getattr(cfg, "routing_dem_path", "") or \
            self.pre_path or self.post_path
        if getattr(cfg, "routing_dem_path", ""):
            print(f"[proc] computing upstream drainage area (D8 routing on "
                  f"{os.path.basename(routing_dem)}) ...")
        else:
            print("[proc] computing upstream drainage area (D8 flow routing) ...")
        router = FlowRouter(routing_dem, ws_layer, cfg, work_wkt)

        # ---- per-transect measurement + metrics ----------------------------
        records = []
        wgs = QgsCoordinateReferenceSystem("EPSG:4326")
        to_wgs = QgsCoordinateTransform(work_crs, wgs, QgsProject.instance())

        # precompute trimline intersections once per transect (cached)
        hits_list = [trimline_hits(tr, trimlines, cfg, post_sampler)
                     for tr in transects]

        prev_elev = None
        prev_s = None
        for ti, tr in enumerate(transects):
            wse, continuity, uL_map, uR_map = self._resolve_wse(
                tr, ti, hits_list, transects, cfg)
            if wse is None or not np.isfinite(wse):
                continue

            # post-event (or single-DEM) cross-section
            xs_post = measure_cross_section(tr, wse, uL_map, uR_map,
                                            post_sampler, cfg)
            if xs_post is None:
                continue
            A_us = router.upstream_area(tr["x"], tr["y"])

            # bed change for usability + volume
            dz_bed = float("nan")
            if dod_sampler is not None:
                dz_bed = dod_sampler.sample(tr["x"], tr["y"])

            m_post = compute_qstar(
                xs_post["A_xs"], xs_post["wetted_perimeter"],
                xs_post["hydraulic_radius"], xs_post["hydraulic_depth"],
                xs_post["max_depth"], A_us, i30, cfg)

            # slope from thalweg bed elevation gradient
            bed_elev = post_sampler.sample(tr["x"], tr["y"])
            slope = float("nan")
            if prev_elev is not None and np.isfinite(bed_elev) and \
                    np.isfinite(prev_elev) and (tr["s"] - prev_s) > 0:
                slope = abs(prev_elev - bed_elev) / (tr["s"] - prev_s)
            prev_elev, prev_s = bed_elev, tr["s"]

            pt = to_wgs.transform(QgsPointXY(tr["x"], tr["y"]))
            usable_post, reason_post, _ = usability(
                xs_post["max_depth"], dz_bed, continuity, cfg)

            rec = dict(
                transect_id=tr["idx"], s=tr["s"], lat=pt.y(), lon=pt.x(),
                slope=slope, continuity=continuity,
                # transect line trimmed to the inundation boundary (mapped
                # trimline crossings, or equal-elevation crossings in
                # discontinuous sections) rather than the full casting width
                transect_geom=QgsGeometry.fromPolylineXY(
                    [QgsPointXY(*xs_post["left_xy"]),
                     QgsPointXY(*xs_post["right_xy"])]),
                left_xy=xs_post["left_xy"], right_xy=xs_post["right_xy"],
                A_us=A_us, dz_bed=dz_bed, wse=wse,
                # post/single metrics
                depth_m=xs_post["max_depth"], width_m=xs_post["top_width"],
                area_xs_m2=xs_post["A_xs"],
                wetted_perimeter_m=xs_post["wetted_perimeter"],
                hydraulic_radius_m=xs_post["hydraulic_radius"],
                velocity=m_post["velocity"], Q_peak=m_post["Q_peak"],
                Q_fluv=m_post["Q_fluv"], Qstar=m_post["Qstar"],
                flow_type=m_post["flow_type"],
                usable=usable_post, flag_reason=reason_post,
            )

            # pre-event cross-section (same WSE + inundation footprint, pre bed)
            if self.have_both and pre_sampler is not None:
                xs_pre = measure_cross_section(
                    tr, wse, xs_post["uL"], xs_post["uR"], pre_sampler, cfg)
                if xs_pre is not None:
                    m_pre = compute_qstar(
                        xs_pre["A_xs"], xs_pre["wetted_perimeter"],
                        xs_pre["hydraulic_radius"], xs_pre["hydraulic_depth"],
                        xs_pre["max_depth"], A_us, i30, cfg)
                    up_pre, rs_pre, _ = usability(
                        xs_pre["max_depth"], dz_bed, continuity, cfg)
                    rec.update(
                        pre_depth_m=xs_pre["max_depth"],
                        pre_width_m=xs_pre["top_width"],
                        pre_area_xs_m2=xs_pre["A_xs"],
                        pre_wetted_perimeter_m=xs_pre["wetted_perimeter"],
                        pre_hydraulic_radius_m=xs_pre["hydraulic_radius"],
                        pre_velocity=m_pre["velocity"], pre_Q_peak=m_pre["Q_peak"],
                        pre_Q_fluv=m_pre["Q_fluv"], pre_Qstar=m_pre["Qstar"],
                        pre_flow_type=m_pre["flow_type"],
                        pre_usable=up_pre, pre_flag_reason=rs_pre)
            records.append(rec)

        print(f"[proc] measured {len(records)} cross-sections")

        # ---- inundation & reach polygons (volume) --------------------------
        reaches = self._build_reaches(records, dod_sampler, work_crs)
        self._attach_polygons(records, work_crs)
        self._attach_reach_volume(records, reaches)

        # ---- write workbook ------------------------------------------------
        rows_pre, rows_post = self._rows(records)
        out_xlsx = os.path.join(cfg.output_dir,
                                f"{basin or 'basin'}_dominantflowtype.xlsx")
        write_workbook(rows_pre, rows_post, cfg, out_xlsx, rain_meta)

        # ---- QGIS layers ---------------------------------------------------
        if cfg.add_layers_to_qgis:
            add_qgis_layers(self._render_records(records), reaches, work_crs,
                            cfg, cfg.output_dir, self.have_both)

        # ---- IDFVA bridge (volume delegated to Guido's IDFVA, separate venv) --
        if cfg.write_idfva_bridge and self._change_raster:
            from .bridge import export_bridge
            # IDFVA hydrology needs a full-watershed DEM; use the routing DEM when
            # available (the SfM DTM covers only the channel), else the post DTM.
            hydrology_dem = self._hydrology_dem or self.post_path
            export_bridge(
                os.path.join(cfg.output_dir, "bridge_bundle"),
                thal_layer, hydrology_dem, self._change_raster, work_crs,
                meta={"basin_id": basin, "from_point_clouds": self._from_clouds,
                      "hydrology_dem": os.path.basename(hydrology_dem),
                      "change_from": "SfM/LiDAR M3C2" if self._from_clouds else "DoD"},
                epsg_hint=(cfg.work_epsg or None))

        print("\nDONE. Outputs in:", cfg.output_dir)
        return records, reaches, rain_meta

    # -- helpers -------------------------------------------------------------
    def _download_optional_dem(self, ws_layer):
        """Handle the optional [1] DEM download and wire it in per its purpose.

        Purpose ('pre' | 'routing' | 'save') was chosen at [1]. 'pre' sets it as
        the pre-event DEM (so users lacking pre-event topography can obtain it);
        'routing' sets ``cfg.routing_dem_path`` (full-coverage A_us only, pre/post
        untouched); 'save' just leaves the file on disk. A failed download degrades
        gracefully (single-DEM mode if it was to be the pre-event DEM).
        """
        cfg = self.cfg
        src = getattr(self, "_refdem_source", "none")
        if src == "none":
            return
        purpose = getattr(self, "_refdem_purpose", "routing")
        path = None
        if src == "3dep":
            bbox = _aoi.load_aoi(ws_layer, pad_deg=0.01)        # (W,S,E,N)
            dl = os.path.join(cfg.output_dir,
                              f"3dep_{cfg.dem3dep_resolution}m.tif")
            print(f"\n[proc] downloading USGS 3DEP {cfg.dem3dep_resolution} m DEM ...")
            path = dem_3dep.get_dem(bbox, resolution=cfg.dem3dep_resolution,
                                    out_path=dl, dst_epsg=cfg.dem3dep_epsg)
        elif src == "opentopo":
            bbox = self._aoi_wgs84(ws_layer)                    # (S,W,N,E)
            opentopo_search(bbox, cfg.opentopo_api_key)
            dl = os.path.join(cfg.output_dir,
                              f"opentopo_{cfg.opentopo_demtype}.tif")
            path = opentopo_download(bbox, cfg.opentopo_api_key,
                                     cfg.opentopo_demtype, dl)
        if not path:
            print("       [1] DEM download failed; continuing without it.")
            if purpose == "pre":
                self._pre_from_dl = False       # fall back to single-DEM on post
            return
        if purpose == "pre":
            self.pre_path = path
            self.have_both = self.post_path is not None
            print(f"       -> using downloaded DEM as the PRE-event DEM.")
        elif purpose == "routing":
            cfg.routing_dem_path = path
            print(f"       -> using downloaded DEM for A_us ROUTING only.")
        else:                                    # save
            print(f"       -> saved DEM (not wired into the analysis): {path}")

    def _aoi_wgs84(self, layer):
        ext = layer.extent()
        wgs = QgsCoordinateReferenceSystem("EPSG:4326")
        xform = QgsCoordinateTransform(layer.crs(), wgs, QgsProject.instance())
        ll = xform.transform(QgsPointXY(ext.xMinimum(), ext.yMinimum()))
        ur = xform.transform(QgsPointXY(ext.xMaximum(), ext.yMaximum()))
        return (ll.y(), ll.x(), ur.y(), ur.x())     # s, w, n, e

    def _read_basin(self, *layers):
        for lyr in layers:
            names = [f.name() for f in lyr.fields()]
            if "Basin" in names:
                for f in lyr.getFeatures():
                    if f["Basin"]:
                        return str(f["Basin"])
        return ""

    def _resolve_wse(self, tr, ti, hits_list, transects, cfg):
        """Determine water-surface elevation + inundation edges for a transect.

        - both banks mapped  -> WSE = lower of the two (paper);   'certain'
        - one bank mapped     -> WSE = that bank; other edge inferred at same
                                 elevation (constant depth);       'inferred_1'
        - neither mapped      -> WSE interpolated from nearest mapped transects;
                                 both edges inferred;              'inferred_2'
        Returns (wse, continuity, uL_map, uR_map) where u*_map is the mapped
        along-transect station on that side or None (=> inferred edge).
        """
        hits = hits_list[ti]
        L = hits["L"]
        R = hits["R"]
        if L is not None and R is not None and np.isfinite(L[3]) and np.isfinite(R[3]):
            wse = min(L[3], R[3])
            return wse, "certain", L[0], R[0]
        if L is not None and np.isfinite(L[3]):
            return L[3], "inferred_1side", L[0], None
        if R is not None and np.isfinite(R[3]):
            return R[3], "inferred_1side", None, R[0]
        # neither: interpolate WSE from nearest mapped neighbours
        wse = self._interp_wse(tr, hits_list, transects)
        return wse, "inferred_2side", None, None

    def _interp_wse(self, tr, hits_list, transects):
        """Interpolate WSE at a gap transect from nearest up/down mapped ones."""
        s0 = tr["s"]
        up = down = None
        for other, hh in zip(transects, hits_list):
            elevs = [h[3] for h in (hh["L"], hh["R"])
                     if h is not None and np.isfinite(h[3])]
            if not elevs:
                continue
            e = min(elevs)
            if other["s"] <= s0 and (up is None or other["s"] > up[0]):
                up = (other["s"], e)
            if other["s"] >= s0 and (down is None or other["s"] < down[0]):
                down = (other["s"], e)
        if up and down and down[0] != up[0]:
            t = (s0 - up[0]) / (down[0] - up[0])
            return up[1] + t * (down[1] - up[1])
        if up:
            return up[1]
        if down:
            return down[1]
        return float("nan")

    def _build_reaches(self, records, dod_sampler, crs):
        """Trapezoidal inundation polygons between consecutive transects; DoD
        zonal mean * area = volume (post-pre; +deposition, -erosion)."""
        reaches = []
        for i in range(len(records) - 1):
            a, b = records[i], records[i + 1]
            try:
                poly = QgsGeometry.fromPolygonXY([[
                    QgsPointXY(*a["left_xy"]), QgsPointXY(*b["left_xy"]),
                    QgsPointXY(*b["right_xy"]), QgsPointXY(*a["right_xy"]),
                    QgsPointXY(*a["left_xy"]),
                ]])
            except Exception:
                poly = None
            area = poly.area() if poly else float("nan")
            mean_dz = float("nan")
            volume = float("nan")
            if dod_sampler is not None and poly is not None:
                mean_dz = self._zonal_mean(poly, dod_sampler)
                if np.isfinite(mean_dz):
                    volume = mean_dz * area
            reaches.append(dict(reach_id=i, polygon=poly, area=area,
                                mean_dz=mean_dz, volume=volume,
                                t_from=a["transect_id"], t_to=b["transect_id"]))
        return reaches

    def _zonal_mean(self, poly, sampler):
        """Mean raster value inside a polygon (grid sampling)."""
        bb = poly.boundingBox()
        step = max(sampler.xres, 0.5)
        xs = np.arange(bb.xMinimum(), bb.xMaximum() + step, step)
        ys = np.arange(bb.yMinimum(), bb.yMaximum() + step, step)
        eng = QgsGeometry(poly)
        vals = []
        for y in ys:
            for x in xs:
                if eng.contains(QgsGeometry.fromPointXY(QgsPointXY(x, y))):
                    v = sampler.sample(x, y)
                    if np.isfinite(v):
                        vals.append(v)
        return float(np.mean(vals)) if vals else float("nan")

    def _attach_polygons(self, records, crs):
        """Give each transect an inundation polygon (half-reach up + half down)."""
        for i, r in enumerate(records):
            # simple ribbon: connect to the next transect's waterline points
            if i + 1 < len(records):
                nxt = records[i + 1]
                poly = QgsGeometry.fromPolygonXY([[
                    QgsPointXY(*r["left_xy"]), QgsPointXY(*nxt["left_xy"]),
                    QgsPointXY(*nxt["right_xy"]), QgsPointXY(*r["right_xy"]),
                    QgsPointXY(*r["left_xy"])]])
            else:
                poly = None
            r["polygon"] = poly

    def _attach_reach_volume(self, records, reaches):
        vol_by_from = {rc["t_from"]: rc["volume"] for rc in reaches}
        for r in records:
            r["volume_m3"] = vol_by_from.get(r["transect_id"], float("nan"))

    def _rows(self, records):
        rows_pre, rows_post = [], []
        b = self.cfg.basin_id
        for r in records:
            if self.have_both and "pre_area_xs_m2" in r:
                rows_pre.append([
                    b, r["transect_id"], r["lat"], r["lon"], r["s"], r["slope"],
                    r["pre_depth_m"], r["pre_width_m"], r["pre_area_xs_m2"],
                    r["pre_wetted_perimeter_m"], r["pre_hydraulic_radius_m"],
                    r["A_us"], r["pre_velocity"], r["pre_Q_peak"],
                    r["pre_Q_fluv"], r["pre_Qstar"], r["pre_flow_type"],
                    r["continuity"], r["pre_usable"], r["pre_flag_reason"]])
            post_row = [
                b, r["transect_id"], r["lat"], r["lon"], r["s"], r["slope"],
                r["depth_m"], r["width_m"], r["area_xs_m2"],
                r["wetted_perimeter_m"], r["hydraulic_radius_m"], r["A_us"],
                r["velocity"], r["Q_peak"], r["Q_fluv"], r["Qstar"],
                r["flow_type"], r["continuity"]]
            # POST_EXTRA inserted before usable
            if self.have_both:
                post_row += [r["dz_bed"], r.get("volume_m3", float("nan"))]
                rows_post.append(post_row + [r["usable"], r["flag_reason"]])
            else:
                # single DEM: no volume; treat this as the pre list too
                rows_pre.append([
                    b, r["transect_id"], r["lat"], r["lon"], r["s"], r["slope"],
                    r["depth_m"], r["width_m"], r["area_xs_m2"],
                    r["wetted_perimeter_m"], r["hydraulic_radius_m"], r["A_us"],
                    r["velocity"], r["Q_peak"], r["Q_fluv"], r["Qstar"],
                    r["flow_type"], r["continuity"], r["usable"],
                    r["flag_reason"]])
        return rows_pre, rows_post

    def _render_records(self, records):
        out = []
        for r in records:
            out.append(dict(
                transect_id=r["transect_id"], polygon=r.get("polygon"),
                transect_geom=r.get("transect_geom"),
                Qstar_post=r.get("Qstar"), flow_type_post=r.get("flow_type"),
                area_xs_m2_post=r.get("area_xs_m2"),
                depth_m_post=r.get("depth_m"), width_m_post=r.get("width_m"),
                usable_post=r.get("usable"), Qstar=r.get("Qstar"),
                flow_type=r.get("flow_type"), area_xs_m2=r.get("area_xs_m2"),
                depth_m=r.get("depth_m"), width_m=r.get("width_m"),
                usable=r.get("usable")))
        return out


def run(cfg=None, answers=None):
    """Interactive entry point. In the QGIS console simply call: run().

    Cancelling any prompt raises ``PromptCancelled``, which is caught here and
    turned into a clean ``None`` return — it must never escape as ``SystemExit``,
    which would terminate the QGIS application.
    """
    prompter = Prompter(answers=answers) if answers else Prompter()
    try:
        return Pipeline(cfg=cfg, prompter=prompter).run()
    except PromptCancelled as ex:
        print(f"\n[dominantflowtype] cancelled — {ex}. No changes made.")
        return None
