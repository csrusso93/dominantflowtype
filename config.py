# -*- coding: utf-8 -*-
"""Configuration for the dominantflowtype workflow."""
from __future__ import annotations

from dataclasses import dataclass, asdict, field


@dataclass
class Config:
    """All tunable parameters. Defaults follow Cavagnaro et al. (2024)."""

    # --- physical constants / method assumptions -----------------------------
    g: float = 9.81                       # gravitational acceleration [m/s^2]
    froude: float = 1.0                   # Froude-critical assumption (Fr=1)
    velocity_scale: str = "hydraulic_radius"   # 'hydraulic_radius' (paper 'h'),
                                               # 'hydraulic_depth' (A/T) or 'max_depth'
    qstar_threshold: float = 1.0          # Q* > threshold => debris flow

    # --- geometry ------------------------------------------------------------
    transect_spacing: float = 5.0         # cast transects every N metres
    transect_halfwidth: float = 50.0      # half-length of each transect [m]
    profile_step: float = 0.5             # DEM sampling step along transect [m]
    tangent_delta: float = 2.5            # +/- distance for tangent estimate [m]

    # --- trimline handling ---------------------------------------------------
    # Name of an OPTIONAL attribute flagging trimline reliability. Recognised
    # values (case-insensitive) that mark a MAPPED/observed trimline:
    trimline_flag_field: str = "type"     # optional; auto-ignored if absent
    trimline_certain_values = ("certain", "mapped", "observed", "1", "true", "yes")
    trimline_snap_tol: float = 2.0        # trimline<->transect snap tolerance [m]

    # --- usability (Cavagnaro Fig 1d/1e: incision/deposition) ----------------
    # A cross-section is UNUSABLE when |bed change| is large relative to depth.
    # Discontinuity of the trimline does NOT make a section unusable.
    max_dz_to_depth_ratio: float = 0.5    # |dz|/depth above this => unusable
    min_depth: float = 0.05               # ignore sections shallower than this [m]

    # --- flow accumulation (A_us) --------------------------------------------
    routing_dem_path: str = ""            # optional full-coverage DEM for A_us
                                          # (set automatically if 3DEP/OpenTopo
                                          # is downloaded); falls back to pre/post
    flow_accum_res: float = 10.0          # resample DEM to this cell size [m]
    flow_accum_buffer: float = 50.0       # buffer watershed before clipping [m]
    snap_radius: float = 15.0             # snap station to max-accum cell within [m]
    use_grass_if_available: bool = False  # prefer bundled numpy router (portable)

    # --- rainfall ------------------------------------------------------------
    # I30 source, in order of preference. 'mrms' = gauge-free NOAA MRMS radar
    # (no account; see dominantflowtype.mrms); 'synoptic' = ground gauges (needs a
    # Synoptic account with data access); 'constant' = the Cavagnaro fallback.
    rainfall_source: str = "synoptic"     # 'synoptic' | 'mrms' | 'constant'
    rain_event_date: str = ""             # 'YYYY-MM-DD' storm day (MRMS/Synoptic)
    rain_window_duration_min: int = 30    # I30 window [min] (15/30/60 for MRMS)
    default_i30_mm_hr: float = 9.7        # Cavagnaro optimised constant fallback
    # Synoptic (ground gauges)
    synoptic_token: str = ""              # public token (per-request)
    synoptic_api_key: str = ""            # master api key (optional)
    # MRMS (radar): how to reduce the AOI peak-intensity field to one I30 number
    mrms_reduce: str = "areal_max"        # 'areal_max' | 'areal_mean' | 'point'

    # --- opentopography / 3DEP ----------------------------------------------
    opentopo_api_key: str = ""
    opentopo_demtype: str = "USGS10m"     # e.g. USGS1m / USGS10m / COP30 / SRTMGL1
    # USGS 3DEP (The National Map) reference/routing DEM (see dem_3dep.py)
    dem3dep_resolution: int = 10          # metres (1/3/5/10/30/60)
    dem3dep_epsg: int = 5070              # output CRS (5070 = CONUS Albers, m)

    # --- DEM vertical units --------------------------------------------------
    # Vertical unit of each input DEM: 'm' (metres), 'ft' (US/international
    # survey feet), or 'auto'. Feet are converted to metres before differencing
    # so bed change is always reported in metres. 'auto' inspects the raster CRS
    # (vertical axis unit, else the horizontal linear unit as a heuristic) and
    # converts when it looks like feet; set 'm'/'ft' explicitly to override.
    pre_vertical_unit: str = "auto"       # 'm' | 'ft' | 'auto'
    post_vertical_unit: str = "auto"      # 'm' | 'ft' | 'auto'

    # --- Stage 0: point clouds & M3C2 (see pointcloud.py) --------------------
    # QGIS 4.x PDAL algorithm ids (verified on QGIS 4.2). Override to pin a
    # different build; run dominantflowtype.pointcloud.diagnose() to list ids.
    dtm_resolution: float = 1.0           # cloud -> raster cell size [m]
    # Common-CRS reproject before M3C2 — pdal:reproject (0 = leave clouds as-is)
    work_epsg: int = 0                    # e.g. 6339; reprojects pre+post to EPSG:<n>
    reproject_alg_id: str = "pdal:reproject"
    reproject_params: dict = field(default_factory=dict)
    # SfM noise pre-pass — pdal:filternoisestatistical (statistical outlier removal)
    denoise_sfm: bool = True              # run on an SfM post cloud before ground
    denoise_alg_id: str = "pdal:filternoisestatistical"
    denoise_mean_k: int = 8               # -> MEAN_K (neighbours)
    denoise_multiplier: float = 3.0       # -> MULTIPLIER (std-dev)
    denoise_params: dict = field(default_factory=dict)
    # M3C2 — pdal:compare (Lague et al., 2013)
    m3c2_alg_id: str = "pdal:compare"
    m3c2_normal_scale: float = 0.0        # -> NORMAL_RADIUS [m] (0 = alg default)
    m3c2_cyl_radius: float = 0.0          # -> CYLINDER_RADIUS [m] (0 = alg default)
    m3c2_subsample: float = 0.0           # -> SUBSAMPLING_CELL_SIZE [m] (0 = full)
    m3c2_registration_error: float = 0.0  # -> REGISTRATION_ERROR (0 = alg default)
    m3c2_params: dict = field(default_factory=dict)   # extra/override params
    # Ground classification — pdal:classifyground (morphological; QGIS has no CSF)
    ground_alg_id: str = "pdal:classifyground"
    classify_sfm_ground: bool = True      # classify an SfM post cloud before DTM
    ground_params: dict = field(default_factory=dict)  # CELL_SIZE/SLOPE/WINDOW/...
    dtm_ground_only: bool = True          # rasterize only ground-class points
    ground_filter_expression: str = "Classification == 2"
    # Rasterise — pdal:exportraster (attribute rasteriser, used for m3c2_distance)
    export_raster_alg_id: str = "pdal:exportraster"
    export_raster_params: dict = field(default_factory=dict)
    # DTM rasterisation method: 'tin' = pdal:exportrastertin (Delaunay, gap-free,
    # like IDFVA's interpolation) or 'binning' = pdal:exportraster ATTRIBUTE=Z.
    dtm_method: str = "tin"               # 'tin' | 'binning'
    exportrastertin_alg_id: str = "pdal:exportrastertin"
    dtm_max_edge_length: float = 0.0      # TIN max triangle edge [m] (0 = default)
    # 3DEP LiDAR acquisition — clip an EPT/COPC to the watershed (pdal:clip)
    clip_alg_id: str = "pdal:clip"
    acquire_params: dict = field(default_factory=dict)
    pointcloud_dry_run: bool = False      # print processing.run calls, don't execute

    # --- volume / IDFVA bridge ----------------------------------------------
    # Q* is always computed. Volume is delegated to Guido's IDFVA (separate venv):
    # when True and a change surface exists (cloud M3C2 or DEM DoD), export a
    # bridge bundle (flow path + DTM + change raster) for IDFVA to consume.
    # (The legacy in-package DoD reach-volume still runs for the two-DEM path;
    # gating it behind a flag is Phase 2.)
    write_idfva_bridge: bool = True       # export flow-path + DTM + change raster

    # --- i/o -----------------------------------------------------------------
    output_dir: str = ""                  # set at runtime; defaults next to inputs
    basin_id: str = ""                    # auto-read from data if blank
    add_layers_to_qgis: bool = True

    def as_dict(self):
        d = asdict(self)
        d["velocity_scale"] = self.velocity_scale
        return d
