# -*- coding: utf-8 -*-
r"""Stage 0 — point-cloud preparation & differencing, inside QGIS.

Everything here drives QGIS Processing algorithms (``processing.run``) so no
external tools are needed:

* **Clip** — ``pdal:clip`` trims both clouds to a polygon AOI (watershed) first,
  to shrink huge surveys for a fast pass.
* **Reproject** — ``pdal:reproject`` brings pre + post into a common CRS
  (``Config.work_epsg``) so ``pdal:compare`` can difference them.
* **Denoise** — ``pdal:filternoisestatistical`` removes SfM outlier/ghost points.
* **Ground classification** — ``pdal:classifyground`` (morphological/PMF-style;
  QGIS 4 ships this rather than CSF) marks bare-earth points as class 2 for SfM.
* **M3C2 change** — ``pdal:compare`` = QGIS 4.0 native "Compare Point Clouds"
  (Lague et al., 2013): reference (pre) INPUT vs INPUT_COMPARE (post) -> point
  cloud with an ``m3c2_distance`` scalar (signed, along the reference normal).
* **Export to raster** — **DTM** via ``pdal:exportrastertin`` (TIN, gap-free) or
  ``pdal:exportraster`` (binning); ``m3c2_distance`` -> **change raster** always
  via ``pdal:exportraster`` (attribute rasteriser). DTM filtered to Class==2.

Because Processing algorithm ids and parameter names differ across QGIS point
releases, this module **discovers** the algorithms from the processing registry
at runtime and prints exactly what it will run. Anything can be pinned explicitly
through ``Config`` (``m3c2_alg_id`` / ``m3c2_params`` etc.) to override discovery.

Run ``dominantflowtype.pointcloud.diagnose()`` in the QGIS Python Console to list
the point-cloud algorithms and their parameter names on YOUR build — paste the
output back and the exact ids/params can then be pinned in Config.

References
----------
Lague, Brodu & Leroux (2013), ISPRS J. Photogramm. Remote Sens. 82, 10-26,
  doi:10.1016/j.isprsjprs.2013.04.009 (M3C2).
Zhang et al. (2016), Remote Sensing 8(6), 501, doi:10.3390/rs8060501 (CSF).
Lutra Consulting (2025), "M3C2 point-cloud comparison in QGIS 4.0".
"""
from __future__ import annotations

import os


# --- processing framework access --------------------------------------------
def _processing():
    """Return the QGIS ``processing`` module, or raise a clear error."""
    try:
        import processing  # noqa: F401  (registered in the QGIS console)
        return processing
    except Exception:
        try:
            from qgis import processing  # noqa: F401
            return processing
        except Exception as ex:
            raise RuntimeError(
                "QGIS 'processing' module is unavailable — run inside the QGIS "
                f"Python Console (Stage 0 needs QGIS Processing). [{ex}]")


def _registry():
    from qgis.core import QgsApplication
    return QgsApplication.processingRegistry()


# --- algorithm & parameter discovery ----------------------------------------
def list_pointcloud_algorithms():
    """Return [(id, displayName, group)] for likely point-cloud algorithms."""
    out = []
    for alg in _registry().algorithms():
        hay = f"{alg.id()} {alg.displayName()} {alg.group()}".lower()
        if "pointcloud" in hay or alg.id().startswith("pdal:") or "point cloud" in hay:
            out.append((alg.id(), alg.displayName(), alg.group()))
    return sorted(out)


def describe_algorithm(alg_id):
    """Return (inputs, outputs) parameter-name lists for an algorithm id."""
    alg = _registry().algorithmById(alg_id)
    if alg is None:
        raise ValueError(f"no algorithm with id {alg_id!r}")
    inputs = [(p.name(), p.__class__.__name__, p.description())
              for p in alg.parameterDefinitions()]
    outputs = [(o.name(), o.__class__.__name__) for o in alg.outputDefinitions()]
    return inputs, outputs


def diagnose():
    """Print point-cloud algorithms and their parameters (run this in QGIS).

    Paste the output back so the exact ids/params can be pinned in Config.
    """
    algs = list_pointcloud_algorithms()
    if not algs:
        print("[pointcloud] No point-cloud algorithms found. Is this QGIS 4.x "
              "with the PDAL provider enabled?")
        return
    print(f"[pointcloud] {len(algs)} point-cloud algorithm(s):")
    for aid, name, group in algs:
        print(f"\n  • {aid}   ({group} → {name})")
        try:
            inputs, outputs = describe_algorithm(aid)
            for pn, ptype, desc in inputs:
                print(f"        in  {pn:24s} {ptype:28s} {desc}")
            for on, otype in outputs:
                print(f"        out {on:24s} {otype}")
        except Exception as ex:
            print(f"        (could not introspect: {ex})")


def _find_alg(token_groups, override_id=""):
    """Resolve an algorithm id. ``token_groups`` is a list of token lists; an
    algorithm matches if, for ANY group, every token is a substring of its
    id+displayName+group. Explicit ``override_id`` wins."""
    if override_id:
        if _registry().algorithmById(override_id) is None:
            raise ValueError(f"pinned algorithm id {override_id!r} not found")
        return override_id
    cands = []
    for aid, name, group in list_pointcloud_algorithms():
        hay = f"{aid} {name} {group}".lower()
        for tokens in token_groups:
            if all(t in hay for t in tokens):
                cands.append(aid)
                break
    if not cands:
        raise LookupError(
            "could not find a matching point-cloud algorithm for "
            f"{token_groups}. Run dominantflowtype.pointcloud.diagnose() and pin "
            "the id in Config.")
    if len(cands) > 1:
        print(f"  [pointcloud] multiple matches for {token_groups}: {cands}; "
              f"using {cands[0]}. Pin one in Config to be explicit.")
    return cands[0]


def _param_names(alg_id):
    return [p.name() for p in _registry().algorithmById(alg_id).parameterDefinitions()]


def _pick(names, *token_alts, default=None):
    """First parameter name matching any token-alternative (each may be several
    tokens that must all appear). Case-insensitive substring match."""
    low = {n: n.lower() for n in names}
    for alt in token_alts:
        toks = [alt] if isinstance(alt, str) else list(alt)
        for n in names:
            if all(t.lower() in low[n] for t in toks):
                return n
    return default


def _run(processing, alg_id, params, dry_run, label):
    # QGIS 4.2 PDAL cloud-output algorithms (clip/reproject/classifyground/
    # filternoise*/...) expose a VPC_OUTPUT_FORMAT enum whose default is invalid in
    # some builds -> "Incorrect parameter value". Inject a valid index (0) when the
    # algorithm has that parameter and the caller didn't set it.
    try:
        if ("VPC_OUTPUT_FORMAT" in _param_names(alg_id)
                and "VPC_OUTPUT_FORMAT" not in params):
            params = dict(params, VPC_OUTPUT_FORMAT=0)
    except Exception:
        pass
    print(f"  [pointcloud] {label}: processing.run({alg_id!r}, {params})")
    if dry_run:
        print("     (dry run — not executed)")
        return None
    return processing.run(alg_id, params)


# --- algorithm ids verified on QGIS 4.2 (override via Config) ----------------
def _resolve(cfg_id, fallback_tokens):
    """Return an algorithm id: the pinned ``cfg_id`` if present in the registry,
    else discover it from ``fallback_tokens``."""
    if cfg_id and _registry().algorithmById(cfg_id) is not None:
        return cfg_id
    if cfg_id:
        print(f"  [pointcloud] pinned id {cfg_id!r} not found on this build; "
              f"falling back to discovery.")
    return _find_alg(fallback_tokens)


# --- Stage-0 operations ------------------------------------------------------
def classify_ground(cloud_path, out_path, cfg):
    """Classify bare-earth (ground) points via ``pdal:classifyground``.

    A morphological ground filter (PMF-style); QGIS 4 ships this rather than CSF
    (Zhang et al., 2016). Output keeps ALL points with ground set to class 2 —
    rasterise with ``FILTER_EXPRESSION="Classification == 2"`` for bare earth.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "ground_alg_id", ""),
                      [["classify", "ground"], ["ground"]])
    params = {"INPUT": cloud_path, "OUTPUT": out_path}
    params.update(getattr(cfg, "ground_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               "classify ground (pdal:classifyground)")
    return out_path if res is not None else cloud_path


def clip_cloud(cloud_path, overlay_path, out_path, cfg):
    """Clip a point cloud to a polygon AOI (watershed) via ``pdal:clip``.

    Run this first to shrink a huge cloud to the reach of interest — essential for
    a fast first pass (e.g. 324 M-point SfM surveys). Returns the clipped path.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "clip_alg_id", ""), [["clip"], ["pdal", "clip"]])
    params = {"INPUT": cloud_path, "OVERLAY": overlay_path, "OUTPUT": out_path}
    params.update(getattr(cfg, "acquire_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               "clip cloud to AOI")
    return out_path if res is not None else cloud_path


def reproject_cloud(cloud_path, out_path, epsg, cfg):
    """Reproject a point cloud to ``EPSG:<epsg>`` via ``pdal:reproject``.

    Ensures pre and post clouds share one CRS before M3C2 (``pdal:compare``
    requires it). A no-op-ish pass when the cloud is already in the target CRS.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "reproject_alg_id", ""), [["reproject"]])
    params = {"INPUT": cloud_path, "CRS": f"EPSG:{epsg}", "OUTPUT": out_path}
    params.update(getattr(cfg, "reproject_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               f"reproject to EPSG:{epsg}")
    return out_path if res is not None else cloud_path


def denoise_cloud(cloud_path, out_path, cfg):
    """Statistical outlier removal via ``pdal:filternoisestatistical`` (SfM pre-pass).

    Removes sparse noise points (SfM stray/ghost points) before ground
    classification. ``MEAN_K`` neighbours, ``MULTIPLIER`` std-dev threshold.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "denoise_alg_id", ""),
                      [["noise", "statistical"], ["filter", "noise"], ["noise"]])
    params = {"INPUT": cloud_path, "OUTPUT": out_path, "REMOVE_NOISE_POINTS": True}
    if getattr(cfg, "denoise_mean_k", 0):
        params["MEAN_K"] = cfg.denoise_mean_k
    if getattr(cfg, "denoise_multiplier", 0):
        params["MULTIPLIER"] = cfg.denoise_multiplier
    params.update(getattr(cfg, "denoise_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               "denoise (statistical outlier removal)")
    return out_path if res is not None else cloud_path


def compare_m3c2(reference_cloud, compare_cloud, out_cloud, cfg):
    """M3C2 change via ``pdal:compare`` (reference=pre, INPUT_COMPARE=post).

    Output cloud carries an ``m3c2_distance`` scalar (signed, along the reference
    normal). Returns the output path.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "m3c2_alg_id", ""),
                      [["compare", "point"], ["m3c2"]])
    params = {"INPUT": reference_cloud, "INPUT_COMPARE": compare_cloud,
              "OUTPUT": out_cloud}
    if getattr(cfg, "m3c2_normal_scale", 0.0):
        params["NORMAL_RADIUS"] = cfg.m3c2_normal_scale
    if getattr(cfg, "m3c2_cyl_radius", 0.0):
        params["CYLINDER_RADIUS"] = cfg.m3c2_cyl_radius
    if getattr(cfg, "m3c2_subsample", 0.0):
        params["SUBSAMPLING_CELL_SIZE"] = cfg.m3c2_subsample
    if getattr(cfg, "m3c2_registration_error", 0.0):
        params["REGISTRATION_ERROR"] = cfg.m3c2_registration_error
    params.update(getattr(cfg, "m3c2_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               "M3C2 compare point clouds")
    return out_cloud if res is not None else None


def export_to_raster(cloud_path, out_tif, cfg, attribute="Z", filter_expression=None):
    """Rasterise a point-cloud attribute via ``pdal:exportraster``.

    ``attribute="Z"`` -> a **DTM**; ``attribute="m3c2_distance"`` -> the **change
    raster**. ``filter_expression`` (e.g. ``"Classification == 2"``) restricts to
    bare-earth points. Cell size = ``cfg.dtm_resolution``.
    """
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "export_raster_alg_id", ""),
                      [["export", "raster"]])
    params = {"INPUT": cloud_path, "OUTPUT": out_tif, "ATTRIBUTE": attribute}
    if getattr(cfg, "dtm_resolution", 0.0):
        params["RESOLUTION"] = cfg.dtm_resolution
    if filter_expression:
        params["FILTER_EXPRESSION"] = filter_expression
    params.update(getattr(cfg, "export_raster_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               f"export '{attribute}' to raster")
    return out_tif if res is not None else None


def export_dtm(cloud_path, out_tif, cfg, filter_expression=None):
    """Rasterise a bare-earth **DTM** from a cloud.

    ``cfg.dtm_method='tin'`` uses ``pdal:exportrastertin`` (Delaunay interpolation
    → gap-free, akin to IDFVA's las2ras interpolation); ``'binning'`` falls back to
    ``pdal:exportraster`` on the ``Z`` attribute (may leave nodata gaps). Both
    honour ``filter_expression`` (e.g. bare-earth ``"Classification == 2"``).
    """
    if getattr(cfg, "dtm_method", "tin").lower() != "tin":
        return export_to_raster(cloud_path, out_tif, cfg, attribute="Z",
                                filter_expression=filter_expression)
    processing = _processing()
    alg_id = _resolve(getattr(cfg, "exportrastertin_alg_id", ""),
                      [["export", "raster", "tin"], ["raster", "triangulation"]])
    params = {"INPUT": cloud_path, "OUTPUT": out_tif}
    if getattr(cfg, "dtm_resolution", 0.0):
        params["RESOLUTION"] = cfg.dtm_resolution
    if getattr(cfg, "dtm_max_edge_length", 0.0):
        params["MAX_EDGE_LENGTH"] = cfg.dtm_max_edge_length
    if filter_expression:
        params["FILTER_EXPRESSION"] = filter_expression
    params.update(getattr(cfg, "export_raster_params", {}) or {})
    res = _run(processing, alg_id, params, getattr(cfg, "pointcloud_dry_run", False),
               "export DTM (TIN interpolation)")
    return out_tif if res is not None else None


def prepare_from_clouds(pre_cloud, post_cloud, post_is_sfm, out_dir, cfg,
                        clip_overlay=None):
    r"""End-to-end Stage 0.

    ``clip_overlay`` (a polygon AOI, e.g. the watershed) clips both clouds first —
    strongly recommended for a fast first pass on huge SfM surveys.

    Returns ``dict(post_dtm=..., change_raster=..., m3c2_cloud=...)`` with
    ``change_raster``/``m3c2_cloud`` = ``None`` when there is no pre-event cloud
    (post-only ⇒ Q\* pathway only, no differencing/volume).
    """
    os.makedirs(out_dir, exist_ok=True)
    result = {"post_dtm": None, "change_raster": None, "m3c2_cloud": None}

    def _p(name):
        return os.path.join(out_dir, name)

    # 0a) clip both clouds to the AOI (optional; do this first to shrink data)
    pre_ready = pre_cloud
    post_ready = post_cloud
    if clip_overlay:
        if pre_cloud:
            pre_ready = clip_cloud(pre_cloud, clip_overlay, _p("pre_clip.laz"), cfg)
        post_ready = clip_cloud(post_cloud, clip_overlay, _p("post_clip.laz"), cfg)

    # 0b) reproject pre + post to a common CRS (optional; needed if they differ)
    epsg = getattr(cfg, "work_epsg", 0)
    if epsg:
        if pre_ready:
            pre_ready = reproject_cloud(pre_ready, _p("pre_reproj.laz"), epsg, cfg)
        post_ready = reproject_cloud(post_ready, _p("post_reproj.laz"), epsg, cfg)

    # 1) SfM post cloud: denoise, then classify ground
    if post_is_sfm and getattr(cfg, "denoise_sfm", True):
        post_ready = denoise_cloud(post_ready, _p("post_denoised.laz"), cfg)
    if post_is_sfm and getattr(cfg, "classify_sfm_ground", True):
        post_ready = classify_ground(post_ready, _p("post_ground.laz"), cfg)

    # bare-earth filter for the DTM (only when we have ground classes)
    dtm_filter = None
    if getattr(cfg, "dtm_ground_only", True):
        if post_is_sfm and not getattr(cfg, "classify_sfm_ground", True):
            print("  [pointcloud] dtm_ground_only set but SfM ground classification "
                  "is off — rasterising all points (set classify_sfm_ground=True).")
        else:
            dtm_filter = getattr(cfg, "ground_filter_expression", "Classification == 2")

    # post DTM (always — needed for Q*); TIN by default for a gap-free surface
    result["post_dtm"] = export_dtm(
        post_ready, _p("post_dtm.tif"), cfg, filter_expression=dtm_filter)

    if not pre_cloud:
        print("  [pointcloud] no pre-event cloud → differencing/volume SKIPPED "
              "(Q* pathway only).")
        return result

    # M3C2 change (reference = pre, comparison = post) -> change raster.
    # Output MUST be COPC/LAS 1.4: M3C2 adds 7 dimensions (m3c2_distance, ...) and
    # a LAS 1.2 writer (inherited from a 1.2 input) rejects them
    # ("LAS version 1.2 only supports point formats 0-5").
    m3c2_cloud = compare_m3c2(
        pre_ready, post_ready, _p("m3c2.copc.laz"), cfg)
    result["m3c2_cloud"] = m3c2_cloud
    if m3c2_cloud:
        result["change_raster"] = export_to_raster(
            m3c2_cloud, _p("change.tif"), cfg, attribute="m3c2_distance")
    return result
