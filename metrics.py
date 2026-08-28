# -*- coding: utf-8 -*-
"""Velocity, discharge, Q* and cross-section usability (Cavagnaro et al., 2024)."""
from __future__ import annotations

import math

import numpy as np


def velocity(hydraulic_radius, hydraulic_depth, max_depth, cfg):
    """Froude-critical velocity v = Fr*sqrt(g*h) (Cavagnaro et al., 2024)."""
    if cfg.velocity_scale == "hydraulic_depth":
        h = hydraulic_depth
    elif cfg.velocity_scale == "max_depth":
        h = max_depth
    else:
        h = hydraulic_radius
    if not np.isfinite(h) or h <= 0:
        return float("nan")
    return cfg.froude * math.sqrt(cfg.g * h)


def compute_qstar(A_xs, wetted_perimeter, hydraulic_radius, hydraulic_depth,
                  max_depth, A_us, I_mm_hr, cfg):
    """Return dict of v, Q_peak, Q_fluv, Q*, flow_type for one cross-section/bed."""
    v = velocity(hydraulic_radius, hydraulic_depth, max_depth, cfg)
    Q_peak = v * A_xs if np.isfinite(v) else float("nan")
    I_ms = (I_mm_hr / 1000.0) / 3600.0             # mm/hr -> m/s
    Q_fluv = A_us * I_ms if (A_us and np.isfinite(A_us)) else float("nan")
    Qstar = (Q_peak / Q_fluv) if (Q_fluv and Q_fluv > 0) else float("nan")
    if np.isfinite(Qstar):
        ftype = "debris flow" if Qstar > cfg.qstar_threshold else "flood"
    else:
        ftype = "undetermined"
    return dict(velocity=v, Q_peak=Q_peak, Q_fluv=Q_fluv, Qstar=Qstar,
                flow_type=ftype)


def usability(depth, dz_bed, continuity, cfg):
    """Flag whether a cross-section is usable per Cavagnaro Fig 1d/1e.

    Unusable ONLY when incision/deposition is large relative to flow depth.
    Trimline discontinuity does NOT disqualify a section (constant-depth
    inference is applied instead).
    """
    reasons = []
    if not np.isfinite(depth) or depth < cfg.min_depth:
        reasons.append("depth below minimum / no valid section")
    if dz_bed is not None and np.isfinite(dz_bed) and np.isfinite(depth) and depth > 0:
        ratio = abs(dz_bed) / depth
        if ratio > cfg.max_dz_to_depth_ratio:
            kind = "incision" if dz_bed < 0 else "deposition"
            reasons.append(f"heavy {kind} (|dz|/depth={ratio:.2f} "
                           f">{cfg.max_dz_to_depth_ratio})")
    usable = len(reasons) == 0
    return usable, ("; ".join(reasons) if reasons else "usable"), continuity
