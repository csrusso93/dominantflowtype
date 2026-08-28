# -*- coding: utf-8 -*-
"""Inundated cross-section geometry (A_xs, wetted perimeter, R, depth, width)."""
from __future__ import annotations

import math

import numpy as np

from ._compat import _trapz


def measure_cross_section(transect, wse, left_station, right_station,
                          bed_sampler, cfg):
    """Compute geometry of the inundated cross-section on ONE bed surface.

    Parameters
    ----------
    wse : float
        Water-surface (peak-flow) elevation for this transect.
    left_station, right_station : float
        Signed along-transect distances (+ toward left endpoint) of the
        inundation boundary on each bank. If None on a side, the boundary is
        found where the bed rises to `wse` (inferred / constant-depth edge).
    bed_sampler : RasterSampler
        DEM used as the channel bed (pre- or post-event).

    Returns dict with A_xs, wetted_perimeter, hydraulic_radius, top_width,
    max_depth, and the resolved left/right boundary stations & xy.
    """
    st = (transect["x"], transect["y"])
    lx, ly = transect["lx"], transect["ly"]     # unit left-normal

    # Dense bed profile across the full transect, parameterised by signed
    # station u (+ = toward left endpoint).
    L = cfg.transect_halfwidth
    n = max(int(math.ceil(2 * L / cfg.profile_step)) + 1, 5)
    u = np.linspace(-L, L, n)                    # +L is left endpoint
    xs = st[0] + lx * u
    ys = st[1] + ly * u
    z = np.array([bed_sampler.sample(px, py) for px, py in zip(xs, ys)])

    # index of channel bed (deepest point near centre)
    center_mask = np.abs(u) <= max(cfg.transect_spacing, 5.0)
    if center_mask.any() and np.isfinite(z[center_mask]).any():
        cz = np.where(center_mask & np.isfinite(z), z, np.inf)
        i_bed = int(np.argmin(cz))
    else:
        i_bed = int(np.nanargmin(z)) if np.isfinite(z).any() else n // 2

    def _edge(direction):
        """Walk outward from bed until bed >= wse or a mapped station; return u."""
        i = i_bed
        while 0 < i < n - 1:
            i += direction
            if not np.isfinite(z[i]):
                break
            if z[i] >= wse:
                # linear crossing between i-direction and i
                z0, z1 = z[i - direction], z[i]
                u0, u1 = u[i - direction], u[i]
                if np.isfinite(z0) and z1 != z0:
                    t = (wse - z0) / (z1 - z0)
                    return u0 + t * (u1 - u0)
                return u[i]
        return u[i]

    # resolve boundary stations: prefer mapped trimline stations
    uL = left_station if left_station is not None else _edge(+1)
    uR = right_station if right_station is not None else _edge(-1)
    if uL is None:
        uL = _edge(+1)
    if uR is None:
        uR = _edge(-1)
    if uL < uR:
        uL, uR = uR, uL                          # ensure uL (left) > uR (right)

    # integrate depth over [uR, uL]
    inside = (u >= uR) & (u <= uL)
    if inside.sum() < 2:
        return None
    ui = u[inside]
    zi = z[inside]
    depth = wse - zi
    depth = np.where(np.isfinite(depth), depth, 0.0)
    depth = np.clip(depth, 0.0, None)

    A_xs = float(_trapz(depth, ui))
    top_width = float(uL - uR)
    max_depth = float(np.nanmax(depth)) if np.isfinite(depth).any() else 0.0

    # wetted perimeter: bed length where submerged
    du = np.diff(ui)
    dz = np.diff(np.where(np.isfinite(zi), zi, wse))
    seg = np.hypot(du, dz)
    submerged = (depth[:-1] > 0) | (depth[1:] > 0)
    wetted_perimeter = float(np.sum(seg[submerged]))
    if wetted_perimeter <= 0:
        wetted_perimeter = top_width

    hydraulic_radius = A_xs / wetted_perimeter if wetted_perimeter > 0 else 0.0
    hydraulic_depth = A_xs / top_width if top_width > 0 else 0.0

    return dict(
        A_xs=A_xs, wetted_perimeter=wetted_perimeter,
        hydraulic_radius=hydraulic_radius, hydraulic_depth=hydraulic_depth,
        top_width=top_width, max_depth=max_depth,
        uL=uL, uR=uR,
        left_xy=(st[0] + lx * uL, st[1] + ly * uL),
        right_xy=(st[0] + lx * uR, st[1] + ly * uR),
        bed_min=float(np.nanmin(zi)) if np.isfinite(zi).any() else float("nan"),
    )
