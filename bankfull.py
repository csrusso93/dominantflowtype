# -*- coding: utf-8 -*-
r"""Bank-full water-surface estimation for channels with no mapped trimline.

Why this exists
---------------
:mod:`cross_section` needs a water-surface elevation (WSE) for every transect.
Where a debris flow has been field-mapped, the WSE comes from the **trimline**
and the resulting Q\* is quantitative (Cavagnaro et al., 2024). Most channels
have no trimline. This module estimates a **bank-full** WSE from the bed profile
alone, so Q\* can be run as a *susceptibility screen* on un-mapped channels.

**Read this before using the output.** A previous implementation of this idea
(since lost; described in the project's ``METHODS_SUMMARY``) was validated at
BC-E10 against trimline-based Q\* on the same transects:

* **flow-type classification agreed 100 %** (debris flow vs flood), but
* it **over-predicted Q\* magnitude by ~140x** (up to ~700x with a wider bank
  search), because at this channel scale the outward search ran past the real
  banks and locked onto **valley walls** -- the resulting section was ~4x wider
  and ~10x deeper than the mapped debris flow.

So: bank-full Q\* is a **class proxy, not a quantitative Q\***. This rewrite
attacks the specific failure mode above by (a) capping the outward search at a
channel scale rather than the full transect half-width, (b) detecting the bank
*crest* by slope reversal instead of walking until the bed rises to some level,
and (c) **reporting when the cap was hit**, so valley-scale sections can be
filtered out instead of silently inflating Q\*. Those flags are the point --
do not discard them.

The prior 100 %/140x validation applies to the *lost* implementation. This one
**needs re-validating** against the BC-E10 and SNC-W8 trimlined transects
before its magnitudes are trusted; :func:`calibrate_against_trimlines` exists
for exactly that.

Method
------
For each transect, on the bed DEM:

1. Sample the bed profile (same parameterisation as :mod:`cross_section`).
2. Locate the channel bed: deepest point near the transect centre.
3. Walk outward on each side to the **first bank crest** -- the first local
   maximum in the profile, i.e. where the bed stops rising and turns over --
   subject to ``bankfull_max_halfwidth`` and a minimum rise of
   ``bankfull_min_rise``.
4. **WSE = the lower of the two crest elevations.** This follows Cavagnaro et
   al. (2024) §2.3, who use the lower of two flow-depth indicators when the
   banks differ, and is conservative: water would spill at the lower bank.

Returns the WSE plus diagnostics; feed the WSE to
:func:`cross_section.measure_cross_section` with ``left_station=None,
right_station=None`` so the section is closed at that elevation.
"""
from __future__ import annotations

import math

import numpy as np


def _bed_profile(transect, bed_sampler, cfg):
    """Dense bed profile across a transect. Mirrors cross_section.measure_cross_section."""
    st = (transect["x"], transect["y"])
    lx, ly = transect["lx"], transect["ly"]
    L = cfg.transect_halfwidth
    n = max(int(math.ceil(2 * L / cfg.profile_step)) + 1, 5)
    u = np.linspace(-L, L, n)
    xs = st[0] + lx * u
    ys = st[1] + ly * u
    z = np.array([bed_sampler.sample(px, py) for px, py in zip(xs, ys)])
    return u, z


def _bed_index(u, z, cfg):
    """Index of the channel bed: deepest finite point near the transect centre."""
    n = len(u)
    center_mask = np.abs(u) <= max(cfg.transect_spacing, 5.0)
    if center_mask.any() and np.isfinite(z[center_mask]).any():
        cz = np.where(center_mask & np.isfinite(z), z, np.inf)
        return int(np.argmin(cz))
    return int(np.nanargmin(z)) if np.isfinite(z).any() else n // 2


def _find_crest(u, z, i_bed, direction, cfg):
    """Walk outward from the bed to the first bank crest on one side.

    A crest is the first point that is a local maximum -- the bed rose by at
    least ``bankfull_min_rise`` above the thalweg and then turned over (the next
    ``bankfull_turnover_pts`` samples do not exceed it).

    Returns (u_crest, z_crest, status) where status is one of:
      'crest'   - a genuine turnover was found (trustworthy)
      'capped'  - hit bankfull_max_halfwidth still rising (VALLEY-WALL RISK:
                  the section is bounded by the cap, not by a real bank)
      'nodata'  - ran into missing DEM data
      'flat'    - never rose bankfull_min_rise above the bed
    """
    n = len(u)
    z_bed = z[i_bed]
    max_hw = getattr(cfg, "bankfull_max_halfwidth", 15.0)
    min_rise = getattr(cfg, "bankfull_min_rise", 0.25)
    turn_pts = int(getattr(cfg, "bankfull_turnover_pts", 3))

    best_i = None
    i = i_bed
    while 0 < i < n - 1:
        i += direction
        if not np.isfinite(z[i]):
            return (u[i - direction], z[i - direction], "nodata")
        if abs(u[i] - u[i_bed]) > max_hw:
            # Ran out of allowed channel width while still climbing.
            return (u[i], z[i], "capped")
        if z[i] - z_bed < min_rise:
            continue                                  # still in the channel floor
        # local maximum test: no following sample within turn_pts exceeds this one
        j0 = i + direction
        j1 = i + direction * (turn_pts + 1)
        lo, hi = sorted((j0, j1))
        lo = max(lo, 0)
        hi = min(hi, n)
        ahead = z[lo:hi]
        ahead = ahead[np.isfinite(ahead)]
        if ahead.size == 0:
            return (u[i], z[i], "nodata")
        if np.nanmax(ahead) <= z[i]:
            return (u[i], z[i], "crest")
        best_i = i
    if best_i is not None:
        return (u[best_i], z[best_i], "capped")
    return (u[i], z[i], "flat")


def bankfull_wse(transect, bed_sampler, cfg):
    """Estimate the bank-full WSE for one transect.

    Returns
    -------
    dict or None
        ``wse``            estimated bank-full water-surface elevation [m]
        ``z_bed``          thalweg elevation [m]
        ``bankfull_depth`` wse - z_bed [m]
        ``u_left``/``u_right``, ``z_left``/``z_right``  crest stations/elevations
        ``status_left``/``status_right``  per-side crest status (see _find_crest)
        ``limiting_side``  which bank set the WSE ('left' or 'right')
        ``reliable``       True only if BOTH sides found a real crest
        ``width``          crest-to-crest top width [m]
        ``note``           human-readable caveat
        ``None`` if no usable bed profile.
    """
    u, z = _bed_profile(transect, bed_sampler, cfg)
    if not np.isfinite(z).any():
        return None
    i_bed = _bed_index(u, z, cfg)
    z_bed = float(z[i_bed])
    if not np.isfinite(z_bed):
        return None

    uL, zL, sL = _find_crest(u, z, i_bed, +1, cfg)     # toward left endpoint
    uR, zR, sR = _find_crest(u, z, i_bed, -1, cfg)     # toward right endpoint
    if not (np.isfinite(zL) and np.isfinite(zR)):
        return None

    # Cavagnaro et al. (2024): where the two banks differ, use the LOWER.
    if zL <= zR:
        wse, limiting = float(zL), "left"
    else:
        wse, limiting = float(zR), "right"

    reliable = (sL == "crest" and sR == "crest")
    depth = wse - z_bed
    if depth <= 0:
        return None

    notes = []
    if not reliable:
        notes.append(
            f"bank crest not resolved on both sides (left={sL}, right={sR}); "
            "section may be valley-scale -- treat as class proxy only")
    if "capped" in (sL, sR):
        notes.append(
            f"outward search hit bankfull_max_halfwidth="
            f"{getattr(cfg, 'bankfull_max_halfwidth', 15.0)} m -- the known "
            "valley-wall failure mode; Q* magnitude will be over-predicted")

    return dict(
        wse=wse, z_bed=z_bed, bankfull_depth=float(depth),
        u_left=float(uL), u_right=float(uR),
        z_left=float(zL), z_right=float(zR),
        status_left=sL, status_right=sR,
        limiting_side=limiting, reliable=bool(reliable),
        width=float(uL - uR),
        note="; ".join(notes) if notes else "both bank crests resolved",
    )


def calibrate_against_trimlines(paired, robust=True):
    """Compare bank-full sections against mapped-trimline sections.

    The prior implementation over-predicted Q\\* magnitude ~140x. Rather than
    assume a correction, derive one where BOTH a trimline and a bank-full
    section exist for the same transect (BC-E10, SNC-W8).

    Parameters
    ----------
    paired : sequence of (A_xs_trimline, A_xs_bankfull)
        Inundated cross-sectional areas, same transects, same bed DEM.
    robust : bool
        Use the median ratio (default) rather than the mean; area ratios are
        heavy-tailed and a handful of valley-wall sections would dominate a mean.

    Returns
    -------
    dict with ``n``, ``ratio`` (bankfull/trimline area), ``area_factor``
    (multiply a bank-full A_xs by this to approximate a trimline A_xs),
    ``qstar_factor`` (the implied Q\\* over-prediction, since Q\\* is linear in
    A_xs at fixed velocity), and the ratio spread.

    Note Q\\* is *not* strictly linear in A_xs: velocity scales with the
    hydraulic radius, so a corrected area implies a corrected velocity too.
    Use this to judge the size of the bias, not as a blind multiplier.
    """
    a = np.asarray([p[0] for p in paired], dtype="float64")
    b = np.asarray([p[1] for p in paired], dtype="float64")
    ok = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    a, b = a[ok], b[ok]
    if a.size == 0:
        return dict(n=0, ratio=float("nan"), area_factor=float("nan"),
                    qstar_factor=float("nan"), ratio_iqr=float("nan"))
    r = b / a
    ratio = float(np.median(r) if robust else np.mean(r))
    q1, q3 = np.percentile(r, [25, 75])
    return dict(
        n=int(a.size),
        ratio=ratio,
        area_factor=float(1.0 / ratio) if ratio > 0 else float("nan"),
        qstar_factor=ratio,
        ratio_iqr=float(q3 - q1),
        ratio_min=float(np.min(r)), ratio_max=float(np.max(r)),
    )
