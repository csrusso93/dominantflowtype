# -*- coding: utf-8 -*-
r"""Regression test: filled depressions must not swallow flow accumulation.

Background
----------
`_priority_flood_fill` originally raised pits to *exactly* the spill elevation,
producing perfectly flat regions. `_d8_accumulation` assigns a receiver only
where the steepest slope is ``> 0``, so every cell in such a flat became a
terminal sink and all upstream flow stopped there.

Severity depends on accumulation resolution, because coarsening averages small
depressions away before they can become flats. On a real 0.445 km2 burned
catchment A_us reached only 17 % of the basin at 0.5 m and 58 % at 2 m, but was
unaffected at the 10 m default. Since ``Q_fluv = A_us * I``, a fine-resolution
run understates A_us and inflates Q\* by the same factor.

Run:  python tests_hydrology_flats.py
"""
import math
import numpy as np

from dominantflowtype.hydrology import _priority_flood_fill, _d8_accumulation

NB = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]


def terminal_cells(filled, valid, cell):
    """Valid interior cells with no downslope neighbour."""
    ny, nx = filled.shape
    z = np.where(valid, filled, np.inf)
    best = np.full((ny - 2, nx - 2), -np.inf)
    for di, dj in NB:
        d = math.hypot(di, dj) * cell
        with np.errstate(invalid="ignore"):
            s = (z[1:-1, 1:-1] - z[1 + di:ny - 1 + di, 1 + dj:nx - 1 + dj]) / d
        best = np.maximum(best, s)
    out = np.zeros((ny, nx), bool)
    out[1:-1, 1:-1] = best <= 0
    return out & valid


def _ramp_with_pit(n=60, cell=2.0, pit=(30, 30), r=6, depth=5.0):
    """A V-shaped valley draining to one edge, with a closed depression carved
    into the valley floor.

    The valley matters: on a tilted *plane* flow is parallel, nothing converges,
    and max accumulation is only ever one column -- so a "most of the grid
    reaches the outlet" assertion would be meaningless. Here every cell drains
    to the axis and then out the +y edge.
    """
    y, x = np.mgrid[0:n, 0:n]
    cx = (n - 1) / 2.0
    z = 100.0 - y * 0.20 + np.abs(x - cx) * 0.60  # cross-slope >> down-valley
    #                                              so flow converges to the axis
    yy, xx = np.ogrid[0:n, 0:n]
    m = (yy - pit[0]) ** 2 + (xx - pit[1]) ** 2 <= r * r
    z[m] -= depth                # carve a pit ON the valley axis
    return z, cell


def test_no_flats_after_fill():
    z, cell = _ramp_with_pit()
    valid = np.isfinite(z)
    filled = _priority_flood_fill(z)
    flats = terminal_cells(filled, valid, cell).sum()
    assert flats == 0, f"filled DEM still has {flats} terminal flat cells"
    print(f"  no-flats ................ OK (0 terminal cells)")


def test_fill_is_monotone_and_minimal():
    """Filling must never lower terrain, and must barely raise it."""
    z, cell = _ramp_with_pit()
    filled = _priority_flood_fill(z)
    assert np.all(filled >= z - 1e-9), "fill lowered the surface somewhere"
    raised = filled - z
    outside_pit = raised[raised < 1e-6]
    assert outside_pit.size > 0
    assert raised.max() <= 5.0 + 1e-6, "raised more than the pit depth"
    print(f"  monotone/minimal ........ OK (max raise {raised.max():.3f} m)")


def test_accumulation_reaches_outlet():
    """Flow from the whole grid must reach the outlet edge despite the pit."""
    z, cell = _ramp_with_pit()
    valid = np.isfinite(z)
    acc = _d8_accumulation(_priority_flood_fill(z), cell, valid)
    n = z.size
    frac = acc.max() / n
    # The ceiling is well below 100 % by construction, not because of flats:
    # _d8_accumulation assigns receivers only to interior cells ([1:-1, 1:-1]),
    # so the border ring never routes, and D8 converges imperfectly on a
    # synthetic valley. 0.70 sits comfortably above the ~0.58 seen when a pit
    # swallows flow, and below the ~0.76 this geometry actually achieves.
    assert frac > 0.70, (
        f"max accumulation is only {frac:.1%} of the grid -- flow is being "
        f"lost, most likely in a flat")
    print(f"  accumulation ............ OK ({frac:.1%} of grid reaches outlet)")


def test_pit_does_not_reduce_accumulation():
    """A pit must not REDUCE accumulation. (It may raise it: a depression is a
    convergence point, so once the epsilon fill lets flow out again it gathers
    more than the plain valley -- ~1.6x here. The failure mode being guarded
    against is a ratio well below 1.)"""
    z_pit, cell = _ramp_with_pit()
    z_flat, _ = _ramp_with_pit(depth=0.0)
    valid = np.isfinite(z_pit)
    a_pit = _d8_accumulation(_priority_flood_fill(z_pit), cell, valid).max()
    a_flat = _d8_accumulation(_priority_flood_fill(z_flat), cell, valid).max()
    ratio = a_pit / a_flat
    assert ratio > 0.95, (
        f"pit costs {(1-ratio):.1%} of accumulation (was ~42 % before the "
        f"epsilon fill)")
    print(f"  pit vs no-pit ........... OK (ratio {ratio:.3f})")


if __name__ == "__main__":
    print("hydrology flats regression")
    test_no_flats_after_fill()
    test_fill_is_monotone_and_minimal()
    test_accumulation_reaches_outlet()
    test_pit_does_not_reduce_accumulation()
    print("all passed")
