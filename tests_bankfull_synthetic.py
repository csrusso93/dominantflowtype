"""Synthetic tests for bankfull.py — channels with a KNOWN correct answer."""
import sys
import numpy as np

# run from the repo parent: python -m dominantflowtype.tests_bankfull_synthetic
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dominantflowtype.bankfull import bankfull_wse, calibrate_against_trimlines
from dominantflowtype.config import Config
from dominantflowtype import cross_section


class ProfileSampler:
    """Fake RasterSampler: bed elevation is a 1-D function of station u."""
    def __init__(self, fn):
        self.fn = fn

    def sample(self, px, py):
        return self.fn(px)          # transect laid along +x so px == station


TRANSECT = dict(x=0.0, y=0.0, lx=1.0, ly=0.0)   # unit left-normal along +x
cfg = Config()

print(f"cfg: bankfull_max_halfwidth={cfg.bankfull_max_halfwidth} m, "
      f"min_rise={cfg.bankfull_min_rise} m, transect_halfwidth={cfg.transect_halfwidth} m\n")

# ---------------------------------------------------------------- case 1
# Simple symmetric channel: 6 m wide, 1.5 m deep, banks then flat floodplain.
def chan_sym(u):
    a = abs(u)
    if a <= 3.0:
        return 1.5 * (a / 3.0) ** 2      # parabolic bed, 0 at centre -> 1.5 at bank
    return 1.5                           # flat floodplain beyond the crest

r = bankfull_wse(TRANSECT, ProfileSampler(chan_sym), cfg)
print("CASE 1  symmetric 6 m x 1.5 m channel, flat beyond banks")
print(f"  wse={r['wse']:.3f} (expect ~1.5)  depth={r['bankfull_depth']:.3f} (expect ~1.5)")
print(f"  width={r['width']:.2f} m (expect ~6)  reliable={r['reliable']}")
print(f"  status L/R = {r['status_left']}/{r['status_right']}   {r['note']}\n")

# ---------------------------------------------------------------- case 2
# Asymmetric: left bank crest at 1.0 m, right bank crest at 2.0 m.
# Cavagnaro rule -> WSE must be the LOWER crest (1.0), not the higher.
def chan_asym(u):
    if u >= 0:
        return min(1.0 * (u / 3.0), 1.0) if u <= 3.0 else 1.0
    a = -u
    return min(2.0 * (a / 3.0), 2.0) if a <= 3.0 else 2.0

r = bankfull_wse(TRANSECT, ProfileSampler(chan_asym), cfg)
print("CASE 2  asymmetric banks: left crest 1.0 m, right crest 2.0 m")
print(f"  wse={r['wse']:.3f}  -> EXPECT 1.0 (the LOWER bank, per Cavagnaro)")
print(f"  limiting_side={r['limiting_side']} (expect 'left')")
print(f"  VERDICT: {'PASS' if abs(r['wse'] - 1.0) < 0.1 else 'FAIL'}\n")

# ---------------------------------------------------------------- case 3
# THE FAILURE MODE: incised channel inside a wide V-shaped valley.
# True banks at +/-3 m (1 m deep). Valley walls climb forever beyond.
# Old behaviour: walk out to the valley wall -> huge section.
def chan_in_valley(u):
    a = abs(u)
    if a <= 3.0:
        return 1.0 * (a / 3.0)           # the real channel, 1 m deep
    return 1.0 + 0.30 * (a - 3.0)        # valley wall, climbs 0.3 m/m forever

r = bankfull_wse(TRANSECT, ProfileSampler(chan_in_valley), cfg)
print("CASE 3  1 m channel incised in a V-valley (the ~140x failure mode)")
print(f"  wse={r['wse']:.3f}  depth={r['bankfull_depth']:.3f}  width={r['width']:.2f} m")
print(f"  reliable={r['reliable']}  status L/R = {r['status_left']}/{r['status_right']}")
print(f"  note: {r['note']}")
capped = "capped" in (r["status_left"], r["status_right"])
print(f"  VERDICT: {'PASS - flagged, not silently inflated' if capped and not r['reliable'] else 'FAIL - not flagged'}")
print(f"  (search capped at {cfg.bankfull_max_halfwidth} m => depth bounded ~"
      f"{1.0 + 0.30 * (cfg.bankfull_max_halfwidth - 3.0):.1f} m, vs "
      f"{1.0 + 0.30 * (cfg.transect_halfwidth - 3.0):.1f} m if it ran the full transect)\n")

# how much does the cap actually save?
cfg_uncapped = Config()
cfg_uncapped.bankfull_max_halfwidth = 50.0
r2 = bankfull_wse(TRANSECT, ProfileSampler(chan_in_valley), cfg_uncapped)
print(f"  same profile, cap raised to 50 m: depth={r2['bankfull_depth']:.2f} m, "
      f"width={r2['width']:.1f} m")
print(f"  -> capping cuts the section depth {r2['bankfull_depth'] / r['bankfull_depth']:.1f}x\n")

# ---------------------------------------------------------------- case 4
# End-to-end: feed the estimated WSE into measure_cross_section.
r = bankfull_wse(TRANSECT, ProfileSampler(chan_sym), cfg)
geom = cross_section.measure_cross_section(
    TRANSECT, r["wse"], None, None, ProfileSampler(chan_sym), cfg)
print("CASE 4  hand the bank-full WSE to cross_section.measure_cross_section")
print(f"  A_xs={geom['A_xs']:.3f} m^2   R={geom['hydraulic_radius']:.3f} m   "
      f"top_width={geom['top_width']:.2f} m")
# analytic area for parabola z=1.5(a/3)^2 under wse=1.5 over [-3,3] = 2/3*b*h = 6 m^2
print(f"  analytic A_xs for this parabola = {2/3*6*1.5:.3f} m^2  "
      f"-> {'PASS' if abs(geom['A_xs'] - 6.0) < 0.3 else 'FAIL'}\n")

# ---------------------------------------------------------------- case 5
print("CASE 5  calibrate_against_trimlines (synthetic 3x over-prediction)")
paired = [(2.0, 6.0), (3.0, 9.5), (1.5, 4.2), (4.0, 12.5), (2.5, 60.0)]  # last = valley wall
c = calibrate_against_trimlines(paired)
print(f"  n={c['n']}  median ratio={c['ratio']:.2f} (expect ~3, robust to the 24x outlier)")
print(f"  area_factor={c['area_factor']:.3f}  range {c['ratio_min']:.1f}-{c['ratio_max']:.1f}")
print(f"  VERDICT: {'PASS' if 2.5 < c['ratio'] < 3.5 else 'FAIL'}")
