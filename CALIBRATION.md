# Calibration record — bank-full Q\*

This file holds the site-specific validation behind `bankfull.py`. The module
itself is site-agnostic; the numbers below are what it produced on one particular
dataset and are recorded so the limits are auditable. **They are not defaults,
thresholds, or anything the code depends on.**

Test site: three burned catchments on the 2024 Bear Fire, Sierra County, CA.
Referred to below as basins A (a 0.20 km² catchment), B (0.07 km²) and C
(0.18 km²), each with field-mapped debris-flow trimlines and thalwegs.

## What was measured

Trimline-derived Q\* and bank-full-derived Q\* were computed on the *same*
transects, on the same bed DEM, at 5 m spacing. Any difference is therefore due
to the water-surface estimate alone.

### Bank-full over-predicts magnitude, on every dataset tried

| bed used | n | area × | Q\* × | width × | depth × |
|---|---|---|---|---|---|
| pre-event lidar | 155 | 61.4 | 175 | 4.81 | 9.5 |
| post-event, per basin | 154 | 26.3 | 54.5 | 4.71 | 5.21 |

An earlier, since-lost implementation of the same idea over-predicted by ~140×
(up to ~700× with a wider search) while agreeing **100 % on flow-type class**
(debris flow vs flood). That class agreement is the reason the module exists.

### The failure mode

In incised channels the bed rises *monotonically* out of the thalweg into the
valley wall. There is no local maximum, so crest detection never fires and every
transect terminates at the outward search cap. The over-prediction is then
essentially a function of that cap:

| `bankfull_max_halfwidth` | 6 m | 8 m | 10 m | 15 m | 25 m |
|---|---|---|---|---|---|
| median Q\* ratio | 9.6× | 17.4× | 23.5× | 36.7× | 53.6× |

Consequences, all measured:

* the `reliable` flag almost never fires — 5 % of transects at a 15 m cap, 0 %
  at 10 m, 18 % at 25 m. Too rare to screen on.
* `bankfull_min_rise` had **no effect at all** on this terrain (identical results
  from 0.10 to 1.00 m).

Tuning the cap until the ratio approaches 1 would fit an arbitrary parameter to
one site's answer and would not transfer to channels of a different size. The
module therefore reports when the cap was hit rather than hiding it.

### Hydraulic geometry does not rescue it

The proposed fix was to predict bank-full width and depth from upstream drainage
area (`A_us`) instead of searching the profile. Fitted across all three basins,
`y = a·A_us^b`, with `A_us` spanning 1.41 decades (0.020–0.503 km²):

| | R² (pre-event bed) | R² (post-event bed) | median err | constant-value baseline |
|---|---|---|---|---|
| width | 0.01 | 0.01 | 13.4 % | 17.7 % |
| depth | −0.06 | 0.13 | 29.2 % | 34.2 % |
| area | −0.12 | 0.09 | 40.6 % | 53.3 % |

On a correctly-vintaged bed the depth fit becomes positive, but it still barely
beats predicting the median. **Drainage area does not predict channel geometry at
this scale**, and depth is the term that drives `A_xs`.

### DEM vintage dominates — check this before trusting any of the above

Measuring identical transects on a pre-event and a post-event DEM changed depth
by 1.20×, 2.30× and 5.78× in basins A, B and C respectively. Width was unchanged
(it is set by the trimline stations), which confirms only the bed moved. A 6×
apparent between-basin depth spread collapsed to 1.6× once every basin was
measured on its own post-event surface.

If you compare basins measured on surfaces of different ages, that difference
will dominate anything the method reports.

## Bottom line

Use bank-full Q\* for **flow-type classification**. Do not quote its magnitude.
The `status`, `limiting_side` and `reliable` fields are the honest part of the
output — do not discard them.

Raw paired transects: `calibration_bankfull_vs_trimline.csv`.
