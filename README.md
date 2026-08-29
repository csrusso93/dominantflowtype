# dominantflowtype

Post-event, DEM-based classification of the **dominant flow type** — runoff-generated
**debris flow** vs. **flood** — in steep, burned channels, implementing the dimensionless
discharge metric **Q\*** of Cavagnaro et al. (2024, *GRL*, doi:10.1029/2024GL109768):

```
Q* = Q_peak / Q_fluv = (v · A_xs) / (I · A_us)      Q* > 1 → debris flow ; Q* ≤ 1 → flood
```

| term | meaning | source |
|------|---------|--------|
| `A_xs` | inundated cross-sectional area | trimlines + DEM bed |
| `P`, `R` | wetted perimeter, hydraulic radius `R = A_xs/P` | DEM bed |
| `v = √(g·R)` | Froude-critical velocity (Fr = 1) | Cavagnaro §2.3 |
| `A_us` | upstream drainage area | D8 flow routing on the DEM, clipped to the watershed |
| `I` | 30-min rainfall intensity `I30` | **MRMS radar** *(no account)*, **Synoptic** gauges, or 9.7 mm/hr |

### Verified against the source paper (2026-08-28)

The implementation was checked line-by-line against Cavagnaro et al. (2024) and is
faithful on every parameter that matters:

* **`h` is the hydraulic radius** `A_xs / wetted_perimeter`, *not* flow depth —
  `metrics.velocity` defaults to `velocity_scale="hydraulic_radius"`, matching the
  paper's `h`. (`hydraulic_depth` and `max_depth` remain available as explicit
  opt-ins for sensitivity testing.)
* `Fr = 1.0`, `g = 9.81`, `qstar_threshold = 1.0` — paper §2.3.
* `I` is the peak **30-minute** intensity — the paper's stated compromise between
  the concentration time of small headwater basins and rainfall-data availability.
* The **9.7 mm/hr** fallback is the paper's own optimum when `I30` is held constant
  across all observations (84.1 % accuracy, vs. 91.5 % event-specific).
* Where both banks are mapped, the water-surface elevation is the **lower of the
  two**, per the paper.
* Cross-sections are rejected where incision or deposition is large relative to
  flow depth, per the paper's Figure 1d/1e (`metrics.usability`).

For reference, Q\* separates 642 observed flows (270 floods, 372 debris flows) to
**91.5 ± 0.9 %** accuracy, with mean Q\* ~42× higher for debris flows than floods
(~63× post-fire).

**Scope note.** This package computes **Q\***. The rainfall anomaly **R\***
(Cavagnaro et al., 2025) is *not* implemented here — it is produced by the companion
[`stormscape`](https://github.com/scottwmccoy/stormscape) tool. Likewise, the M3C2
step here yields a **change raster**, not incremental along-path volumes; see
Scheip & Wegmann (2022) and Guido (2025) under Reference.

## Data import (methods ported from `stormscape`)

The DEM and rainfall importers follow **[`stormscape`](https://github.com/scottwmccoy/stormscape)**
(S. W. McCoy, MIT), adapted to this package's GDAL/QGIS stack (no
rasterio/xarray). All sources are account-free public US data except Synoptic:

| module | what it imports | source | account |
|--------|-----------------|--------|---------|
| `mrms.py` | peak i15/i30/i60 rainfall-intensity fields + scalar `mrms_i30()` | NOAA **MRMS** radar (public S3, GDAL GRIB) | **none** |
| `rainfall.py` | gauge `I30` (derived precip, 1-min interpolation, `(i16+i14)/2` estimator) | **Synoptic / MesoWest** | token |
| `dem_3dep.py` | bare-earth DEM for an AOI at 1/3/5/10/30 m | USGS **3DEP** (py3dep *or* the GDAL ImageServer fallback) | **none** |
| `aoi.py` | normalise a bbox / vector / QGIS layer to a WGS84 bounding box | — | — |

> **MRMS solves the Bear Fire rainfall blocker.** The Synoptic account tied to
> this project returns HTTP 403 (no data access). Set `rainfall_source="mrms"`
> and give a storm date to get a gauge-free `I30` straight from public radar —
> no account, no extra installs (QGIS 4.x GDAL already has the GRIB driver).

## Requirements

Runs inside a **QGIS Python environment** (uses `qgis.core` + GDAL). QGIS 3.34+ /
4.x ship every dependency it needs (`numpy`, `scipy`, `pandas`, `openpyxl`,
`requests`, GDAL). MRMS radar import additionally needs GDAL's **GRIB** driver —
present in the QGIS 4.x GDAL build (verified: GDAL 3.13, `GRIB` driver `YES`). The
faithful py3dep 3DEP path is optional (`pip install py3dep rioxarray`); without it
`dem_3dep` falls back to the account-free GDAL ImageServer path automatically.

## Run it in the QGIS Python Console

Open **QGIS → Plugins → Python Console**, then:

The line you add to `sys.path` must be the folder that **contains**
`dominantflowtype` — *not* the package folder itself. Set `ROOT` to wherever you
cloned the repository:

```python
import sys
ROOT = r"/path/to/parent/of/dominantflowtype"   # the folder CONTAINING the package
sys.path.insert(0, ROOT)

import dominantflowtype as dft
records, reaches, rain = dft.run()      # capture the return so you can inspect it
```

Alternatively, run `__main__.py` directly and skip the `sys.path` line entirely —
it locates the package itself, from its own `__file__`, from the
`DOMINANTFLOWTYPE_ROOT` environment variable, or by walking up from the working
directory:

```python
exec(open(r"/path/to/dominantflowtype/__main__.py", encoding="utf-8").read())
```

From a shell, `cd` anywhere at or below the folder containing the package and run:

```bash
python -m dominantflowtype
```

`run()` returns `(records, reaches, rain)`:

* **`records`** — one dict per transect (`Qstar`, `flow_type`, `depth_m`,
  `width_m`, `A_us`, `velocity`, ...);
* **`reaches`** — inundation polygons between transects (`area`, `mean_dz`,
  `volume`);
* **`rain`** — the I30 metadata (`i30_mm_hr`, `source`, `reduce`, `peak_utc`, ...).

Print a quick summary after it finishes:

```python
import numpy as np
qs = np.array([r["Qstar"] for r in records if np.isfinite(r.get("Qstar", np.nan))])
df = sum(r.get("flow_type") == "debris flow" for r in records)
fl = sum(r.get("flow_type") == "flood"       for r in records)
print(f"I30            = {rain['i30_mm_hr']:.1f} mm/hr ({rain.get('reduce','')})")
print(f"transects      = {len(records)}")
print(f"Q* min/med/max = {qs.min():.2f} / {np.median(qs):.2f} / {qs.max():.2f}")
print(f"classification = {df} debris flow / {fl} flood / {len(records)-df-fl} undetermined")
```

`run()` prompts you, one field at a time, for:

1. *(optional, skippable)* **Download a DEM** (USGS 3DEP, no key / OpenTopography).
   Answer **no** to bypass if you already have both DEMs. If yes, choose what it
   is *for*: `pre` (use it as the pre-event DEM), `routing` (full-coverage A_us
   routing DEM only — pre/post untouched), or `save` (just download a copy).
2. **Pre-event DEM** (type `skip` to omit; auto-filled if [1] purpose was `pre`)
3. **Post-event DEM** (type `skip` to omit)
4. **Trimlines** (`.shp` / `.gpkg`)
5. **Thalweg** (`.shp` / `.gpkg`)
6. **Watershed** outline (`.shp` / `.gpkg`)
7. **Rainfall source** — `mrms` (radar; asks for a storm date only), `synoptic`
   (asks for an API key + optional station id + date), or `constant`

> If **one** of the two DEMs is skipped, volume (DEM-of-difference) products are
> skipped automatically and only Q\* is computed on the single DEM.

### Point-cloud pathway (Stage 0) + IDFVA volume bridge

`run()` now offers an optional **point-cloud front-end** ([0] in the prompts):
supply pre/post LiDAR or SfM clouds and it runs **QGIS 4.0 native M3C2**
(*Compare Point Clouds*; Lague et al., 2013) → rasterises the post cloud to a
**DTM** for Q\* and the `m3c2_distance` field to a signed **change raster**. SfM
post clouds get a **CSF ground-filter** first. **Post-cloud-only ⇒ Q\* only**
(nothing to difference; volume is skipped).

Volume itself is **not** computed here — it is delegated to Guido's **IDFVA**
(separate `pip install idfva` venv). When a change surface exists (cloud M3C2 or a
two-DEM DoD), the run writes a **bridge bundle** (`bridge_bundle/`) with
`flow_path.shp` (thalweg), `post_dtm.tif`, `change.tif`, and `manifest.json` for
IDFVA to consume. See `../INTEGRATION_PLAN.md`.

The Stage-0 algorithm ids/params are auto-discovered from the QGIS Processing
registry. **Run this once in the QGIS console to print the exact ids on your
build** (then pin them via `Config`, e.g. `m3c2_alg_id=...`):

```python
import dominantflowtype as dft
dft.pointcloud_diagnose()      # lists point-cloud algorithms + their parameters
```

### Smoke test (fast end-to-end validation)

Before a full-resolution run, validate the whole chain in minutes on a small clip:

```python
import dominantflowtype as dft
# one-time: a 200 m AOI centred in the cloud (BC-E10 extent shown; use yours)
dft.make_test_aoi((738405.27, 4387968.70, 739533.81, 4388812.37),
                  "smoketest_aoi.gpkg", epsg=6339, size_m=200)
dft.run(cfg=dft.smoke_config())        # coarse subsample/DTM; actually runs
```
At the prompts pick the SfM cloud, `sfm`, the 3DEP DEM `[PC-4]`, and the
`smoketest_aoi.gpkg` at `[PC-5]`. `smoke_config()` uses a coarse M3C2 subsample and
2 m TIN DTM so it finishes quickly — use it to confirm the **M3C2 sign**
(deposition positive), CRS handling and the bridge, then rerun with a full
`Config(work_epsg=6339, dtm_resolution=0.5, m3c2_subsample=0.5, …)`.

### Scripting the importers directly

```python
import dominantflowtype as dft
# gauge-free I30 over a watershed for one storm day (returns mm/hr + meta)
i30, meta = dft.mrms_i30("watershed.gpkg", "2024-02-04", reduce="areal_max")
# a full-coverage 10 m routing DEM from USGS 3DEP (path to a GeoTIFF)
dem = dft.get_dem_3dep("watershed.gpkg", resolution=10, out_path="routing_10m.tif")
```

### Re-running after edits

The QGIS console caches imported modules. **`importlib.reload()` is not enough** —
it re-executes only the top package, leaving the edited *submodules*
(`rainfall`, `pipeline`, ...) stale, which surfaces as confusing `ImportError`s or
the old prompt order. Either **restart the console** (trash/restart icon in the
Python Console toolbar), or purge every submodule from `sys.modules` first:

```python
import sys
for _m in [m for m in list(sys.modules)
           if m == "dominantflowtype" or m.startswith("dominantflowtype.")]:
    del sys.modules[_m]
import dominantflowtype as dft
records, reaches, rain = dft.run()
```

You are on the fresh code when prompt **[1]** reads *"(optional) Download a DEM
(USGS 3DEP / OpenTopography)?"*.

## Outputs

* `<basin>_dominantflowtype.xlsx` — a single workbook with sheets **pre_event**,
  **post_event** (adds `dz_bed_m`, `volume_m3`), **usable_only** and
  **rainfall_I30**; unusable cross-sections flagged with a reason. (Per-sheet CSVs
  are written only as a fallback if the workbook write fails.)
* `DoD_post_minus_pre.tif` — DEM-of-Difference (post − pre).
* QGIS layers: **DFT_inundation_Qstar** (green Q\*<1 → white Q\*=1 → purple
  Q\*>1), **DFT_transects**, and **DFT_erosion_deposition** (red erosion → white →
  blue deposition). GeoPackage copies are saved to the output directory.

## Configuration

All parameters live in `dominantflowtype.Config` and can be overridden:

```python
from dominantflowtype import Config, run
cfg = Config(
    transect_spacing=5.0,          # transect interval [m]
    velocity_scale="hydraulic_radius",  # or "hydraulic_depth" / "max_depth"
    max_dz_to_depth_ratio=0.5,     # incision/deposition usability threshold
    flow_accum_res=10.0,           # DEM resample for A_us routing [m]
    rainfall_source="mrms",        # 'mrms' | 'synoptic' | 'constant'
    rain_event_date="2024-02-04",  # storm day for MRMS/Synoptic
    mrms_reduce="areal_max",       # 'areal_max' | 'areal_mean' | 'point'
    dem3dep_resolution=10,         # USGS 3DEP routing DEM [m]
    pre_vertical_unit="auto",      # 'm' | 'ft' | 'auto' (feet -> metres for DoD)
    post_vertical_unit="auto",     # 'm' | 'ft' | 'auto'
)
run(cfg=cfg)
```

## Scripted / non-interactive use

Pass a list of answers (same order as the prompts) to replay a session — used by
the test harness:

```python
from dominantflowtype import run
# order: [1] download? ("no" bypasses), preDEM, postDEM, trimlines, thalweg,
#        watershed, out_dir, rainfall_source [, source-specific answers]
records, reaches, rain = run(answers=[
    "no", r"pre.tif", r"post.tif", r"trim.gpkg", r"thal.shp",
    r"ws.shp", r"C:\out", "mrms", "2024-02-04", "areal_max"])

# with the [1] 3DEP download enabled (answer "yes" -> source, purpose, resolution):
# records, reaches, rain = run(answers=[
#     "yes", "3dep", "routing", "10", r"pre.tif", r"post.tif", r"trim.gpkg",
#     r"thal.shp", r"ws.shp", r"C:\out", "mrms", "2025-08-27", "point"])
```

## Install as a package (optional)

From the folder containing `pyproject.toml`:

```bash
pip install -e .
```

(GDAL and QGIS are provided by the QGIS environment, not by pip.)

## Reference

**Implemented here**

Cavagnaro, D. B., McCoy, S. W., Kean, J. W., Thomas, M. A., Lindsay, D. N.,
McArdell, B. W., & Hirschberg, J. (2024). *A Robust Quantitative Method to
Distinguish Runoff-Generated Debris Flows From Floods.* Geophysical Research
Letters, 51, e2024GL109768. doi:10.1029/2024GL109768 — **the Q\* metric.**

Lague, D., Brodu, N., & Leroux, J. (2013). *Accurate 3D comparison of complex
topography with terrestrial laser scanner: application to the Rangitikei canyon
(N-Z).* ISPRS J. Photogramm. Remote Sens., 82, 10–26.
doi:10.1016/j.isprsjprs.2013.04.009 — **M3C2**, used by `pointcloud.py`.

O'Callaghan, J. F., & Mark, D. M. (1984); Barnes, R., Lehman, C., & Mulla, D.
(2014) — **D8 flow routing** and **priority-flood** pit filling, used by
`hydrology.py`.

**Related, not implemented here**

Cavagnaro, D. B., McCoy, S. W., Thomas, M. A., Kostelnik, J., & Lindsay, D. N.
(2025). *Improved Prediction of Postfire Debris Flows Through Rainfall Anomaly
Maps.* Geophysical Research Letters, 52, e2025GL114791.
doi:10.1029/2025GL114791 — the **R\*** anomaly (R\*₁₅ₘᵢₙ = peak i15 ÷ 1-yr
Atlas-14 i15, R\*crit = 1) and the **M1 & R\*** model. Produced by `stormscape`,
not by this package.

Marc, O., Gosset, M., Saito, H., Uchida, T., & Malet, J.-P. (2019). *Spatial
Patterns of Storm-Induced Landslides and Their Relation to Rainfall Anomaly Maps.*
Geophysical Research Letters, 46, 11,167–11,177. doi:10.1029/2019GL083173 — origin
of the rainfall-anomaly concept, for *landslides*, normalizing event rainfall by
the **10-year** return period. Cavagnaro et al. (2025) adapt it to short-duration
i15 and the 1-year RI for runoff-generated debris flows.

Scheip, C., & Wegmann, K. (2022). *Insights on the growth and mobility of debris
flows from repeat high-resolution lidar.* Landslides, 19(6), 1297–1319 —
incremental along-path volume estimation from repeat topography.

Guido, L. E. (2025). *IncrementalDebrisFlowVolumeAnalyzer: A Python tool for
estimating volumes along a debris flow path.* Journal of Open Source Software —
**submitted 21 Aug 2025; draft, not yet published** (cite as in review). Builds on
Scheip & Wegmann (2022). A candidate downstream consumer of this package's change
raster.
