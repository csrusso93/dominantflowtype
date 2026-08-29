# -*- coding: utf-8 -*-
"""
dominantflowtype
================
Post-event, DEM-based estimation of the dominant flow type (runoff-generated
*debris flow* vs. *flood*) in steep, burned channels, following:

    Cavagnaro, D. B., McCoy, S. W., Kean, J. W., Thomas, M. A., Lindsay, D. N.,
    McArdell, B. W., & Hirschberg, J. (2024). A Robust Quantitative Method to
    Distinguish Runoff-Generated Debris Flows From Floods. Geophysical Research
    Letters, 51, e2024GL109768. https://doi.org/10.1029/2024GL109768

Core diagnostic (Cavagnaro et al., 2024):

    Q*  =  Q_peak / Q_fluv  =  (v * A_xs) / (I * A_us)

        Q_peak = v * A_xs                      observed peak discharge
        Q_fluv = A_us * I                      theoretical max clear-water discharge
        v      = Fr * sqrt(g * R),  Fr = 1     Froude-critical velocity
        R      = A_xs / P                       hydraulic radius (paper: "h")
        A_xs   = inundated cross-sectional area (from trimlines + DEM bed)
        P      = wetted perimeter
        A_us   = upstream drainage area (DEM flow routing, clipped to watershed)
        I      = 30-min rainfall intensity, I30 (MRMS radar / Synoptic gauges /
                 9.7 mm/hr constant)

    Q* > 1  -> debris flow   |   Q* <= 1 -> flood   (~92% accuracy; threshold=1)

Run interactively inside the QGIS Python console.  Add the folder that
*contains* this package to ``sys.path`` first -- set ``ROOT`` to wherever you
cloned it::

    >>> import sys
    >>> ROOT = r"/path/to/parent/of/dominantflowtype"
    >>> sys.path.insert(0, ROOT)
    >>> import dominantflowtype as dft
    >>> dft.run()

or, to run the single script directly (any OS)::

    >>> import os
    >>> exec(open(os.path.join(ROOT, "dominantflowtype", "__main__.py"),
    ...           encoding="utf-8").read())

``__main__.py`` resolves the package location on its own -- from its own
``__file__``, from ``$DOMINANTFLOWTYPE_ROOT``, or by walking up from the working
directory -- so running it needs no path edits.

You will be prompted, one field at a time, for an (optional) OpenTopography
key, pre- and post-event DEMs (either skippable), trimlines, thalweg and
watershed layers, and Synoptic rainfall parameters.  If ONE DEM is skipped the
volume (DEM-of-difference) products are skipped and only Q* is computed.

Package layout
--------------
    _compat.py       shared GDAL / QGIS imports + shims
    config.py        Config (all tunable parameters)
    io_prompts.py    Prompter (interactive console input)
    raster.py        RasterSampler, build_dod()
    geometry.py      vector/CRS helpers, cast_transects(), trimline_hits()
    cross_section.py measure_cross_section()  (A_xs, P, R, depth, width)
    hydrology.py     FlowRouter  (D8 A_us, priority-flood fill)
    aoi.py           load_aoi()  (bbox/vector/QGIS layer -> WGS84 bounds)
    rainfall.py      rainfall_i30() dispatcher; synoptic_i30(), synoptic_i30_bbox()
    mrms.py          mrms_i30(), i_storm_day()  (NOAA MRMS radar; no account)
    dem_3dep.py      get_dem()  (USGS 3DEP; py3dep or GDAL ImageServer)
    opentopo.py      opentopo_search(), opentopo_download()
    metrics.py       velocity(), compute_qstar(), usability()
    outputs.py       write_workbook(), add_qgis_layers(), styling
    pipeline.py      Pipeline, run()  (orchestration incl. reach volumes)
"""
from __future__ import annotations

__version__ = "0.1.0"

from ._compat import _HAS_QGIS
from .config import Config
from .io_prompts import Prompter
from .pipeline import Pipeline, run

# selected internals re-exported for scripting / testing
from .raster import RasterSampler, build_dod
from .cross_section import measure_cross_section
from .metrics import velocity, compute_qstar, usability
from .rainfall import (
    synoptic_i30, synoptic_i30_bbox, synoptic_token, rainfall_i30,
)
from .hydrology import FlowRouter
from .geometry import cast_transects, trimline_hits
from . import aoi, mrms, dem_3dep, pointcloud, bridge, acquire
from .mrms import mrms_i30
from .dem_3dep import get_dem as get_dem_3dep
from .pointcloud import diagnose as pointcloud_diagnose, prepare_from_clouds
from .bridge import export_bridge
from .acquire import list_3dep_projects, acquire_3dep_cloud
from .smoketest import smoke_config, make_test_aoi

__all__ = [
    "Config", "Prompter", "Pipeline", "run", "_HAS_QGIS",
    "RasterSampler", "build_dod", "measure_cross_section",
    "velocity", "compute_qstar", "usability",
    "synoptic_i30", "synoptic_i30_bbox", "synoptic_token", "rainfall_i30",
    "FlowRouter", "cast_transects", "trimline_hits",
    # stormscape-derived import methods
    "aoi", "mrms", "dem_3dep", "mrms_i30", "get_dem_3dep",
    # Stage 0 point clouds + IDFVA bridge
    "pointcloud", "bridge", "acquire", "pointcloud_diagnose",
    "prepare_from_clouds", "export_bridge",
    "list_3dep_projects", "acquire_3dep_cloud",
    "smoke_config", "make_test_aoi",
]
