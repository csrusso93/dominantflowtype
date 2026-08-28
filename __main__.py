# -*- coding: utf-8 -*-
"""Launch the interactive workflow.

Works three ways, on PC / Mac / Linux, with no path edits:

    python -m dominantflowtype                     # from package_parent/
    exec(open(r".../dominantflowtype/__main__.py", encoding="utf-8").read())
    (run this file directly)

When launched as a module (``python -m``) the package is already importable and
we just call ``run()``.  When exec'd or run as a loose script the package is not
yet on ``sys.path``; we locate the folder that *contains* ``dominantflowtype``
and add it, so the same file bootstraps itself on any machine.
"""
from __future__ import annotations


def _bootstrap_run():
    import os
    import sys

    try:
        # Fast path: imported as part of the package (python -m dominantflowtype)
        from . import run
    except ImportError:
        # Exec'd / run directly: put the folder CONTAINING the package on sys.path.
        _here = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None
        _candidates = []
        if _here:
            _candidates.append(os.path.dirname(_here))  # parent of the package
        _candidates += [
            r"/path/to/parent/of/dominantflowtype",                     # PC
            os.path.expanduser("~/path/to/parent/of/dominantflowtype"),  # Mac
            os.path.expanduser("~/path/to/parent/of/dominantflowtype"),                    # Box Drive / Linux
        ]
        for _root in _candidates:
            if _root and os.path.isdir(os.path.join(_root, "dominantflowtype")):
                if _root not in sys.path:
                    sys.path.insert(0, _root)
                break
        else:
            raise FileNotFoundError(
                "Could not locate the folder containing 'dominantflowtype' on this machine"
            )
        from dominantflowtype import run

    return run()


if __name__ == "__main__":
    # Covers `python -m dominantflowtype` AND exec(open(__main__.py).read())
    # in the QGIS console (its global namespace has __name__ == "__main__").
    _bootstrap_run()
