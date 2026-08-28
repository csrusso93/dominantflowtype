# -*- coding: utf-8 -*-
"""Interactive, one-field-at-a-time prompts.

Three back-ends, chosen automatically:

* **scripted** — if ``answers`` is supplied, replay them in order (testing /
  headless runs). Takes precedence over everything.
* **Qt dialogs** — when a QApplication exists (i.e. running in the QGIS Python
  Console / QGIS GUI). ``input()`` does NOT work there ("lost sys.stdin"), so we
  use QFileDialog / QInputDialog / QMessageBox instead.
* **terminal** — plain ``input()`` when there is no GUI (a real shell).
"""
from __future__ import annotations

import os

RASTER_FILTER = "Raster DEM (*.tif *.tiff *.vrt *.img *.dem);;All files (*.*)"
VECTOR_FILTER = "Vector (*.shp *.gpkg);;All files (*.*)"
CLOUD_FILTER = "Point cloud (*.las *.laz *.copc *.copc.laz *.ept.json);;All files (*.*)"


class PromptCancelled(Exception):
    """User cancelled/aborted an interactive prompt.

    A plain ``Exception`` (NOT ``SystemExit``): raising ``SystemExit`` inside the
    QGIS Python Console propagates into the embedded interpreter and shuts down
    the whole QGIS application. This is caught by ``run()`` and returns cleanly.
    """


def _sb(qmessagebox, name):
    """Return QMessageBox.StandardButton.<name> across PyQt5 (QGIS 3) / PyQt6 (QGIS 4)."""
    holder = getattr(qmessagebox, "StandardButton", qmessagebox)
    return getattr(holder, name)


class Prompter:
    def __init__(self, answers=None, use_gui=None):
        self._scripted = list(answers) if answers else None
        self._gui = self._detect_gui() if use_gui is None else use_gui

    # -- back-end detection ---------------------------------------------------
    @staticmethod
    def _detect_gui():
        try:
            from qgis.PyQt.QtWidgets import QApplication
            return QApplication.instance() is not None
        except Exception:
            return False

    @staticmethod
    def _parent():
        try:
            from qgis.utils import iface
            if iface is not None:
                return iface.mainWindow()
        except Exception:
            pass
        return None

    def _next_scripted(self, label):
        val = self._scripted.pop(0)
        print(label + "  >> " + repr(val))
        return val

    # -- yes / no -------------------------------------------------------------
    def ask_yes_no(self, prompt, default=True):
        if self._scripted is not None:
            v = str(self._next_scripted(prompt)).strip().lower()
            return default if v == "" else v in ("y", "yes", "true", "1")
        if self._gui:
            from qgis.PyQt.QtWidgets import QMessageBox
            yes, no = _sb(QMessageBox, "Yes"), _sb(QMessageBox, "No")
            res = QMessageBox.question(self._parent(), "dominantflowtype",
                                       prompt, yes | no, yes if default else no)
            return res == yes
        d = "Y/n" if default else "y/N"
        while True:
            v = input(f"{prompt} [{d}]: ").strip().lower()
            if v == "":
                return default
            if v in ("y", "yes"):
                return True
            if v in ("n", "no"):
                return False
            print("  (please answer y or n)")

    # -- free text ------------------------------------------------------------
    def ask_text(self, prompt, default=None, allow_blank=True):
        if self._scripted is not None:
            v = str(self._next_scripted(prompt)).strip()
            if v == "" and default is not None:
                return default
            return v
        if self._gui:
            from qgis.PyQt.QtWidgets import QInputDialog
            text, ok = QInputDialog.getText(self._parent(), "dominantflowtype",
                                            prompt, text=str(default or ""))
            if not ok:
                if allow_blank:
                    return default if default is not None else ""
                raise PromptCancelled(f"cancelled at: {prompt}")
            text = text.strip()
            if text == "" and default is not None:
                return default
            return text
        suffix = f" [{default}]" if default not in (None, "") else ""
        while True:
            val = input(f"{prompt}{suffix}: ").strip()
            if val == "" and default is not None:
                return default
            if val == "" and allow_blank:
                return ""
            if val != "":
                return val
            print("  (a value is required)")

    # -- file path ------------------------------------------------------------
    def ask_path(self, prompt, must_exist=True, skippable=False,
                 skip_word="skip", file_filter="All files (*.*)"):
        if self._scripted is not None:
            v = str(self._next_scripted(prompt)).strip().strip('"').strip("'")
            if skippable and v.lower() in (skip_word, "", "none", "n/a"):
                return None
            return v
        if self._gui:
            from qgis.PyQt.QtWidgets import QFileDialog, QMessageBox
            if skippable:
                yes, no = _sb(QMessageBox, "Yes"), _sb(QMessageBox, "No")
                res = QMessageBox.question(
                    self._parent(), "dominantflowtype",
                    f"{prompt}\n\nProvide this file?   (No = skip)",
                    yes | no, yes)
                if res != yes:
                    return None
            path, _ = QFileDialog.getOpenFileName(
                self._parent(), prompt, "", file_filter)
            if not path:
                if skippable:
                    return None
                raise PromptCancelled(f"no file selected for: {prompt}")
            return path
        while True:
            extra = f" (or type '{skip_word}')" if skippable else ""
            val = input(f"{prompt}{extra}: ").strip().strip('"').strip("'")
            if skippable and val.lower() in (skip_word, "", "none", "n/a"):
                return None
            if not must_exist or os.path.exists(val):
                return val
            print(f"  (path not found: {val})")

    # -- directory ------------------------------------------------------------
    def ask_dir(self, prompt, default=""):
        if self._scripted is not None:
            v = str(self._next_scripted(prompt)).strip().strip('"').strip("'")
            return v if v else default
        if self._gui:
            from qgis.PyQt.QtWidgets import QFileDialog
            d = QFileDialog.getExistingDirectory(self._parent(), prompt,
                                                 default or "")
            return d or default
        val = input(f"{prompt} [{default}]: ").strip().strip('"').strip("'")
        return val or default

    # -- number ---------------------------------------------------------------
    def ask_float(self, prompt, default=None):
        if self._scripted is not None:
            v = str(self._next_scripted(prompt)).strip()
            return float(v) if v != "" else float(default)
        if self._gui:
            from qgis.PyQt.QtWidgets import QInputDialog
            val, ok = QInputDialog.getDouble(
                self._parent(), "dominantflowtype", prompt,
                float(default) if default is not None else 0.0)
            return float(val) if ok else float(default)
        while True:
            val = input(
                f"{prompt}" + (f" [{default}]" if default is not None else "") + ": "
            ).strip()
            if val == "" and default is not None:
                return float(default)
            try:
                return float(val)
            except ValueError:
                print("  (enter a number)")
