"""Qt window for interactive core/boundary/surface selection."""

from pymol import cmd
from pymol.Qt import QtCore, QtWidgets

from . import sasa

Qt = QtCore.Qt

# Slider ranges per unit: (maximum, label)
UNIT_RANGE = {"abs": (100, "Å²"), "rel": (100, "%"), "iso": (100, "%")}


class CoreToSurfWindow(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(CoreToSurfWindow, self).__init__(parent, Qt.Window)
        self.setWindowTitle("Core / Boundary / Surface")
        self.unit = "abs"
        # Remember cutoffs per unit so switching units restores them.
        self.cutoffs = {u: [float(x) for x in v]
                        for u, v in sasa.DEFAULT_CUTOFFS.items()}
        self._build_ui()
        self.refresh_objects()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # Object picker
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Object:"))
        self.obj_combo = QtWidgets.QComboBox()
        self.obj_combo.setMinimumWidth(160)
        row.addWidget(self.obj_combo, 1)
        btn = QtWidgets.QPushButton("Refresh list")
        btn.clicked.connect(self.refresh_objects)
        row.addWidget(btn)
        layout.addLayout(row)

        # Units
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("SASA units:"))
        self.abs_radio = QtWidgets.QRadioButton("Absolute (Å²)")
        self.rel_radio = QtWidgets.QRadioButton("Relative (% of max ASA)")
        self.rel_radio.setToolTip(
            "SASA as % of the theoretical maximum ASA of the residue type "
            "(Tien et al. 2013).")
        self.iso_radio = QtWidgets.QRadioButton("Relative (% of isolated)")
        self.iso_radio.setToolTip(
            "SASA as % of the SASA of the same residue alone, in its own "
            "conformation. Works for any residue type.")
        self.abs_radio.setChecked(True)
        self.unit_radios = {"abs": self.abs_radio, "rel": self.rel_radio,
                            "iso": self.iso_radio}
        for radio in self.unit_radios.values():
            # Only react to the newly checked button, not the unchecked one.
            radio.toggled.connect(
                lambda checked: checked and self._on_unit_changed())
            row.addWidget(radio)
        row.addStretch(1)
        layout.addLayout(row)

        # Sliders, each with an edit box. The edit box accepts any
        # non-negative value, also beyond the slider range.
        grid = QtWidgets.QGridLayout()
        self.core_slider, self.core_edit, self.core_unit = self._make_slider()
        self.bound_slider, self.bound_edit, self.bound_unit =             self._make_slider()
        grid.addWidget(QtWidgets.QLabel("Core cutoff:"), 0, 0)
        grid.addWidget(self.core_slider, 0, 1)
        grid.addWidget(self.core_edit, 0, 2)
        grid.addWidget(self.core_unit, 0, 3)
        grid.addWidget(QtWidgets.QLabel("Boundary cutoff:"), 1, 0)
        grid.addWidget(self.bound_slider, 1, 1)
        grid.addWidget(self.bound_edit, 1, 2)
        grid.addWidget(self.bound_unit, 1, 3)
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)

        self.core_slider.valueChanged.connect(self._on_core_changed)
        self.bound_slider.valueChanged.connect(self._on_bound_changed)
        self.core_slider.sliderReleased.connect(self.update_selection)
        self.bound_slider.sliderReleased.connect(self.update_selection)
        self.core_edit.editingFinished.connect(self._on_core_edited)
        self.bound_edit.editingFinished.connect(self._on_bound_edited)

        # Live update checkbox
        self.live_check = QtWidgets.QCheckBox(
            "Live update (update while dragging)")
        self.live_check.setChecked(True)
        self.live_check.toggled.connect(self._on_live_toggled)
        layout.addWidget(self.live_check)

        # Buttons
        row = QtWidgets.QHBoxLayout()
        apply_btn = QtWidgets.QPushButton("Apply")
        apply_btn.clicked.connect(self.update_selection)
        recompute_btn = QtWidgets.QPushButton("Recompute SASA")
        recompute_btn.clicked.connect(self.recompute)
        reset_btn = QtWidgets.QPushButton("Reset cutoffs")
        reset_btn.setToolTip("Restore the default cutoffs: %s" %
                             self._defaults_text())
        reset_btn.clicked.connect(self.reset_cutoffs)
        row.addWidget(apply_btn)
        row.addWidget(reset_btn)
        row.addWidget(recompute_btn)
        layout.addLayout(row)

        # Export
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Format:"))
        self.format_combo = QtWidgets.QComboBox()
        self.format_combo.addItem("PyMOL selection", "pymol")
        self.format_combo.addItem("Rosetta ResidueIndexSelector", "rosetta")
        row.addWidget(self.format_combo, 1)
        copy_btn = QtWidgets.QPushButton("Print + copy to clipboard")
        copy_btn.clicked.connect(self.print_and_copy)
        row.addWidget(copy_btn)
        layout.addLayout(row)

        # Status line
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        self.obj_combo.currentIndexChanged.connect(
            lambda _i: self.update_selection())

        self._apply_unit_to_sliders()

    def _make_slider(self):
        slider = QtWidgets.QSlider(Qt.Horizontal)
        slider.setMinimumWidth(220)
        edit = QtWidgets.QLineEdit()
        edit.setFixedWidth(60)
        edit.setAlignment(Qt.AlignRight)
        unit = QtWidgets.QLabel()
        unit.setMinimumWidth(20)
        return slider, edit, unit

    # ------------------------------------------------------------- helpers
    def _apply_unit_to_sliders(self):
        maximum, _ = UNIT_RANGE[self.unit]
        for slider in (self.core_slider, self.bound_slider):
            slider.blockSignals(True)
            slider.setRange(0, maximum)
            slider.blockSignals(False)
        self._sync_widgets()

    @staticmethod
    def _defaults_text():
        return ", ".join(
            "%g / %g %s" % (core, bound, UNIT_RANGE[unit][1])
            for unit, (core, bound) in sasa.DEFAULT_CUTOFFS.items())

    def _sync_widgets(self):
        """Show the stored cutoffs in the sliders and edit boxes.

        Values outside the slider range put the slider at its end.
        """
        _, suffix = UNIT_RANGE[self.unit]
        core, bound = self.cutoffs[self.unit]
        for slider, edit, unit, value in (
                (self.core_slider, self.core_edit, self.core_unit, core),
                (self.bound_slider, self.bound_edit, self.bound_unit, bound)):
            slider.blockSignals(True)
            slider.setValue(int(round(value)))  # QSlider clamps to its range
            slider.blockSignals(False)
            edit.setText("%g" % value)
            unit.setText(suffix)

    def _maybe_update(self, slider):
        # With live update off, wait for the slider to be released. Changes
        # made with the keyboard or mouse wheel (slider not held down) still
        # apply immediately.
        if self.live_check.isChecked() or not slider.isSliderDown():
            self.update_selection()

    def _set_cutoff(self, which, value):
        """Set the core (0) or boundary (1) cutoff, keeping core <= boundary."""
        cutoffs = self.cutoffs[self.unit]
        cutoffs[which] = float(value)
        if cutoffs[0] > cutoffs[1]:
            # Move the other cutoff along with the one being changed.
            cutoffs[1 - which] = cutoffs[which]
        self._sync_widgets()

    @staticmethod
    def _parse(edit):
        """Return the edit box value, or None if it is not a valid number."""
        try:
            value = float(edit.text().strip().replace(",", "."))
        except ValueError:
            return None
        if value != value or value < 0:  # NaN or negative
            return None
        return value

    # -------------------------------------------------------------- slots
    def _on_core_changed(self, value):
        self._set_cutoff(0, value)
        self._maybe_update(self.core_slider)

    def _on_bound_changed(self, value):
        self._set_cutoff(1, value)
        self._maybe_update(self.bound_slider)

    def _on_edited(self, which, edit):
        value = self._parse(edit)
        if value is None:
            self.status.setText("Invalid cutoff '%s'." % edit.text())
            self._sync_widgets()  # restore the previous value
            return
        if value == self.cutoffs[self.unit][which]:
            self._sync_widgets()  # normalize the text, nothing changed
            return
        self._set_cutoff(which, value)
        self.update_selection()

    def _on_core_edited(self):
        self._on_edited(0, self.core_edit)

    def _on_bound_edited(self):
        self._on_edited(1, self.bound_edit)

    def _on_unit_changed(self):
        self.unit = next(u for u, radio in self.unit_radios.items()
                         if radio.isChecked())
        self._apply_unit_to_sliders()
        self.update_selection()

    def _on_live_toggled(self, _checked):
        self.update_selection()

    # ------------------------------------------------------------ actions
    def refresh_objects(self):
        current = self.obj_combo.currentText()
        names = [n for n in cmd.get_names("objects")
                 if cmd.get_type(n) == "object:molecule"]
        self.obj_combo.blockSignals(True)
        self.obj_combo.clear()
        self.obj_combo.addItems(names)
        if current in names:
            self.obj_combo.setCurrentIndex(names.index(current))
        self.obj_combo.blockSignals(False)
        self.update_selection()

    def reset_cutoffs(self):
        # Restore the defaults for both units, then show the current one.
        self.cutoffs = {u: [float(x) for x in v]
                        for u, v in sasa.DEFAULT_CUTOFFS.items()}
        self._apply_unit_to_sliders()
        self.update_selection()

    def recompute(self):
        obj = self.obj_combo.currentText()
        if obj:
            sasa.clear_cache(obj)
        self.update_selection()

    def update_selection(self):
        obj = self.obj_combo.currentText()
        if not obj:
            self.status.setText("No molecular object loaded. "
                                "Load one and press 'Refresh list'.")
            return
        if obj not in cmd.get_names("objects"):
            self.status.setText("Object '%s' no longer exists." % obj)
            return
        core, bound = self.cutoffs[self.unit]
        try:
            counts = sasa.apply_layers(obj, core, bound, self.unit)
        except Exception as e:  # report errors in the window, not a traceback
            self.status.setText("Error: %s" % e)
            return
        self.status.setText(
            "core: %(core)d   boundary: %(boundary)d   surface: %(surface)d "
            "residues" % counts)

    def print_and_copy(self):
        # Make sure the exported residues match the current slider values.
        self.update_selection()
        fmt = self.format_combo.currentData()
        try:
            text, warnings = sasa.format_layers(fmt=fmt)
        except Exception as e:
            self.status.setText("Error: %s" % e)
            return
        print(text)
        for w in warnings:
            print(" CoreToSurf warning: " + w)
        QtWidgets.QApplication.clipboard().setText(text)
        msg = "Selections printed and copied to the clipboard (%s)." % \
            self.format_combo.currentText()
        if warnings:
            msg += " Some residues were omitted, see the PyMOL log."
        self.status.setText(msg)


_window = None


def show_window():
    global _window
    if _window is None:
        # No parent: an owned window gets no taskbar button on Windows.
        _window = CoreToSurfWindow()
        main = None
        try:
            from pymol import gui
            main = gui.get_qtwindow()
        except Exception:
            pass
        if main is not None:
            _window.setWindowIcon(main.windowIcon())
            # Without a parent the window would outlive PyMOL's main window.
            main.destroyed.connect(_window.close)
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(_window.close)
    else:
        _window.refresh_objects()
    _window.show()
    _window.raise_()
    _window.activateWindow()
    return _window
