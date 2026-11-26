"""CoreToSurf: interactive core / boundary / surface residue selection.

Ways to use it:

* ``import coretosurf`` in PyMOL (or in ``.pymolrc``) adds the ``coretosurf``
  command; ``coretosurf`` without arguments opens the window.
* Optional one-time registration in the Plugin menu: run
  ``coretosurf install-plugin`` in a shell, or ``coretosurf.install_plugin()``
  at the PyMOL prompt.
* Classic install: Plugin > Plugin Manager > Install New Plugin with this
  folder (or a zip of it).
"""

__version__ = "0.1.0"

try:
    from pymol import cmd
except ImportError:  # e.g. `coretosurf --version` without PyMOL
    cmd = None


def coretosurf(selection="", core=None, boundary=None, unit="abs",
               state=-1, recompute=0, quiet=0):
    """
DESCRIPTION

    Classify polymer residues into core, boundary and surface by SASA and
    create the selections "core" (tv_orange), "boundary" (tv_green) and
    "surface" (marine). SASA is computed once per object and cached.

USAGE

    coretosurf [object [, core [, boundary [, unit [, state [, recompute ]]]]]]

ARGUMENTS

    object = molecular object (default: first object)
    core = residues with SASA below this are core
    boundary = residues with SASA below this (and >= core) are boundary
    unit = abs (Å², default 20/40), rel (% of max ASA, default 20/50) or
           iso (% of the SASA of the residue in isolation, default 20/50)
    recompute = 1 to discard the cached SASA first

    With no arguments, opens the GUI.
    """
    from . import sasa
    if not selection:
        if core is None and boundary is None:
            open_gui()
            return
        names = [n for n in cmd.get_names("objects")
                 if cmd.get_type(n) == "object:molecule"]
        if not names:
            raise ValueError("no molecular object loaded")
        selection = names[0]
    unit = str(unit).lower()
    if unit not in sasa.DEFAULT_CUTOFFS:
        raise ValueError("unit must be 'abs', 'rel' or 'iso'")
    dcore, dbound = sasa.DEFAULT_CUTOFFS[unit]
    core = dcore if core is None else float(core)
    boundary = dbound if boundary is None else float(boundary)
    if int(recompute):
        sasa.clear_cache(selection)
    counts = sasa.apply_layers(selection, core, boundary, unit, int(state))
    if not int(quiet):
        print(" CoreToSurf: core %(core)d, boundary %(boundary)d, "
              "surface %(surface)d residues" % counts)
    return counts


def copy_to_clipboard(text):
    """Copy text to the clipboard. Returns False if no clipboard is available."""
    try:
        from pymol.Qt import QtWidgets
        if QtWidgets.QApplication.instance() is None:
            return False
        QtWidgets.QApplication.clipboard().setText(text)
        return True
    except Exception:
        return False


def print_sel(selection="sele", format="pymol", copy=1):
    """
DESCRIPTION

    Print the residues of any selection with PDB numbering and chain,
    and copy the result to the clipboard.

USAGE

    print_sel [ selection [, format [, copy ]]]

ARGUMENTS

    selection = any PyMOL selection {default: sele}
    format = pymol (default): (chain A and resi 1-5+9) or (chain B and resi 3)
             rosetta: 1A-5A,9A,3B  (ResidueIndexSelector resnums; residues
             with insertion codes, negative numbers or without a
             single-letter chain id are omitted with a warning)
    copy = 0 to only print {default: 1}

EXAMPLES

    print_sel core
    print_sel visible and (core or surface), rosetta
    print_sel chain A and boundary, copy=0
    """
    from . import sasa
    text, skipped = sasa.format_selection(selection, format)
    print(text)
    if skipped:
        print(" print_sel warning: %d residue(s) with insertion codes, "
              "negative numbers or without a single-letter chain id "
              "omitted: %s" % (len(skipped), ", ".join(skipped)))
    if int(copy) and not copy_to_clipboard(text):
        print(" print_sel: clipboard not available")
    return text


def open_gui():
    from . import gui
    return gui.show_window()


def install_plugin(plugdir=None, quiet=False):
    """Register CoreToSurf in PyMOL's Plugin menu (one-time, optional).

    When called inside a running PyMOL GUI, the menu item is also added
    right away.
    """
    from . import registration
    stub = registration.install_plugin(plugdir, quiet=True)
    added_now = False
    try:
        __init_plugin__()
        added_now = True
    except Exception:
        pass
    if not quiet:
        print(" CoreToSurf: plugin loader written to %s" % stub)
        if added_now:
            print(" CoreToSurf: 'Core Boundary Surface' added to the "
                  "Plugin menu.")
    return stub


def uninstall_plugin(quiet=False):
    """Remove the Plugin menu registration (takes effect after restart)."""
    from . import registration
    return registration.uninstall_plugin(quiet=quiet)


_menu_added = False


def __init_plugin__(app=None):
    global _menu_added
    if _menu_added:
        return
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Core Boundary Surface", open_gui)
    _menu_added = True


if cmd is not None:
    cmd.extend("coretosurf", coretosurf)
    cmd.auto_arg[0]["coretosurf"] = cmd.auto_arg[0]["zoom"]
    cmd.auto_arg[3]["coretosurf"] = [
        lambda: cmd.Shortcut(["abs", "rel", "iso"]), "unit", ", "]
    cmd.extend("print_sel", print_sel)
    cmd.auto_arg[0]["print_sel"] = cmd.auto_arg[0]["zoom"]
    cmd.auto_arg[1]["print_sel"] = [
        lambda: cmd.Shortcut(["pymol", "rosetta"]), "format", ", "]
