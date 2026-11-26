# CoreToSurf

A PyMOL 3 plugin for interactive classification of residues into **core**
(tv_orange), **boundary** (tv_green) and **surface** (marine) based on solvent
accessible surface area (SASA).

![CoreToSurf: moving the cutoff sliders updates the core, boundary and surface selections](assets/coretosurf_readme_animation.gif)

## Installation

### Quick install from the PyMOL prompt

Paste this line into the PyMOL command line. It installs the package into
PyMOL's own Python and adds *Plugin → Core Boundary Surface* to the menu:

```
import sys, subprocess, importlib; subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "https://github.com/ajasja/CoreToSurf/archive/refs/heads/main.zip"]); importlib.invalidate_caches(); import coretosurf; coretosurf.install_plugin()
```

The menu item appears immediately and on every later start. Type `coretosurf`
(or use the menu) to open the window. Run the same line again to upgrade. Git
does not need to be installed.

### Install from a shell

Install with the **same Python environment that PyMOL uses**. PyMOL itself is
not a pip dependency (it is installed with conda or the Schrödinger installer).

```
# from a clone of this repository
python -m pip install .
# or straight from git
python -m pip install git+https://github.com/ajasja/CoreToSurf.git
```

For a conda install of PyMOL, activate its environment first, or call its
Python directly, e.g. `C:\path\to\conda\envs\pymol\python -m pip install .`.
For development use `pip install -e .`.

### Using it

**Option A: load on demand.** At the PyMOL prompt (or in `~/.pymolrc`):

```
import coretosurf
coretosurf          # opens the window
```

**Option B: add it to the Plugin menu (one time).** Register the plugin once:

```
coretosurf install-plugin           # from a shell
python -m coretosurf install-plugin # same, if the Scripts dir is not on PATH
```

or, at the PyMOL prompt, `import coretosurf; coretosurf.install_plugin()`.
This writes a small loader (`coretosurf.py`) into PyMOL's user plugin
directory (`%APPDATA%\pymol\startup` on Windows, `~/.pymol/startup`
elsewhere, or the first writable directory already on PyMOL's plugin path)
and adds that directory to PyMOL's plugin path. *Plugin → Core Boundary
Surface* is then available on every start. Because the loader only imports the
installed package, upgrading with pip does not require registering again.
Remove it with `coretosurf uninstall-plugin`.

**Option C: classic plugin install without pip.** *Plugin → Plugin Manager →
Install New Plugin* and choose a zip of the `src/coretosurf` folder.

### Uninstalling

Paste this line into the PyMOL command line. It removes the Plugin menu entry
and then uninstalls the package from PyMOL's Python:

```
import sys, subprocess, coretosurf; coretosurf.uninstall_plugin(); subprocess.check_call([sys.executable, "-m", "pip", "uninstall", "-y", "coretosurf"])
```

Restart PyMOL afterwards. From a shell, the same is
`coretosurf uninstall-plugin` followed by `python -m pip uninstall coretosurf`.
Remove the plugin entry first: if only the package is removed, the leftover
loader in PyMOL's startup directory cannot find it when PyMOL starts.
If you installed with Option C, remove it in *Plugin → Plugin Manager*.

## Window

- **Object**: the molecular object to classify (*Refresh list* after loading new ones).
- **SASA units** (each unit remembers its own cutoffs):
  - *Absolute* (Å², slider 0–100).
  - *Relative (% of max ASA)*: % of the theoretical maximum ASA of the residue
    type from Tien et al. 2013 (slider 0–100). Residue types without a
    reference value use 200 Å².
  - *Relative (% of isolated)*: % of the SASA of the same residue on its own,
    in its own conformation, with only its own atoms occluding (slider 0–100).
    It needs no reference table, so it works for any polymer residue
    (nucleotides, modified residues) and accounts for the actual side-chain
    conformation.
- **Core cutoff**: residues with SASA below this value are core.
- **Boundary cutoff**: residues with SASA from the core cutoff up to this value are
  boundary. Residues at or above it are surface. The two cutoffs are kept so that core ≤ boundary.
- Next to each slider is an edit box. Type a value and press Enter (or leave the box)
  to apply it. Typed values can have decimals and can be larger than the slider
  maximum, in which case the slider stays at its end.
- **Live update**: when checked, the selections update while a slider is dragged;
  when unchecked, they update when the slider is released.
- **Reset cutoffs**: restore the default cutoffs (for both units).
- **Recompute SASA**: discard the cached SASA and recalculate.
- **Format** + **Print + copy to clipboard**: print the three layers to the
  PyMOL log and copy them to the clipboard, using PDB residue numbers and chains:
  - *PyMOL selection* (default), ready to paste back into PyMOL:
    `select core, (chain A and resi 1-5+9) or (chain B and resi 3+5)`
  - *Rosetta ResidueIndexSelector* (`resnums` format):
    `core: 1A-5A,9A,3B,5B`.
    This format cannot express residues with insertion codes, negative residue numbers or
    blank or multi-letter chain ids. Such residues are omitted, and a warning
    in the PyMOL log lists them.

The plugin creates three selections named `core`, `boundary` and `surface`, and
colors their carbon atoms. The representation is not changed.

Default cutoffs: 20 / 40 Å² (Rosetta LayerSelector SASA mode) and 20 / 50 %
for both relative units.

## PyMOL command

```
coretosurf                      # open the GUI
coretosurf 1ubq, 20, 40         # absolute Å²
coretosurf 1ubq, 20, 50, rel    # % of max ASA (Tien et al. 2013)
coretosurf 1ubq, 20, 50, iso    # % of the isolated residue
coretosurf 1ubq, recompute=1    # force SASA recalculation
```

### `print_sel`: residues of any selection

`print_sel selection [, format [, copy]]` prints the residues of any PyMOL
selection with PDB numbering and chain, and copies the result to the clipboard.

```
print_sel core                                  # (chain A and resi 3+5+13-17) or ...
print_sel visible and (core or surface)         # any selection expression
print_sel chain A and boundary, rosetta         # 7A,21A,36A,44A-45A
print_sel core, copy=0                          # print only
```

The format is `pymol` (default) or `rosetta`, as described above. Residues
from several objects that share a chain and residue number appear once.

## How SASA is computed

- Whole-residue SASA is calculated once per object and state with `get_area`
  (`dot_solvent=1`, `dot_density=3`) on a temporary copy of the object.
- Waters and hydrogens are removed before the calculation. Ligands and other
  hetero atoms are kept so they occlude the polymer, but only polymer residues
  are classified.
- For the *% of isolated* unit, the isolated SASA of every polymer residue is
  computed in the same pass: ligands are removed from the temporary copy, the
  residues are moved rigidly apart on a grid so that none occludes another,
  and `get_area` is run again with the same settings.
- The result is cached. Moving a slider only reclassifies the cached values.
- The cache is invalidated automatically when the object's atom count changes.
  After other edits, such as moving atoms, changing coordinates or renumbering
  residues, press *Recompute SASA*.
