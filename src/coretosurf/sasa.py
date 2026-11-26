"""SASA calculation, caching and core/boundary/surface classification.

This module has no Qt dependency, so it can be used from scripts and from the
``coretosurf`` PyMOL command as well as from the GUI.
"""

import math

from pymol import cmd

# Theoretical maximum ASA per residue (Å²), Tien et al. 2013, PLoS ONE 8:e80635.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
# Common non-standard / force-field residue names mapped to their parent.
ALIASES = {
    "MSE": "MET", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS",
    "HSE": "HIS", "HSP": "HIS", "CYX": "CYS", "CYM": "CYS", "ASH": "ASP",
    "GLH": "GLU", "LYN": "LYS", "SEP": "SER", "TPO": "THR", "PTR": "TYR",
}
DEFAULT_MAX_ASA = 200.0  # fallback for residues without a reference value

LAYERS = ("core", "boundary", "surface")
COLORS = {"core": "tv_orange", "boundary": "tv_green", "surface": "marine"}

# SASA units:
#   abs - absolute SASA in Å²
#   rel - % of the theoretical maximum ASA of the residue type (MAX_ASA)
#   iso - % of the SASA of the same residue in isolation: the residue in its
#         own conformation with only its own atoms occluding (as sasa_gemmi's
#         sasa_pct). Needs no reference table, so it works for any residue
#         type (nucleotides, modified residues) and accounts for the actual
#         side-chain conformation.
UNITS = ("abs", "rel", "iso")

# Default cutoffs per unit: (core, boundary).
# Absolute values follow Rosetta's LayerSelector SASA mode.
DEFAULT_CUTOFFS = {"abs": (20.0, 40.0), "rel": (20.0, 50.0),
                   "iso": (20.0, 50.0)}

# Upper bound of an atom's probe-extended radius (vdW + 1.4 Å probe), used to
# keep isolated residues far enough apart not to see each other.
_MAX_EXT_RADIUS = 5.0

# Cache: (object, state) -> {"natoms": int, "residues": [Residue, ...]}
_cache = {}


class Residue(object):
    __slots__ = ("segi", "chain", "resi", "resn", "sasa", "sasa_iso",
                 "rel", "iso")

    def __init__(self, segi, chain, resi, resn, sasa, sasa_iso):
        self.segi, self.chain, self.resi, self.resn = segi, chain, resi, resn
        self.sasa = sasa
        self.sasa_iso = sasa_iso
        ref = MAX_ASA.get(ALIASES.get(resn, resn), DEFAULT_MAX_ASA)
        self.rel = 100.0 * sasa / ref
        self.iso = 100.0 * sasa / sasa_iso if sasa_iso > 0.0 else 0.0

    def value(self, unit):
        if unit == "rel":
            return self.rel
        if unit == "iso":
            return self.iso
        return self.sasa


def _spread_residues(obj):
    """Move the residues of ``obj`` (state 1) apart on a cubic grid.

    Afterwards no atom is within reach of an atom of another residue, so one
    ``get_area`` call gives the isolated SASA of every residue. Each residue
    is translated rigidly, so its conformation is unchanged.
    """
    coords, order = {}, []

    def add(segi, chain, resi, x, y, z):
        key = (segi, chain, resi)
        if key not in coords:
            coords[key] = []
            order.append(key)
        coords[key].append((x, y, z))

    cmd.iterate_state(1, obj, "add(segi, chain, resi, x, y, z)",
                      space={"add": add})
    if not order:
        return

    centers, rmax = {}, 0.0
    for key, pts in coords.items():
        n = float(len(pts))
        c = (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n,
             sum(p[2] for p in pts) / n)
        centers[key] = c
        for p in pts:
            rmax = max(rmax, math.sqrt((p[0] - c[0]) ** 2 +
                                       (p[1] - c[1]) ** 2 +
                                       (p[2] - c[2]) ** 2))
    spacing = 2.0 * (rmax + _MAX_EXT_RADIUS) + 1.0
    side = max(1, int(math.ceil(len(order) ** (1.0 / 3.0))))
    while side ** 3 < len(order):  # guard against rounding of the cube root
        side += 1

    shift = {}
    for i, key in enumerate(order):
        c = centers[key]
        shift[key] = ((i % side) * spacing - c[0],
                      (i // side % side) * spacing - c[1],
                      (i // (side * side)) * spacing - c[2])

    def move(segi, chain, resi, x, y, z):
        d = shift[(segi, chain, resi)]
        return x + d[0], y + d[1], z + d[2]

    cmd.alter_state(1, obj, "(x, y, z) = move(segi, chain, resi, x, y, z)",
                    space={"move": move})


def _residue_sums(selection):
    """{(segi, chain, resi): [resn, sum of b]} and the key order."""
    sums, order = {}, []

    def add(segi, chain, resi, resn, b):
        key = (segi, chain, resi)
        if key not in sums:
            sums[key] = [resn, 0.0]
            order.append(key)
        sums[key][1] += b

    cmd.iterate(selection, "add(segi, chain, resi, resn, b)",
                space={"add": add})
    return sums, order


def _compute(obj, state):
    """Per-residue SASA of polymer residues in ``obj``.

    Waters and hydrogens are removed first. Ligands and other hetero atoms are
    kept so that they occlude the polymer, but only polymer residues are
    reported. For each residue the isolated SASA (the residue alone, in the
    same conformation) is computed as well, as the reference for the ``iso``
    unit.
    """
    tmp = cmd.get_unused_name("_c2s_tmp")
    try:
        cmd.create(tmp, "(%s) and not solvent and not hydro" % obj,
                   source_state=state, target_state=1, zoom=0)
        cmd.set("dot_solvent", 1, tmp)
        cmd.set("dot_density", 3, tmp)
        cmd.get_area(tmp, state=1, load_b=1)
        sums, order = _residue_sums("%s and polymer" % tmp)

        # Isolated reference: keep only the polymer, pull the residues apart
        # and compute the SASA again with the same settings.
        cmd.remove("%s and not polymer" % tmp)
        _spread_residues(tmp)
        cmd.get_area(tmp, state=1, load_b=1)
        iso, _ = _residue_sums(tmp)
    finally:
        cmd.delete(tmp)

    return [Residue(k[0], k[1], k[2], sums[k][0], sums[k][1],
                    iso.get(k, ("", 0.0))[1])
            for k in order]


def get_residues(obj, state=-1, force=False):
    """Return cached per-residue SASA for ``obj``, computing it if needed.

    The cache is invalidated automatically if the atom count of the object
    changes (e.g. atoms were added or removed).
    """
    if state in (-1, 0):
        state = cmd.get_state()
    key = (obj, state)
    natoms = cmd.count_atoms(obj)
    entry = _cache.get(key)
    if force or entry is None or entry["natoms"] != natoms:
        entry = {"natoms": natoms, "residues": _compute(obj, state)}
        _cache[key] = entry
    return entry["residues"]


def clear_cache(obj=None):
    if obj is None:
        _cache.clear()
    else:
        for key in [k for k in _cache if k[0] == obj]:
            del _cache[key]


def classify(residues, core_cut, boundary_cut, unit="abs"):
    """Split residues into layers. Returns {layer: [Residue, ...]}."""
    out = {name: [] for name in LAYERS}
    for r in residues:
        v = r.value(unit)
        if v < core_cut:
            out["core"].append(r)
        elif v < boundary_cut:
            out["boundary"].append(r)
        else:
            out["surface"].append(r)
    return out


def _escape_resi(resi):
    # A leading minus would be read as a range operator.
    return "\\" + resi if resi.startswith("-") else resi


def _selection_expr(obj, residues):
    """Build one compact selection expression for a list of residues."""
    if not residues:
        return "none"
    groups = {}
    for r in residues:
        groups.setdefault((r.segi, r.chain), []).append(_escape_resi(r.resi))
    parts = []
    for (segi, chain), resis in groups.items():
        parts.append("(segi '%s' and chain '%s' and resi %s)"
                     % (segi, chain, "+".join(resis)))
    return "(%s) and (%s)" % (obj, " or ".join(parts))


FORMATS = ("pymol", "rosetta")

# Result of the most recent apply_layers call: {layer: [Residue, ...]}
last_layers = None


def _split_resi(resi):
    """'52' -> (52, ''), '52A' -> (52, 'A'), '-3' -> (-3, ''). None if odd."""
    i = len(resi)
    while i > 0 and not resi[i - 1].isdigit():
        i -= 1
    try:
        return int(resi[:i]), resi[i:]
    except ValueError:
        return None


def _runs(numbers):
    """Sorted unique ints -> list of (start, end) runs of consecutive ints."""
    runs = []
    for n in numbers:
        if runs and n == runs[-1][1] + 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    return runs


def _group_by_chain(residues):
    """{chain: (sorted plain resi ints, [resi with insertion code / odd])}.

    Chains keep their order of first appearance. Segment ids are ignored.
    """
    chains = {}
    for r in residues:
        plain, other = chains.setdefault(r.chain, (set(), []))
        split = _split_resi(r.resi)
        if split is not None and split[1] == "":
            plain.add(split[0])
        elif r.resi not in other:
            other.append(r.resi)
    return {c: (sorted(p), o) for c, (p, o) in chains.items()}


def _pymol_num(n):
    # A leading minus would be read as a range operator.
    return "\\%d" % n if n < 0 else "%d" % n


def format_pymol(residues):
    """PyMOL selection: (chain A and resi 1-5+9) or (chain B and resi 3)."""
    parts = []
    for chain, (plain, other) in _group_by_chain(residues).items():
        items = []
        for a, b in _runs(plain):
            if a == b:
                items.append(_pymol_num(a))
            elif a + 1 == b:
                items.extend([_pymol_num(a), _pymol_num(b)])
            else:
                items.append("%s-%s" % (_pymol_num(a), _pymol_num(b)))
        items.extend(_escape_resi(x) for x in other)
        chain_sel = "chain %s" % chain if chain.strip() else "chain ''"
        parts.append("(%s and resi %s)" % (chain_sel, "+".join(items)))
    return " or ".join(parts) if parts else "none"


def format_rosetta(residues):
    """Rosetta ResidueIndexSelector PDB numbering: 1A-5A,9A,3B.

    Returns (string, skipped) where skipped lists residues that cannot be
    expressed (insertion codes, missing or multi-letter chain ids).
    """
    items, skipped = [], []
    for chain, (plain, other) in _group_by_chain(residues).items():
        label = "chain %s" % chain if chain.strip() else "no chain"
        if len(chain) != 1 or not chain.strip():
            skipped.extend("%s (%s)" % (n, label) for n in plain)
            skipped.extend("%s (%s)" % (x, label) for x in other)
            continue
        # Negative numbers would be ambiguous with the range dash.
        skipped.extend("%d (%s)" % (n, label) for n in plain if n < 0)
        for a, b in _runs([n for n in plain if n >= 0]):
            if a == b:
                items.append("%d%s" % (a, chain))
            else:
                items.append("%d%s-%d%s" % (a, chain, b, chain))
        skipped.extend("%s (%s)" % (x, label) for x in other)
    return ",".join(items), skipped


class _Res(object):
    __slots__ = ("chain", "resi")

    def __init__(self, chain, resi):
        self.chain, self.resi = chain, resi


def residues_in_selection(selection):
    """Unique (chain, resi) residues of a PyMOL selection, in model order."""
    seen, out = set(), []

    def add(chain, resi):
        if (chain, resi) not in seen:
            seen.add((chain, resi))
            out.append(_Res(chain, resi))

    cmd.iterate(selection, "add(chain, resi)", space={"add": add})
    return out


def format_selection(selection, fmt="pymol"):
    """Format the residues of any PyMOL selection.

    Returns (text, skipped). pymol: (chain A and resi 1-5) or (...);
    rosetta: 1A-5A,3B. ``skipped`` lists residues the Rosetta format
    cannot express (always empty for pymol).
    """
    fmt = fmt.lower()
    if fmt not in FORMATS:
        raise ValueError("format must be one of %s" % ", ".join(FORMATS))
    residues = residues_in_selection(selection)
    if fmt == "pymol":
        return format_pymol(residues), []
    return format_rosetta(residues)


def format_layers(layers=None, fmt="pymol"):
    """Text with one line per layer, plus a list of warnings.

    pymol:   select core, (chain A and resi 1-5) or (chain B and resi 3+5)
    rosetta: core: 1A-5A,3B,5B
    """
    if layers is None:
        layers = last_layers
    if layers is None:
        raise ValueError("no classification yet; run coretosurf first")
    fmt = fmt.lower()
    if fmt not in FORMATS:
        raise ValueError("format must be one of %s" % ", ".join(FORMATS))
    lines, warnings = [], []
    for name in LAYERS:
        if fmt == "pymol":
            lines.append("select %s, %s" % (name, format_pymol(layers[name])))
        else:
            text, skipped = format_rosetta(layers[name])
            lines.append("%s: %s" % (name, text))
            if skipped:
                warnings.append(
                    "%s: %d residue(s) with insertion codes, negative "
                    "numbers or without a single-letter chain id "
                    "omitted: %s"
                    % (name, len(skipped), ", ".join(skipped)))
    return "\n".join(lines), warnings


def apply_layers(obj, core_cut, boundary_cut, unit="abs", state=-1,
                 color=True, prefix=""):
    """Create the core/boundary/surface selections and color them.

    Returns a dict with the number of residues per layer.
    """
    global last_layers
    residues = get_residues(obj, state)
    layers = classify(residues, core_cut, boundary_cut, unit)
    last_layers = layers
    for name in LAYERS:
        sele = prefix + name
        cmd.select(sele, _selection_expr(obj, layers[name]), enable=0)
        if color:
            cmd.color(COLORS[name], "%s and elem C" % sele)
    return {name: len(layers[name]) for name in LAYERS}
