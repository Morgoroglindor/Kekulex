"""
layout2d.py - 2D coordinates for a structure that comes without any (SMILES, a name, Clean up): the skeletons chemists
always draw the same way are placed as they are usually drawn, the rest is laid out by CoordGen. New skeletons only
need a line in TEMPLATES: a SMARTS and the (x, y) of its atoms, in SMARTS order, in bond lengths with y upwards.
"""
from math import cos, hypot, pi, radians, sin, sqrt

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCoordGen, rdDepictor
from rdkit.Geometry import Point2D

R3 = sqrt(3) / 2                                            # half the width of a hexagon standing on a point
SAME_SHAPE = 0.35                                           # RMS, bond lengths: CoordGen's skeleton is the template's


def _hexagon(cx, cy):
    """The corners of a hexagon with vertical sides, clockwise from the top."""
    return [(cx, cy + 1), (cx + R3, cy + 0.5), (cx + R3, cy - 0.5), (cx, cy - 1), (cx - R3, cy - 0.5), (cx - R3, cy + 0.5)]


def _pentagon(a, b):
    """The other three corners of a regular pentagon on the side a-b, on the left going from b to a: from a's neighbour
    round to b's."""
    out, (x, y) = [], a
    dx, dy = a[0] - b[0], a[1] - b[1]
    for _ in range(3):
        c, s = cos(radians(-72)), sin(radians(-72))
        dx, dy = dx * c - dy * s, dx * s + dy * c
        x, y = x + dx, y + dy
        out.append((x, y))
    return out


# Rings drawn side by side: positions used by several templates
_L = _hexagon(0, 0)                                         # top, upper right, lower right, bottom, lower left, upper left
_R = _hexagon(2 * R3, 0)
_C5_TOP, _C5_MID, _C5_BOTTOM = (R3 + 0.951, 0.809), (R3 + 1.539, 0.0), (R3 + 0.951, -0.809)     # a 5-ring right of _L

# C1 ... C9, C11 ... C18 of a triterpene (C10 goes last in the SMARTS)
_RINGS_ABCD = [(0, 1), (-R3, 0.5), (-R3, -0.5), (0, -1), (R3, -0.5), (2 * R3, -1), (3 * R3, -0.5), (3 * R3, 0.5), (2 * R3, 1),
               (2 * R3, 2), (3 * R3, 2.5), (4 * R3, 2), (4 * R3, 1), (5 * R3, 0.5), (6 * R3, 1), (6 * R3, 2), (5 * R3, 2.5)]

TEMPLATES = (
    # Steroid: C1 ... C17, ring A lower left, ring D upper right
    ('steroid', '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](~[#6]~[#6]~[#6]4~[#6]3~[#6]~[#6]~[#6]4)~[#6]12',
     [(0, 1), (-R3, 0.5), (-R3, -0.5), (0, -1), (R3, -0.5), (2 * R3, -1), (3 * R3, -0.5), (3 * R3, 0.5), (2 * R3, 1),
      (2 * R3, 2), (3 * R3, 2.5), (4 * R3, 2), (4 * R3, 1), (4 * R3 + 0.951, 0.691), (4 * R3 + 1.539, 1.5),
      (4 * R3 + 0.951, 2.309), (R3, 0.5)]),
    # Pentacyclic triterpenes, rising from ring A at the lower left to ring E at the upper right, C1 ... C22.
    # Oleanane and ursane (also friedelane, taraxastane): five six-membered rings
    ('oleanane', '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](~[#6]~[#6]~[#6]4~[#6]3~[#6]~[#6]~[#6]5~[#6]4~[#6]~[#6]~[#6]~[#6]5)~[#6]12',
     _RINGS_ABCD + [(5 * R3, 3.5), (6 * R3, 4), (7 * R3, 3.5), (7 * R3, 2.5), (R3, 0.5)]),
    # Lupane: ring E has five members
    ('lupane', '[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](~[#6]~[#6]~[#6]4~[#6]3~[#6]~[#6]~[#6]5~[#6]4~[#6]~[#6]~[#6]5)~[#6]12',
     _RINGS_ABCD + _pentagon((5 * R3, 2.5), (6 * R3, 2)) + [(R3, 0.5)]),
    # Morphinan with its 4,5-epoxy bridge: aromatic ring upper left, N on the right. C1, C2, C3, C4, C12, C11, C10, C9, C14,
    # C13, C5, O, C6, C7, C8, C15, C16, N
    ('morphinan', '[#6]1~[#6]~[#6]~[#6]2~[#6]3~[#6]~1~[#6]~[#6]4~[#6]5~[#6]~3(~[#6](~[#8]~2)~[#6]~[#6]~[#6]~5)~[#6]~[#6]~[#7]~4',
     [(0, 1), (-R3, 0.5), (-R3, -0.5), (0, -1), (R3, -0.5), (R3, 0.5), (2 * R3, 1), (3 * R3, 0.5), (3 * R3, -0.5), (2 * R3, -1),
      (2 * R3, -2), (0.75, -2.0), (3 * R3, -2.5), (4 * R3, -2), (4 * R3, -1), (3 * R3 + 0.35, -1.0), (4 * R3 + 0.6, -0.55),
      (4 * R3 + 0.3, 0.4)]),
    # Penicillin: beta-lactam on the left, thiazolidine on the right with its carboxyl at the lower right
    ('penam', '[#6]1~[#6](=[#8])~[#7]2~[#6]~1~[#16]~[#6]~[#6]~2',
     [(0, 0.5), (0, -0.5), (-0.7, -1.2), (1, -0.5), (1, 0.5), (1.95, 0.81), (2.54, 0), (1.95, -0.81)]),
    # Cephalosporin: the same with a six-membered ring
    ('cephem', '[#6]1~[#6](=[#8])~[#7]2~[#6]~1~[#16]~[#6]~[#6]~[#6]~2',
     [(0, 0.5), (0, -0.5), (-0.7, -1.2), (1, -0.5), (1, 0.5), (1 + R3, 1), (1 + 2 * R3, 0.5), (1 + 2 * R3, -0.5), (1 + R3, -1)]),
    # Purine (and xanthines): six-membered ring on the left, N9 at the lower right; N1, C2, N3, C4, C5, C6, N7, C8, N9
    ('purine', '[#7]1~[#6]~[#7]~[#6]2~[#6](~[#6]~1)~[#7]~[#6]~[#7]~2',
     [_L[5], _L[4], _L[3], _L[2], _L[1], _L[0], _C5_TOP, _C5_MID, _C5_BOTTOM]),
    # Betalain core (betalamic acid): ring N at the bottom with a carboxyl on each side, the conjugated chain going up.
    # N, C2 (the unsaturated side, on the left), its carboxyl C, C3, C4, chain C, C5, C6, its carboxyl C
    ('betalamic', '[#7;R]1~[#6;R;$(*=*),a](-[#6](~[#8])~[#8])~[#6;R]~[#6;R](~[#6;!R])~[#6;R]~[#6;R]~1-[#6](~[#8])~[#8]',
     [(0, -1), (-R3, -0.5), (-2 * R3, -1), (-2 * R3 - R3, -0.5), (-2 * R3, -2), (-R3, 0.5), (0, 1), (0, 2), (R3, 0.5), (R3, -0.5),
      (2 * R3, -1), (3 * R3, -0.5), (2 * R3, -2)]),
    # Flavonoid core: chromene on the left with O at the bottom, the 2-aryl ring on the right
    ('flavone', '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6](~[#8]~2)-[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1',
     [_L[5], _L[4], _L[3], _L[2], _L[1], _L[0], _R[0], _R[1], _R[2], _R[3],
      (4 * R3, -1), (5 * R3, -0.5), (6 * R3, -1), (6 * R3, -2), (5 * R3, -2.5), (4 * R3, -2)]),
    # Quinuclidine: the usual cage
    ('quinuclidine', '[#7]12~[#6]~[#6]~[#6](~[#6]~[#6]~1)~[#6]~[#6]~2',
     [(0, 1.1), (-0.95, 0.55), (-0.95, -0.55), (0, -1.1), (0.95, -0.55), (0.95, 0.55), (0.3, -0.25), (0.3, 0.55)]),
    # Pyranose: ring O at the upper right, C5 with its CH2OH at the upper left
    ('pyranose', '[#6;!R]-[#6;R]1-[#8;R]-[#6;R]-[#6;R]-[#6;R]-[#6;R]-1',
     [(-1.73, 1.0), (-R3, 0.5), (0, 1), (R3, 0.5), (R3, -0.5), (0, -1), (-R3, -0.5)]),
    # Furanose with the base on the right, as in nucleosides: ring O on top
    ('furanose', '[#6;!R]-[#6;R]1-[#8;R]-[#6;R]-[#6;R]-[#6;R]-1',
     [(-1.75, 0.75), (-0.95, 0.31), (0, 0.95), (0.95, 0.31), (0.59, -0.81), (-0.59, -0.81)]),
    # Indole: benzene on the left, NH at the lower right
    ('indole', '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#7]~2',
     [_L[5], _L[4], _L[3], _L[2], _L[1], _L[0], _C5_TOP, _C5_MID, _C5_BOTTOM]),
    # Quinoline: N at the bottom of the right-hand ring
    ('quinoline', '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]~[#7]~2',
     [_L[5], _L[4], _L[3], _L[2], _L[1], _L[0], _R[0], _R[1], _R[2], _R[3]]),
    # Isoquinoline: N on the right
    ('isoquinoline', '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#7]~[#6]~2',
     [_L[5], _L[4], _L[3], _L[2], _L[1], _L[0], _R[0], _R[1], _R[2], _R[3]]),
)
_NEUTRAL_BETALAIN = Chem.MolFromSmarts('[#7;R;X2]=[#6;R]-[#6;R]=[#6;R]-[#6;!R]=[#6;!R]-[#7;X3;+0]')
_IMINIUM_ACID = Chem.MolFromSmarts('[#6;!R]=[#7+;R](-[#6;R]-[#6](~[#8])~[#8])')
_ACID = Chem.MolFromSmarts('[#6]-[#6](=O)-[OX2H1]')
_compiled = [(name, Chem.MolFromSmarts(smarts), coords) for name, smarts, coords in TEMPLATES]


def template_for(mol):
    """(name, {atom index: Point2D}) of the largest skeleton in mol, or None."""
    for name, query, coords in sorted(_compiled, key=lambda t: -len(t[2])):
        matches = mol.GetSubstructMatches(query, uniquify=False, maxMatches=64)
        if matches:
            match = max(matches, key=lambda m: _ring_a_score(mol, m)) if name == 'oleanane' else matches[0]
            return name, {a: Point2D(*coords[k]) for k, a in enumerate(match)}
    return None


def _ring_a_score(mol, match):
    """The oleanane skeleton reads the same from both ends: ring A is the end with an oxygen on C3 and methyls on C4."""
    c3, c4 = mol.GetAtomWithIdx(match[2]), mol.GetAtomWithIdx(match[3])
    return (2 * any(n.GetAtomicNum() == 8 for n in c3.GetNeighbors())
            + sum(n.GetAtomicNum() == 6 and not n.IsInRing() for n in c4.GetNeighbors()))


def usual_form(smiles):
    """A betalain as PubChem may give it (ring N=C, chain amine, as for betanin) in the form it is usually drawn: ring NH,
    iminium N+ on the chain and a carboxylate. Same formula and stereocentres; other molecules come back unchanged."""
    mol = Chem.MolFromSmiles(smiles)
    match = mol.GetSubstructMatch(_NEUTRAL_BETALAIN) if mol is not None else ()
    if not match:
        return smiles
    n1, c2, c3, c4, cc, cb, amine = match
    extended = mol.GetBondBetweenAtoms(cc, cb).GetStereo() == Chem.BondStereo.STEREOE
    rw = Chem.RWMol(mol)
    Chem.Kekulize(rw, clearAromaticFlags=True)
    for a, b, order in ((n1, c2, Chem.BondType.SINGLE), (c2, c3, Chem.BondType.DOUBLE), (c3, c4, Chem.BondType.SINGLE),
                        (c4, cc, Chem.BondType.DOUBLE), (cc, cb, Chem.BondType.SINGLE), (cb, amine, Chem.BondType.DOUBLE)):
        bond = rw.GetBondBetweenAtoms(a, b)
        bond.SetBondType(order)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
    for bond in rw.GetBonds():                              # the old E/Z marks no longer fit these double bonds
        bond.SetBondDir(Chem.BondDir.NONE)
    rw.GetAtomWithIdx(amine).SetFormalCharge(1)
    near = {n.GetIdx() for n in rw.GetAtomWithIdx(amine).GetNeighbors()}
    acids = [m for m in rw.GetSubstructMatches(_ACID) if m[0] == c2 or m[0] in near]
    acids.sort(key=lambda m: m[0] == c2)      # next to the iminium N if there is one, else the ring's
    if acids:
        oxygen = rw.GetAtomWithIdx(acids[0][3])
        oxygen.SetFormalCharge(-1)
        oxygen.SetNoImplicit(False)
    rw.GetAtomWithIdx(n1).SetNoImplicit(False)
    rw.GetAtomWithIdx(n1).SetNumExplicitHs(0)
    try:
        out = rw.GetMol()
        Chem.SanitizeMol(out)
    except (RuntimeError, ValueError):
        return smiles
    if extended:                                            # an E chain stays extended, as ChemSpider has the iminium forms
        Chem.AssignStereochemistry(out, cleanIt=True, force=True)
        bond = out.GetBondBetweenAtoms(c4, cc)
        bond.SetStereoAtoms(c3, cb)
        bond.SetStereo(Chem.BondStereo.STEREOTRANS)
        acid_side = {m[0] for m in acids if m[0] in near}
        other = [a for a in near if a != cb and a not in acid_side]
        if len(other) == 1:
            bond = out.GetBondBetweenAtoms(cb, amine)
            bond.SetStereoAtoms(cc, other[0])
            bond.SetStereo(Chem.BondStereo.STEREOTRANS)
    return Chem.MolToSmiles(out)


def _coordgen(mol):
    rdCoordGen.AddCoords(mol, rdCoordGen.CoordGenParams())


def _face_rings_left(mol):
    """Without a template: a molecule with rings has them on its left, an acid or other end group of a chain on its
    right. Only mirrors (the stereo goes on as wedges after this)."""
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    xs = [conf.GetAtomPosition(i).x for i in range(n)]
    mean = sum(xs) / n
    info = mol.GetRingInfo()
    rings = [i for i in range(n) if info.NumAtomRings(i)]
    if rings:
        off = sum(xs[i] for i in rings) / len(rings) - mean
    else:
        hetero = [i for i in range(n) if mol.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6)]
        off = (mean - sum(xs[i] for i in hetero) / len(hetero)) if hetero else 0.0
    if off > 0.5:
        for i in range(n):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (-p.x, p.y, p.z))


def _align(mol, fixed):
    """Turns (and if needed mirrors) the layout as a whole so that the skeleton lies as the template has it; returns how
    far its atoms stay from the template, RMS in bond lengths."""
    conf = mol.GetConformer()
    atoms = list(fixed)
    have = np.array([[conf.GetAtomPosition(a).x, conf.GetAtomPosition(a).y] for a in atoms])
    want = np.array([[fixed[a].x, fixed[a].y] for a in atoms])
    bond = np.median([_length(conf, b) for b in mol.GetBonds()])
    scale = 1.0 / (bond or 1.0)
    c_have, c_want = have.mean(axis=0), want.mean(axis=0)
    u, _s, vt = np.linalg.svd((have - c_have).T @ (want - c_want))
    r = vt.T @ u.T                                          # may be a mirror image: 2D coordinates carry no stereo yet
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        x, y = r @ ((np.array([p.x, p.y]) - c_have) * scale) + c_want
        conf.SetAtomPosition(i, (float(x), float(y), 0.0))
    moved = (r @ ((have - c_have) * scale).T).T + c_want
    return float(np.sqrt(((moved - want) ** 2).sum(axis=1).mean()))


def _length(conf, b):
    p, q = conf.GetAtomPosition(b.GetBeginAtomIdx()), conf.GetAtomPosition(b.GetEndAtomIdx())
    return hypot(p.x - q.x, p.y - q.y)


def clashes(mol):
    """Pairs of atoms not bonded to each other and closer than 0.45 bond lengths: an unreadable layout."""
    conf = mol.GetConformer()
    xy = conf.GetPositions()[:, :2]
    bond = np.median([_length(conf, b) for b in mol.GetBonds()]) if mol.GetNumBonds() else 1.0
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2) / (bond or 1.0)
    return sum(1 for i in range(len(xy)) for j in range(i + 1, len(xy))
               if d[i, j] < 0.45 and mol.GetBondBetweenAtoms(i, j) is None)


def compute(mol):
    """Gives mol one new 2D conformer. Its atoms and stereo stay as they are; wedges are for the caller to set.

    CoordGen lays out the whole molecule, then it is turned to lie like the template. Only when its skeleton has another
    shape (bridged systems) does the template's own shape go in, with RDKit placing the rest, if that stays readable."""
    mol.RemoveAllConformers()
    try:
        _coordgen(mol)
    except (RuntimeError, ValueError):
        rdDepictor.SetPreferCoordGen(False)
        rdDepictor.Compute2DCoords(mol)
    found = template_for(mol)
    if not found:
        _face_rings_left(mol)
        regularise(mol)
        _keep_double_bond_geometry(mol)
        return
    if _align(mol, found[1]) < SAME_SHAPE:
        conf = mol.GetConformer()
        for a, p in found[1].items():                       # exactly the template's shape, whatever CoordGen made of it
            conf.SetAtomPosition(a, (p.x, p.y, 0.0))
    else:
        trial = Chem.Mol(mol)
        rdDepictor.SetPreferCoordGen(False)
        rdDepictor.Compute2DCoords(trial, coordMap=found[1])
        if clashes(trial) <= clashes(mol):
            mol.RemoveAllConformers()
            mol.AddConformer(trial.GetConformer(), assignId=True)
    regularise(mol, set(found[1]))
    _keep_double_bond_geometry(mol, set(found[1]))
    if found[0] == 'betalamic':
        _betalain_carboxylate_by_the_chain(mol)
        _betalain_turned_over(mol)


def _keep_double_bond_geometry(mol, fixed=()):
    """A double bond with a stated E/Z that the layout drew the other way: one end turned over the bond's axis - the end
    away from the template's skeleton, else the smaller one."""
    conf = mol.GetConformer()
    trans = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
    for bond in mol.GetBonds():
        ends = list(bond.GetStereoAtoms())
        if bond.GetStereo() == Chem.BondStereo.STEREONONE or len(ends) != 2 or bond.IsInRing():
            continue
        if bond.GetStereo() not in trans + (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS):
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        xy = conf.GetPositions()[:, :2]
        axis = xy[b] - xy[a]
        axis = axis / (np.hypot(*axis) or 1.0)
        side = lambda atom, base: axis[0] * (xy[atom] - xy[base])[1] - axis[1] * (xy[atom] - xy[base])[0]
        drawn_trans = side(ends[0], a) * side(ends[1], b) < 0
        if drawn_trans == (bond.GetStereo() in trans):
            continue
        sides = [(_branch(mol, b, a), b), (_branch(mol, a, b), a)]
        moved, pivot = min(sides, key=lambda s: (len(s[0] & set(fixed)), len(s[0])))
        for i in moved:
            rel = xy[i] - xy[pivot]
            x, y = xy[pivot] + 2 * (rel @ axis) * axis - rel
            conf.SetAtomPosition(i, (float(x), float(y), 0.0))


def _betalain_turned_over(mol):
    """Betalains: turned 180 degrees about the page's vertical axis, the ring N still at the bottom and the ring's aliphatic
    carboxyl going into the page (the wedges go on after this, so they follow)."""
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (-p.x, p.y, 0.0))


def _betalain_carboxylate_by_the_chain(mol):
    """Betanin, indicaxanthin ...: the carboxyl next to the iminium N is drawn on the chain's side of the C=N+ bond, as
    in betanin (E). Turns that end over the bond's axis, unless the molecule states the geometry of that double bond."""
    match = mol.GetSubstructMatch(_IMINIUM_ACID)
    if not match:
        return
    chain, nitrogen, acid_side = match[0], match[1], match[2]
    bond = mol.GetBondBetweenAtoms(chain, nitrogen)
    if bond.GetStereo() != Chem.BondStereo.STEREONONE:
        return
    before = [n.GetIdx() for n in mol.GetAtomWithIdx(chain).GetNeighbors() if n.GetIdx() != nitrogen and n.GetAtomicNum() > 1]
    if not before:
        return
    conf = mol.GetConformer()
    xy = conf.GetPositions()[:, :2]
    axis = xy[nitrogen] - xy[chain]
    axis = axis / (np.hypot(*axis) or 1.0)
    side = lambda atom, base: axis[0] * (xy[atom] - xy[base])[1] - axis[1] * (xy[atom] - xy[base])[0]
    if side(before[0], chain) * side(acid_side, nitrogen) > 0:
        return
    for a in _branch(mol, nitrogen, chain):
        rel = xy[a] - xy[nitrogen]
        x, y = xy[nitrogen] + 2 * (rel @ axis) * axis - rel
        conf.SetAtomPosition(a, (float(x), float(y), 0.0))


def _bridged_atoms(mol):
    """Atoms of rings sharing more than one bond with another ring (cages): no regular polygons fit those."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    out = set()
    for i, r in enumerate(rings):
        for q in rings[i + 1:]:
            if len(r & q) > 2:
                out |= r | q
    return out


def _angle_targets(mol, i, xy, ring_size):
    """(neighbour, next neighbour, chord length) around atom i: the distances between neighbours that give its ideal
    angles. Inside a ring the polygon's angle, the rest shared out evenly."""
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()]
    if len(nbrs) < 2:
        return []
    nbrs.sort(key=lambda n: np.arctan2(*(xy[n] - xy[i])[::-1]))
    pairs = [(nbrs[k], nbrs[(k + 1) % len(nbrs)]) for k in range(len(nbrs))]
    fixed, free = {}, []
    for pair in pairs:
        size = ring_size(i, *pair)
        if size:
            fixed[pair] = (size - 2) * pi / size
        else:
            free.append(pair)
    if mol.GetAtomWithIdx(i).IsInRing():                    # what hangs on a ring is placed afterwards (_place_substituents)
        return [(p, q, 2 * sin(t / 2)) for (p, q), t in fixed.items()]
    if len(nbrs) == 2 and not fixed:                        # a chain: 120 degrees on the side it already bends to, a triple bond straight
        linear = any(b.GetBondType() == Chem.BondType.TRIPLE for b in mol.GetAtomWithIdx(i).GetBonds())
        return [(nbrs[0], nbrs[1], 2.0 if linear else 2 * sin(pi / 3))]
    if free:
        share = (2 * pi - sum(fixed.values())) / len(free)
        fixed.update({pair: share for pair in free})
    return [(p, q, 2 * sin(t / 2)) for (p, q), t in fixed.items() if t < 2.97]


def regularise(mol, fixed=(), rounds=400):
    """Every bond the same length and the usual angles, keeping the layout's shape: repeatedly moves pairs of atoms to the
    distances wanted. Atoms in fixed stay put; so do cages. Given up if atoms would come to overlap."""
    if not mol.GetNumBonds():
        return
    conf = mol.GetConformer()
    before = clashes(mol)
    xy = conf.GetPositions()[:, :2].copy()
    original = xy.copy()
    lengths = [np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    if not fixed:
        xy /= float(np.median(lengths)) or 1.0
    frozen = set(fixed) | _bridged_atoms(mol)
    move = np.array([0.0 if i in frozen else 1.0 for i in range(len(xy))])
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]

    def ring_size(i, p, q):
        sizes = [len(r) for r in rings if {i, p, q} <= r]
        return min(sizes) if sizes and min(sizes) <= 8 else 0
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), 1.0) for b in mol.GetBonds()]
    for i in range(len(xy)):
        pairs += _angle_targets(mol, i, xy, ring_size)
    for _ in range(rounds):
        for p, q, want in pairs:
            total = move[p] + move[q]
            if not total:
                continue
            d = xy[q] - xy[p]
            now = float(np.hypot(*d))
            if now < 1e-6:
                continue
            shift = d * (now - want) / now * 0.5
            xy[p] += shift * move[p] * 2 / total
            xy[q] -= shift * move[q] * 2 / total
    _place_substituents(mol, xy, frozen, rings)
    for i, (x, y) in enumerate(xy):
        conf.SetAtomPosition(i, (float(x), float(y), 0.0))
    if clashes(mol) > before:
        for i, (x, y) in enumerate(original):
            conf.SetAtomPosition(i, (float(x), float(y), 0.0))


def _branch(mol, start, away_from):
    """The atoms reached from start without passing through away_from."""
    seen, stack = {start}, [start]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetIdx() != away_from and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return seen


def _place_substituents(mol, xy, frozen, rings):
    """What hangs on a ring goes straight into the widest gap between the atom's other bonds, one bond length away, turned
    as a whole. Works outwards from the main ring system, so that a branch never moves what it hangs on."""
    if not rings:
        return
    start = frozen or max(rings, key=len)
    seen = set(start)
    queue = list(start)
    while queue:
        i = queue.pop(0)
        atom = mol.GetAtomWithIdx(i)
        taken = []
        children = []
        for b in atom.GetBonds():
            n = b.GetOtherAtomIdx(i)
            if b.IsInRing() or n in seen:
                taken.append(n)
                if n not in seen:
                    seen.add(n)
                    queue.append(n)
            else:
                children.append(n)
        if not atom.IsInRing():
            for n in children:
                seen.add(n)
                queue.append(n)
            continue
        if not children or not taken:
            continue
        angle = lambda n: float(np.arctan2(*(xy[n] - xy[i])[::-1])) % (2 * pi)
        have = sorted(angle(n) for n in taken)
        gaps = [((have[(k + 1) % len(have)] - have[k]) % (2 * pi) or 2 * pi, have[k]) for k in range(len(have))]
        gap, low = max(gaps)
        step = min(gap / (len(children) + 1), radians(70))
        middle = low + gap / 2
        # keep their turning order: the first child by angle inside the gap gets the first place
        children.sort(key=lambda n: (angle(n) - low) % (2 * pi))
        for k, child in enumerate(children):
            target = middle + (k - (len(children) - 1) / 2) * step
            turn = target - angle(child)
            c, s_ = cos(turn), sin(turn)
            scale = 1.0 / (float(np.hypot(*(xy[child] - xy[i]))) or 1.0)
            members = sorted(_branch(mol, child, i))
            rel = xy[members] - xy[child]
            first = xy[child] - xy[i]
            new_child = xy[i] + np.array([first[0] * c - first[1] * s_, first[0] * s_ + first[1] * c]) * scale
            xy[members] = new_child + np.column_stack((rel[:, 0] * c - rel[:, 1] * s_, rel[:, 0] * s_ + rel[:, 1] * c))
            seen.add(child)
            queue.append(child)


def with_fusion_hydrogens(mol):
    """Stereocentres at ring fusions show their H, wedged, as chemists draw them: the wedge then need not lie on a ring
    bond. Returns the molecule with those hydrogens added."""
    Chem.AssignStereochemistry(mol, cleanIt=False, force=False)
    atoms = [a.GetIdx() for a in mol.GetAtoms()
             if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and a.GetTotalNumHs() == 1 and a.GetDegree() == 3
             and all(b.IsInRing() for b in a.GetBonds())]
    return Chem.AddHs(mol, onlyOnAtoms=atoms) if atoms else mol
