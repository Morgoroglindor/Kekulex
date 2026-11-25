"""
fragments.py - the mass-spectrum fragments a structure can logically give, from the textbook fragmentation rules:

    pre, frags = fragments.predict(doc, atom_ids, mode='ESI+', steps=2)

mode: 'EI' (the radical cation M•+), 'ESI+' ([M+H]+) or 'ESI-' ([M−H]−). A structure drawn as an ion is used as it is.

EI (odd-electron M•+): α-cleavage next to N, O, S (iminium / oxonium ions) and next to C=O (acylium ions), benzylic and
allylic cleavage (C7H7+ tropylium), inductive cleavage at C–halogen and C–O, loss of HX, the McLafferty rearrangement
(the γ-H moves to the C=O oxygen), retro-Diels–Alder of a six-membered ring with a double bond, and the usual neutral
losses (H2O, CH3•, ketene, NO2•, CO, HCN ...).
ESI (even-electron [M+H]+ / [M−H]−): heterolytic cleavage of C–O, C–N and C–S bonds, either the carbon side as a cation
(acylium / b-type, benzylic, sugar oxonium) or the heteroatom side with the charge and an H moved to it (aglycone Y ions,
y-type amines, [acid+H]+, carboxylates, quinate ...), retro-Diels–Alder (flavonoid 1,3A ions), and neutral losses
(H2O, NH3, CH3OH, HCOOH, CO2, CO, C4H8, ketene, SO3, H3PO4, CH3• from aryl methoxy groups, C-glycoside cross-ring 90 / 120).
With steps > 1 each fragment loses the neutral molecules its own atoms allow (and an acylium ion CO).

Electron bookkeeping: a fragment is a set of atoms, the bonds to the rest broken and capped with H, then marks on some
atoms (charge, unpaired electron, H gained or lost); its formula and m/z follow, and the lost neutral is the precursor
minus the fragment. The same marks draw it (draw_doc), so the drawn fragment has the listed formula.
Likelihood is a rough plausibility from the rule (and the stability of the ion it gives), not a predicted intensity.
"""
from dataclasses import dataclass, field

from rdkit import Chem

import chem

MODES = (('ESI+', 'ESI+  [M+H]+'), ('ESI-', 'ESI−  [M−H]−'), ('EI', 'EI  M•+'))
LIKELY, POSSIBLE = 0.65, 0.40


@dataclass
class Fragment:
    mz: float
    z: int
    odd: bool                                               # odd-electron (radical) ion
    counts: dict
    ion: str                                                # 'C7H7+', 'C3H6O•+'
    lost: str                                               # 'C3H7•', 'H2O + CO' (from the precursor)
    how: str
    score: float
    atoms: frozenset = frozenset()                          # drawing atom ids in the ion
    lost_atoms: frozenset = frozenset()                     # drawing atom ids that leave
    step: int = 1
    also: list = field(default_factory=list)                # other rules giving the same ion
    S: frozenset = None                                     # expanded-drawing atom ids (None: formula only)
    marks: dict = None                                      # {expanded atom id: (charge, radicals, H)}

    @property
    def likelihood(self):
        return 'likely' if self.score >= LIKELY else 'possible' if self.score >= POSSIBLE else 'minor'

    @property
    def formula(self):
        return chem.formula_text(self.counts)


def ion_text(counts, z, odd):
    return chem.formula_text(counts) + ('•' if odd else '') + chem.charge_text(z)


def _sub(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) - v
    return out


def _clean(c):
    return {k: v for k, v in c.items() if v}


# (name, SMARTS, atoms that leave, attach atom, H that leave with them, radical, {mode: likelihood})
# formula-only losses (no atoms): (name, test, formula, radical, {mode: likelihood})
GROUP_LOSSES = (
    ('H2O', '[CX4;!$(C=O)][OX2H1]', (1,), 0, 2, False, {'EI': 0.45, 'ESI+': 0.75, 'ESI-': 0.35}, 'alcohol'),
    ('H2O', '[CX3](=O)[OX2H1]', (2,), 0, 2, False, {'ESI+': 0.55, 'ESI-': 0.25}, 'carboxylic acid'),
    ('CO2', '[#6][CX3](=O)[OX2H1,OX1-]', (1, 2, 3), 0, 0, False, {'EI': 0.15, 'ESI+': 0.2, 'ESI-': 0.9},
     'carboxylic acid'),
    ('HCOOH', '[#6][CX3](=O)[OX2H1]', (1, 2, 3), 0, 2, False, {'ESI+': 0.45, 'EI': 0.15}, 'carboxylic acid'),
    ('NH3', '[CX4][NX3H2;!$(N[#6]=[O,S,N])]', (1,), 0, 3, False, {'EI': 0.2, 'ESI+': 0.75}, 'primary amine'),
    ('NH3', '[CX3](=O)[NX3H2]', (2,), 0, 3, False, {'ESI+': 0.5}, 'primary amide'),
    ('CH3OH', '[CX3](=O)[OX2][CH3]', (2, 3), 0, 4, False, {'ESI+': 0.65, 'ESI-': 0.3}, 'methyl ester'),
    ('CH3OH', '[CX4][OX2][CH3]', (1, 2), 0, 4, False, {'ESI+': 0.4}, 'methyl ether'),
    ('CH3•', 'c[OX2][CH3]', (2,), 1, 3, True, {'EI': 0.55, 'ESI+': 0.55, 'ESI-': 0.6}, 'aryl methoxy group'),
    ('C4H8', '[O,N,S]C([CH3])([CH3])[CH3]', (1, 2, 3, 4), 0, 8, False, {'EI': 0.5, 'ESI+': 0.8, 'ESI-': 0.3},
     'tert-butyl group'),
    ('C2H2O', '[O,N;!$([O,N]C(=O)C(=O))]C(=O)[CH3]', (1, 2, 3), 0, 2, False, {'EI': 0.75, 'ESI+': 0.7, 'ESI-': 0.4},
     'acetyl group, as ketene'),
    ('HCl', '[CX4;!H0,$(C[C;!H0])][Cl]', (1,), 0, 1, False, {'EI': 0.45, 'ESI+': 0.4, 'ESI-': 0.2}, 'alkyl chloride'),
    ('HBr', '[CX4;!H0,$(C[C;!H0])][Br]', (1,), 0, 1, False, {'EI': 0.3, 'ESI+': 0.3}, 'alkyl bromide'),
    ('HF', '[CX4;!H0,$(C[C;!H0])][F]', (1,), 0, 1, False, {'EI': 0.35, 'ESI+': 0.4}, 'alkyl fluoride'),
    ('NO2•', '[#6][$([NX3](=O)=O),$([NX3+](=O)[O-])](~[OX1])~[OX1]', (1, 2, 3), 0, 0, True,
     {'EI': 0.6, 'ESI-': 0.5}, 'nitro group'),
    ('H2S', '[CX4][SX2H1]', (1,), 0, 2, False, {'EI': 0.35, 'ESI+': 0.6}, 'thiol'),
    ('SO3', '[O,N,#6]S(=O)(=O)[OX2H1,OX1-]', (1, 2, 3, 4), 0, 0, False, {'ESI+': 0.6, 'ESI-': 0.85},
     'sulfate / sulfonic acid'),
    ('H3PO4', '[#6]O[PX4](=O)([OX2H1,OX1-])[OX2H1,OX1-]', (1, 2, 3, 4, 5), 0, 3, False, {'ESI+': 0.8, 'ESI-': 0.4},
     'phosphate ester'),
    ('HPO3', '[#6]O[PX4](=O)([OX2H1,OX1-])[OX2H1,OX1-]', (2, 3, 4, 5), 1, 1, False, {'ESI-': 0.6, 'ESI+': 0.3},
     'phosphate ester'),
)

FORMULA_LOSSES = (
    ('CO', '[#6;R][#6X3;R](=O)[#6;R]', {'C': 1, 'O': 1}, False, {'EI': 0.45, 'ESI+': 0.55, 'ESI-': 0.4},
     'ring C=O (from the ring)'),
    ('CO2', '[#6;R](=O)[O,o;R]', {'C': 1, 'O': 2}, False, {'EI': 0.45, 'ESI+': 0.35, 'ESI-': 0.5},
     'lactone / coumarin (from the ring)'),
    ('CO2', 'O=[#6;R]1[#6;R]~[#6;R][O,o;R][#6;R]~[#6;R]1', {'C': 1, 'O': 2}, False, {'ESI+': 0.2, 'ESI-': 0.4},
     'chromone C ring (from the ring)'),
    ('C2H2O', 'O=[#6;R]1[#6;R]~[#6;R][O,o;R][#6;R]~[#6;R]1', {'C': 2, 'H': 2, 'O': 1}, False, {'ESI+': 0.35},
     'chromone C ring, as ketene (from the ring)'),
    ('CO', 'c[OX2H1]', {'C': 1, 'O': 1}, False, {'EI': 0.4, 'ESI+': 0.2, 'ESI-': 0.3}, 'phenol (from the ring)'),
    ('NO•', '[#6][$([NX3](=O)=O),$([NX3+](=O)[O-])]', {'N': 1, 'O': 1}, True, {'EI': 0.4, 'ESI-': 0.35},
     'nitro group (as the nitrite)'),
    ('HCN', '[n,$([NX3H2]c)]', {'C': 1, 'N': 1, 'H': 1}, False, {'EI': 0.35, 'ESI+': 0.25}, 'aromatic N'),
    ('C3H6O3', '[c][CX4;R1]1[OX2;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1', {'C': 3, 'H': 6, 'O': 3}, False,
     {'ESI+': 0.55, 'ESI-': 0.7}, 'C-glycoside (cross-ring, 0,3X)'),
    ('C4H8O4', '[c][CX4;R1]1[OX2;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1', {'C': 4, 'H': 8, 'O': 4},
     False, {'ESI+': 0.55, 'ESI-': 0.75}, 'C-glycoside (cross-ring, 0,2X)'),
)


class _Mol:
    """The drawing's structure as an RDKit molecule, its atoms mapped back to the drawing."""

    def __init__(self, doc, atom_ids=None):
        full, grown = doc.expanded()
        self.full = full
        self.grown = grown
        self.owner = {a: a for a in doc.atoms}
        for a, more in grown.items():
            for g in more:
                self.owner[g] = a
        ids = [a for a in doc.atoms if atom_ids is None or a in atom_ids]
        ids = {a for a in chem.expand_ids(grown, ids) if a in full.atoms}
        self.mol, self.order, self.error = chem.to_mol(full, ids)
        if self.mol is None:
            return
        m = self.mol
        self.n = m.GetNumAtoms()
        self.el = [a.GetSymbol() for a in m.GetAtoms()]
        self.key = []
        for a in m.GetAtoms():
            iso = a.GetIsotope()
            self.key.append(chem.iso_key(a.GetSymbol(), iso) if iso else a.GetSymbol())
        self.mass = [chem.key_mass(k) or chem.key_mass(e) for k, e in zip(self.key, self.el)]
        self.hs = [a.GetTotalNumHs() for a in m.GetAtoms()]
        self.q = [a.GetFormalCharge() for a in m.GetAtoms()]
        self.rad = [a.GetNumRadicalElectrons() for a in m.GetAtoms()]
        self.nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in m.GetAtoms()]
        self.kek = Chem.Mol(m)                              # bond orders as single / double (pyranone rings ...)
        try:
            Chem.Kekulize(self.kek, clearAromaticFlags=False)
        except Exception:
            self.kek = m
        # aromatic as a chemist sees it: rings of aromatic atoms without an exocyclic C=O (benzene, pyridine, furan;
        # not the pyranone ring of a flavone or coumarin, whose bonds break like ordinary ones)
        self.arobond = set()
        aro = set()
        for ring in m.GetRingInfo().AtomRings():
            if not all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                continue
            if any(self.el[j] in ('O', 'S', 'N') and j not in ring and
                   self.kek.GetBondBetweenAtoms(i, j).GetBondType() == Chem.BondType.DOUBLE
                   for i in ring for j in self.nb[i]):
                continue
            aro |= set(ring)
            for k in range(len(ring)):
                self.arobond.add(frozenset((ring[k], ring[(k + 1) % len(ring)])))
        self.aro = [i in aro for i in range(self.n)]
        self._sides = {}

    # ------------------------------------------------------------------------------------------- bonds and pieces
    def bond(self, i, j):
        return self.kek.GetBondBetweenAtoms(i, j)

    def aromatic_bond(self, i, j):
        return frozenset((i, j)) in self.arobond

    def single(self, i, j, ring=False):
        b = self.bond(i, j)
        return (b is not None and b.GetBondType() == Chem.BondType.SINGLE and (ring or not b.IsInRing())
                and not self.aromatic_bond(i, j))

    def double(self, i, j):
        b = self.bond(i, j)
        return b is not None and b.GetBondType() == Chem.BondType.DOUBLE and not self.aromatic_bond(i, j)

    def bond_order(self, i, j):
        if self.aromatic_bond(i, j):
            return None
        b = self.bond(i, j)
        t = b.GetBondType()
        return {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}.get(t)

    def h(self, i):
        return self.hs[i] + sum(1 for j in self.nb[i] if self.el[j] == 'H')

    def heavy(self, i):
        return [j for j in self.nb[i] if self.el[j] != 'H']

    def side(self, i, j, cut=()):
        """The atoms reached from i without crossing bond i-j (and the bonds in cut); None when j is reached."""
        key = (i, j, tuple(sorted(cut)))
        if key in self._sides:
            return self._sides[key]
        blocked = {frozenset((i, j))} | {frozenset(c) for c in cut}
        seen, todo = {i}, [i]
        while todo:
            x = todo.pop()
            for y in self.nb[x]:
                if y not in seen and frozenset((x, y)) not in blocked:
                    seen.add(y); todo.append(y)
        out = None if j in seen else frozenset(seen)
        self._sides[key] = out
        return out

    def cut_orders(self, i, S):
        out = 0
        for j in self.nb[i]:
            if j not in S:
                o = self.bond_order(i, j)
                if o is None:
                    return None                             # an aromatic bond: no capped piece
                out += o
        return out

    def hcap(self, S):
        total = 0
        for i in S:
            c = self.cut_orders(i, S)
            if c is None:
                return None
            total += self.hs[i] + c
        return total

    def counts(self, S, h):
        out = {}
        for i in S:
            out[self.key[i]] = out.get(self.key[i], 0) + 1
        if h:
            out['H'] = out.get('H', 0) + h
        return out

    # ------------------------------------------------------------------------------------------- atom types
    def carbonyl_o(self, c):
        """The O (or S) double-bonded to carbon c, else None."""
        if self.el[c] != 'C' or self.aro[c]:
            return None
        return next((j for j in self.nb[c] if self.el[j] in ('O', 'S') and self.double(c, j)), None)

    def saturated(self, c):
        return self.el[c] == 'C' and not self.aro[c] and all(self.single(c, j, ring=True) for j in self.nb[c])

    def vinyl(self, c):
        return any(self.el[j] == 'C' and self.double(c, j) and not self.aro[j] for j in self.nb[c])

    def anomeric(self, c, x):
        """c is a sugar's anomeric carbon (in a ring, bonded to the ring O) and x an atom outside that ring."""
        if self.el[c] != 'C' or self.aro[c]:
            return False
        ring = self.mol.GetRingInfo()
        for o in self.nb[c]:
            if o != x and self.el[o] == 'O' and self.bond(c, o).IsInRing():
                for r in ring.AtomRings():
                    if c in r and o in r and x not in r and len(r) in (5, 6):
                        return True
        return False

    def cation(self, c, S):
        """0..1: how well carbon c carries a positive charge in piece S."""
        if self.el[c] != 'C':
            return 0.3
        if self.aro[c]:
            return 0.15
        if self.carbonyl_o(c) is not None:
            return 0.85
        if not self.saturated(c):
            return 0.1                                      # vinyl cation
        near = [j for j in self.nb[c] if j in S and self.el[j] != 'H']
        if any(self.el[j] in ('N', 'O', 'S') and self.single(c, j, ring=True) and self.carbonyl_o(j) is None
               for j in near):
            return 0.7                                      # iminium / oxocarbenium
        if any(self.aro[j] for j in near):
            return 0.75                                     # benzylic
        if any(self.vinyl(j) for j in near) or self.vinyl(c):
            return 0.55                                     # allylic
        return {0: 0.05, 1: 0.15, 2: 0.35}.get(len(near), 0.6)


class _Engine:
    def __init__(self, m, mode):
        self.m = m
        self.mode = mode
        allh = sum(m.hs) + sum(1 for e in m.el if e == 'H')
        self.M = m.counts(range(m.n), sum(m.hs))
        drawn_z = sum(m.q)
        drawn_rad = sum(m.rad)
        self.charged = bool(drawn_z)
        if self.charged:                                    # an ion drawn as such: used as it is
            self.z, self.odd, self.pre = drawn_z, bool(drawn_rad % 2), dict(self.M)
            self.label = 'M' + ('•' if self.odd else '') + chem.charge_text(drawn_z)
            self.mode = 'ESI+' if drawn_z > 0 else 'ESI-'
            if self.odd:
                self.mode = 'EI'
        elif mode == 'EI':
            self.z, self.odd, self.pre, self.label = 1, not drawn_rad % 2, dict(self.M), '[M]•+'
        elif mode == 'ESI+':
            self.z, self.odd, self.label = 1, bool(drawn_rad % 2), '[M+H]+'
            self.pre = dict(self.M, H=self.M.get('H', 0) + 1)
        else:
            self.z, self.odd, self.label = -1, bool(drawn_rad % 2), '[M−H]−'
            self.pre = dict(self.M, H=self.M.get('H', 0) - 1)
        self.allh = allh
        self.out = []

    # ------------------------------------------------------------------------------------------- making fragments
    def mz(self, counts, z):
        mass = sum(chem.key_mass(k) * n for k, n in counts.items())
        return (mass - z * chem.ELECTRON) / abs(z)

    def place(self, S, base, site, dq, rad, dh):
        """Marks: base ones plus dq / rad / dh at site; H that the site cannot give comes from its neighbours."""
        m = self.m
        marks = {k: list(v) for k, v in (base or {}).items() if k in S}
        marks.setdefault(site, [0, 0, 0])
        marks[site][0] += dq

        def room(i):
            c = m.cut_orders(i, S)
            return None if c is None else m.hs[i] + c + marks.get(i, [0, 0, 0])[2]
        if rad and dh < 0:                                  # distonic: the unpaired electron on a neighbour, −H
            j = next((j for j in m.nb[site] if j in S and j not in marks and m.el[j] == 'C' and room(j)), None)
            if j is not None:
                marks[j] = [0, rad, -1]
                rad, dh = 0, dh + 1
        marks[site][1] += rad
        r = room(site)
        if r is None:
            return None
        if dh >= 0 or r + dh >= 0:
            marks[site][2] += dh
        else:
            marks[site][2] -= r
            need = -(r + dh)
            seen, todo = {site}, [site]
            while need and todo:
                x = todo.pop(0)
                for y in m.nb[x]:
                    if y in S and y not in seen:
                        seen.add(y); todo.append(y)
                        ry = room(y)
                        if ry:
                            take = min(ry, need, 1)
                            marks.setdefault(y, [0, 0, 0])[2] -= take
                            need -= take
                            if not need:
                                break
            if need:
                return None
        return {k: tuple(v) for k, v in marks.items() if any(v)}

    def acid_site(self, S, avoid=()):
        """The O–H in S most likely to lose its proton in negative mode: a carboxylic acid, then a phenol, then an
        alcohol (None when there is none)."""
        m = self.m
        best = None
        for o in S:
            if o in avoid or m.el[o] != 'O' or not m.hs[o] or len(m.heavy(o)) != 1:
                continue
            c = m.heavy(o)[0]
            rank = 0 if m.carbonyl_o(c) is not None else 1 if m.aro[c] else 2
            if best is None or rank < best[0]:
                best = (rank, o)
        return best[1] if best else None

    def base_site(self, S, avoid=()):
        """The atom in S most likely to hold the proton in positive mode: an amine N, a pyridine-type N, an aniline
        N, a C=O oxygen, another O (None when there is none)."""
        m = self.m
        best = None
        for i in S:
            if i in avoid or m.q[i] or m.el[i] not in ('N', 'O'):
                continue
            hv = m.heavy(i)
            if m.el[i] == 'N':
                if m.aro[i]:
                    rank = 1 if len(hv) == 2 and not m.hs[i] else 9
                elif any(m.carbonyl_o(j) is not None for j in hv) or not all(m.single(i, j, True) for j in m.nb[i]):
                    rank = 9
                else:
                    rank = 2 if any(m.aro[j] for j in hv) else 0
            elif any(m.double(i, j) for j in m.nb[i]):
                rank = 3
            else:
                rank = 4
            if rank < 9 and (best is None or rank < best[0]):
                best = (rank, i)
        return best[1] if best else None

    def make(self, S, marks, odd, how, score, step=1, lost_name=None, parent=None):
        """A fragment from atom set S (molecule indices) and marks, or None when the bookkeeping fails."""
        m = self.m
        S = frozenset(S)
        if not S or len(S) == m.n and not marks:
            return None
        hc = m.hcap(S)
        if hc is None or marks is None:
            return None
        h = hc + sum(v[2] for v in marks.values())
        for i, v in marks.items():
            if i not in S or m.hs[i] + m.cut_orders(i, S) + v[2] < 0:
                return None
        z = sum(m.q[i] for i in S) + sum(v[0] for v in marks.values())
        if z != self.z or h < 0:
            return None
        counts = _clean(m.counts(S, h))
        if parent is not None and (not counts.get('C') or sum(v for k, v in counts.items() if k != 'H') < 2):
            return None                                     # nothing sensible is left
        start = parent.counts if parent else self.pre
        lost = _sub(start, counts)
        if any(v < 0 for v in lost.values()) or not any(lost.values()):
            return None
        lost = _clean(lost)
        was_odd = parent.odd if parent else self.odd
        name = lost_name or (chem.formula_text(lost) + ('•' if was_odd != odd else ''))
        total = name if parent is None else f'{parent.lost} + {name}'
        return self._frag(counts, odd, total, how, score, step, S, marks)

    def formula_only(self, lost, radical, how, score, highlight, parent=None, step=1):
        start = parent.counts if parent else self.pre
        counts = _sub(start, lost)
        if any(v < 0 for v in counts.values()):
            return None
        counts = _clean(counts)
        if not counts.get('C') and not counts.get('N'):
            return None
        odd = (parent.odd if parent else self.odd) != radical
        name = chem.formula_text(lost) + ('•' if radical else '')
        total = name if parent is None else f'{parent.lost} + {name}'
        f = self._frag(counts, odd, total, how, score, step, None, None)
        f.atoms = frozenset(highlight)
        return f

    def _frag(self, counts, odd, lost, how, score, step, S, marks):
        m = self.m
        f = Fragment(self.mz(counts, self.z), self.z, odd, counts, ion_text(counts, self.z, odd), lost, how,
                     max(0.0, min(1.0, score)), step=step)
        if S is not None:
            f.S = frozenset(m.order[i] for i in S)
            f.marks = {m.order[i]: v for i, v in marks.items()}
            f._idx = S
            f._imarks = marks
            f.atoms = frozenset(m.owner[m.order[i]] for i in S)
        else:
            f._idx = None
        return f

    def add(self, f):
        if f is not None:
            self.out.append(f)
        return f

    # ------------------------------------------------------------------------------------------- the rules
    def run(self, steps=2):
        m = self.m
        if self.charged:
            self.charge_remote()
        elif self.mode == 'EI':
            self.ei_cleavages()
            self.mclafferty()
        else:
            self.esi_cleavages()
        if not self.charged:
            self.rda()
        self.group_losses()
        frontier = [f for f in self.out if f._idx is not None]
        for step in range(2, steps + 1):
            frontier = sorted(frontier, key=lambda f: -f.score)[:60]
            new = []
            for f in frontier:
                new += self.next_losses(f, step)
            self.out += new
            frontier = [f for f in new if f._idx is not None]
        return self.merged()

    def ei_cleavages(self):
        m = self.m
        alpha = {'N': 0.85, 'S': 0.7, 'O': 0.65}
        for y in range(m.n):                                # α-cleavage next to a heteroatom
            if m.el[y] not in alpha or m.aro[y] or m.q[y] > 0:
                continue
            acyl = any(m.carbonyl_o(j) is not None for j in m.nb[y])
            kind = {'N': 'iminium', 'O': 'oxonium', 'S': 'sulfonium'}[m.el[y]]
            for ca in m.nb[y]:
                if not m.single(ca, y, ring=True) or not m.saturated(ca):
                    continue
                for r in m.heavy(ca):
                    if r == y or not m.single(ca, r):
                        continue
                    S = m.side(ca, r)
                    if S is None:
                        continue
                    score = alpha[m.el[y]] - (0.3 if acyl else 0) + (0.05 if len(m.side(r, ca)) >= 2 else 0)
                    self.add(self.make(S, {ca: (1, 0, -1)}, False, f'α-cleavage at {m.el[y]} ({kind} ion)', score))
        for c in range(m.n):                                # α-cleavage at C=O: acylium ions
            o = m.carbonyl_o(c)
            if o is None:
                continue
            for r in m.heavy(c):
                if r == o or not m.single(c, r):
                    continue
                S = m.side(c, r)
                if S is None:
                    continue
                if m.el[r] == 'C':
                    score = 0.85
                elif m.el[r] == 'O':
                    score = 0.45 if m.h(r) else 0.75
                elif m.el[r] == 'N':
                    score = 0.55
                else:
                    score = 0.4
                self.add(self.make(S, {c: (1, 0, -1)}, False, 'α-cleavage at C=O (acylium ion)', score))
                if m.el[r] == 'C':
                    S2 = m.side(r, c)
                    self.add(self.make(S2, {r: (1, 0, -1)}, False, 'α-cleavage at C=O (charge on the other side)',
                                       0.6 * m.cation(r, S2)))
        for b in range(m.n):                                # benzylic / allylic cleavage
            if not m.saturated(b):
                continue
            aro = [j for j in m.nb[b] if m.aro[j]]
            allyl = [j for j in m.nb[b] if m.vinyl(j)]
            if not aro and not allyl:
                continue
            for r in m.heavy(b):
                if r in aro or r in allyl or not m.single(b, r):
                    continue
                S = m.side(b, r)
                if S is None:
                    continue
                f = self.make(S, {b: (1, 0, -1)}, False, 'benzylic cleavage' if aro else 'allylic cleavage',
                              0.9 if aro else 0.6)
                if f is not None and f.counts == {'C': 7, 'H': 7}:
                    f.how = 'benzylic cleavage (tropylium ion C7H7+)'
                self.add(f)
        base = {'I': 0.75, 'Br': 0.7, 'Cl': 0.45, 'F': 0.15}
        for x in range(m.n):                                # inductive cleavage at C–X
            for c in m.nb[x]:
                if m.el[c] != 'C' or not m.single(c, x):
                    continue
                S = m.side(c, x)
                if S is None:
                    continue
                if m.el[x] in base:
                    s = base[m.el[x]] * (0.5 + 0.5 * m.cation(c, S)) * (0.6 if m.aro[c] else 1)
                    self.add(self.make(S, {c: (1, 0, -1)}, False, f'inductive cleavage at C–{m.el[x]}', s))
                elif m.el[x] in ('O', 'S') and not m.h(x) and m.saturated(c):
                    s = 0.3 * (0.5 + m.cation(c, S))
                    self.add(self.make(S, {c: (1, 0, -1)}, False, f'inductive cleavage at C–{m.el[x]}', s))
        for a in range(m.n):                                # σ-cleavage of C–C bonds (branch points)
            if not m.saturated(a):
                continue
            for b in m.nb[a]:
                if b < a or not m.saturated(b) or not m.single(a, b):
                    continue
                for c, other in ((a, b), (b, a)):
                    S = m.side(c, other)
                    if S is None or len(S) < 2:
                        continue
                    s = m.cation(c, S)
                    if s >= 0.35:
                        self.add(self.make(S, {c: (1, 0, -1)}, False, 'σ-cleavage of C–C (alkyl cation)',
                                           0.1 + 0.4 * s))

    def mclafferty(self):
        m = self.m
        done = set()
        for c1 in range(m.n):
            ys = [y for y in m.nb[c1] if m.double(c1, y) and m.el[y] in ('O', 'N', 'S', 'C') and not m.aro[y]]
            ortho = [y for y in m.nb[c1] if m.aro[y] and m.el[y] == 'C'] if m.aro[c1] else []
            if m.el[c1] != 'C' or (not ys and not ortho):
                continue
            for xa in m.heavy(c1):
                if xa in ys or m.aro[xa] or m.el[xa] not in ('C', 'O', 'N') or not m.single(c1, xa, ring=True):
                    continue
                if m.el[xa] == 'C' and not m.saturated(xa):
                    continue
                for cb in m.heavy(xa):
                    if cb == c1 or not m.saturated(cb) or not m.single(xa, cb):
                        continue
                    S = m.side(xa, cb)
                    if S is None:
                        continue
                    gs = [g for g in m.nb[cb] if g not in S and m.el[g] == 'C' and m.h(g) and not m.aro[g]]
                    if not gs:
                        continue
                    for y in ys:
                        if (y, xa, cb) in done:
                            continue
                        done.add((y, xa, cb))
                        s = {'O': 0.85, 'S': 0.6, 'N': 0.6, 'C': 0.35}[m.el[y]] - (0.15 if m.el[xa] == 'O' else 0)
                        self.add(self.make(S, {y: (1, 1, 0)}, True, 'McLafferty rearrangement', s))
                        S2 = m.side(cb, xa)
                        self.add(self.make(S2, {cb: (1, 0, -1), gs[0]: (0, 1, -1)}, True,
                                           'McLafferty rearrangement (charge on the alkene)', 0.3))
                    if ortho and not ys and m.el[xa] == 'C':
                        self.add(self.make(S, {ortho[0]: (1, 1, 0)}, True, 'McLafferty-type H shift to the ring',
                                           0.55))

    def rda(self):
        """Retro-Diels–Alder: a six-membered ring with a double (or aromatic fusion) bond a=b splits at the two
        bonds across from it into the diene and the dienophile."""
        m = self.m
        seen = set()
        for ring in m.mol.GetRingInfo().AtomRings():
            if len(ring) != 6:
                continue
            for k in range(6):
                p, a, b, q, u, w = (ring[(k + d) % 6] for d in (-1, 0, 1, 2, 3, 4))
                if not (m.double(a, b) or m.aromatic_bond(a, b)):
                    continue
                if not (m.single(p, a, ring=True) and m.single(b, q, ring=True) and m.single(q, u, ring=True)
                        and m.single(p, w, ring=True)):
                    continue
                if m.aromatic_bond(u, w):
                    continue
                key = frozenset((frozenset((q, u)), frozenset((p, w))))
                if key in seen:
                    continue
                seen.add(key)
                D = m.side(a, -1, cut=((q, u), (p, w)))
                if D is None or u in D or w in D:
                    continue
                E = frozenset(range(m.n)) - D
                if self.mode == 'EI':
                    site = p if m.el[p] == 'C' else q
                    other = q if site == p else p
                    self.add(self.make(D, {site: (1, 0, -1), other: (0, 1, -1)}, True, 'retro-Diels–Alder (diene)',
                                       0.7))
                    self.add(self.make(E, {w: (1, 0, -1), u: (0, 1, -1)}, True, 'retro-Diels–Alder (dienophile)',
                                       0.3))
                else:
                    sign = 1 if self.z > 0 else -1
                    for piece, ends, s, part in ((D, (p, q), 0.65, '1,3A-type diene'),
                                                 (E, (w, u), 0.35, 'dienophile')):
                        if sign > 0:            # the proton on the piece: one end keeps its cap H, the other + and −H
                            c = sorted(ends, key=lambda i: (m.el[i] != 'C', m.carbonyl_o(i) is None))[0]
                            marks = self.place(piece, {}, c, 1, 0, -1)
                        else:                   # an acidic OH gives the proton; the two ends close the diene
                            acid = self.acid_site(piece, ends)
                            if acid is not None:
                                marks = self.place(piece, {}, acid, -1, 0, -1)
                                for e in ends:
                                    marks = self.place(piece, marks, e, 0, 0, -1) if marks else None
                            else:
                                marks = None            # nothing to deprotonate: no [M−H]− ion of this piece
                        self.add(self.make(piece, marks, False, f'retro-Diels–Alder ({part})', s))

    def esi_cleavages(self):
        """Heterolysis of C–O, C–N, C–S bonds outside rings, as in CID of [M+H]+ / [M−H]−."""
        m = self.m
        pos = self.z > 0
        for c in range(m.n):
            if m.el[c] != 'C':
                continue
            for x in m.nb[c]:
                if m.el[x] not in ('N', 'O', 'S') or not m.single(c, x):
                    continue
                Sc, Sx = m.side(c, x), m.side(x, c)
                if Sc is None:
                    continue
                acyl = m.carbonyl_o(c) is not None
                glyco = m.anomeric(c, x) and m.el[x] in ('O', 'N', 'S')
                pendant = len(m.heavy(x)) == 1                  # OH, NH2, SH: their losses are group losses
                x_acyl = any(m.carbonyl_o(j) is not None for j in m.heavy(x) if j != c)
                x_aryl = any(m.aro[j] for j in m.heavy(x) if j != c)
                beta_h = any(m.h(j) for j in m.heavy(c) if j in Sc and j != x) and not m.aro[c]
                if pos and not pendant:                         # the carbon side as the cation, XH leaves
                    if glyco:
                        s, how = 0.55, 'glycosidic cleavage (sugar oxonium ion, B)'
                    elif acyl:
                        s, how = (0.7, 'amide cleavage (acylium, b-type)') if m.el[x] == 'N' else \
                            (0.8, 'ester cleavage (acylium ion)')
                    else:
                        s = 0.9 * m.cation(c, Sc)
                        how = f'C–{m.el[x]} cleavage (carbocation)'
                    if s >= 0.15:
                        self.add(self.make(Sc, {c: (1, 0, -1)}, False, how, s))
                if not beta_h or pendant and not glyco:
                    continue
                if pos:                                         # the heteroatom side takes the proton
                    if glyco:
                        s, how = 0.9, 'glycosidic cleavage (aglycone or sugar loss, Y ion)'
                    elif acyl and m.el[x] == 'N':
                        s, how = 0.65, 'amide cleavage, charge on the amine (y-type)'
                    elif acyl:
                        s, how = 0.35, 'ester cleavage, charge on the alcohol (ketene lost)'
                    elif x_acyl and m.el[x] == 'O':
                        s, how = 0.5, 'ester cleavage: [acid+H]+, alkene lost'
                    elif x_acyl:
                        s, how = 0.35, 'C–N cleavage: [amide+H]+, alkene lost'
                    elif m.el[x] == 'N':
                        s, how = 0.45, 'C–N cleavage: amine, alkene lost'
                    else:
                        s, how = 0.3, f'C–{m.el[x]} cleavage, alkene lost'
                    self.add(self.make(Sx, {x: (1, 0, 1)}, False, how, s))
                else:                                           # the heteroatom side as the anion
                    if glyco:
                        s, how = 0.85, 'glycosidic cleavage (aglycone or sugar loss, Y⁻ ion)'
                    elif x_acyl and m.el[x] == 'O':
                        s, how = 0.7, 'ester cleavage (carboxylate ion)'
                    elif acyl and m.el[x] == 'O':
                        s, how = 0.6, 'ester cleavage (alkoxide ion, ketene lost)'
                    elif x_aryl:
                        s, how = 0.5, 'C–O cleavage (phenoxide ion)'
                    elif acyl:
                        s, how = 0.3, 'amide cleavage (amide anion)'
                    else:
                        s, how = 0.2, f'C–{m.el[x]} cleavage (anion)'
                    self.add(self.make(Sx, {x: (-1, 0, -1)}, False, how, s))

    def charge_remote(self):
        """A structure drawn as an ion: each single bond outside a ring splits, the charged side keeping an H from
        the side that leaves as a neutral molecule."""
        m = self.m
        for i in range(m.n):
            for j in m.heavy(i):
                if j < i or not m.single(i, j):
                    continue
                for a, b in ((i, j), (j, i)):
                    S = m.side(a, b)
                    if S is None or not any(m.q[k] for k in S):
                        continue
                    if not any(m.h(k) for k in m.heavy(b) if k not in S) and not m.h(b):
                        continue
                    self.add(self.make(S, {}, self.odd, 'cleavage, neutral molecule lost (charge-remote)', 0.45))

    def _matches(self):
        if not hasattr(self, '_gm'):
            self._gm = []
            for row in GROUP_LOSSES:
                patt = Chem.MolFromSmarts(row[1])
                self._gm.append((row, self.m.mol.GetSubstructMatches(patt) if patt is not None else ()))
        return self._gm

    def group_losses(self):
        m = self.m
        everything = frozenset(range(m.n))
        for row, matches in self._matches():
            name, _s, rem, att, _lh, radical, scores, what = row
            s = scores.get(self.mode)
            if not s:
                continue
            for match in matches:
                gone = {match[k] for k in rem}
                gone |= {j for g in gone for j in m.nb[g] if m.el[j] == 'H'}
                S = everything - gone
                attach = match[att]
                if attach not in S or m.side(attach, -1, cut=[(attach, g) for g in gone]) != S:
                    continue
                f = self._loss(S, None, attach, radical, name, _lh, s, f'loss of {name} ({what})', None, 1)
                self.add(f)
        for name, smarts, lost, radical, scores, what in FORMULA_LOSSES:
            s = scores.get(self.mode)
            patt = Chem.MolFromSmarts(smarts)
            if not s or patt is None:
                continue
            hits = m.mol.GetSubstructMatches(patt)
            if hits:
                hl = {m.owner[m.order[i]] for i in range(m.n)}
                self.add(self.formula_only(lost, radical, f'loss of {name} ({what})', s, hl))

    def _loss(self, S, base, attach, radical, name, lh, s, how, parent, step):
        """S after a group left from attach: the charge / radical and the H balance placed at attach."""
        m = self.m
        start = parent if parent is not None else None
        odd = (start.odd if start else self.odd) != radical
        lost = {}
        if start is None:
            gone = frozenset(range(m.n)) - S
        else:
            gone = start._idx - S
        for g in gone:
            if m.el[g] != 'H':
                lost[m.key[g]] = lost.get(m.key[g], 0) + 1
        lost['H'] = lh
        target = _sub(start.counts if start else self.pre, lost)
        hc = m.hcap(S)
        if hc is None or any(v < 0 for v in target.values()):
            return None
        have = hc + sum(v[2] for k, v in (base or {}).items() if k in S)
        dh = target.get('H', 0) - have
        if start is None:
            dq = self.z - sum(m.q[i] for i in S)
            rad = 1 if odd else 0
            if dq < 0 and (lh or radical):                  # [M−H]−: an acidic OH elsewhere carries the charge
                acid = self.acid_site(S, (attach,))
                if acid is not None:
                    base = {acid: (-1, 0, -1)}
                    dq = 0
                    dh += 1
            elif dq > 0 and self.mode != 'EI' and (radical or not lh):  # [M+H]+: the proton stays on a basic atom
                site = self.base_site(S, (attach,))
                if site is None:
                    return None
                base = {site: (1, 0, 1)}
                dq = 0
                dh -= 1
        else:                                               # the charge stays where it was
            if self.z > 0 and lh and not radical:           # ... unless the proton leaves with the molecule
                proton = next((k for k, v in base.items() if k in S and v[0] > 0 and v[2] > 0 and m.el[k] in 'NOS'),
                              None)
                if proton is not None:
                    v = base[proton]
                    base = dict(base)
                    base[proton] = (v[0] - 1, v[1], v[2] - 1)
                    dh += 1
            dq = self.z - sum(m.q[i] for i in S) - sum(v[0] for k, v in base.items() if k in S)
            rad = (1 if odd else 0) - sum(v[1] for k, v in base.items() if k in S)
            if rad < 0:                                     # a radical left: the unpaired electron goes with it
                base = {k: list(v) for k, v in base.items()}
                for k, v in base.items():
                    if v[1] > 0 and rad < 0:
                        v[1] -= 1; rad += 1
                base = {k: tuple(v) for k, v in base.items()}
        marks = self.place(S, base, attach, dq, rad, dh)
        if marks is None:
            return None
        return self.make(S, marks, odd, how, s, step, lost_name=name + ('' if not radical or name.endswith('•')
                                                                            else '•'), parent=parent)

    def next_losses(self, f, step):
        m = self.m
        out = []
        S0 = f._idx
        for row, matches in self._matches():
            name, _s, rem, att, _lh, radical, scores, what = row
            s = scores.get(self.mode)
            if not s or (radical and not f.odd):          # even-electron ions lose molecules, not radicals
                continue
            for match in matches:
                gone = {match[k] for k in rem}
                gone |= {j for g in gone for j in m.nb[g] if m.el[j] == 'H'}
                attach = match[att]
                if not gone <= S0 or attach not in S0 or any(g in f._imarks for g in gone):
                    continue
                if any(f._imarks.get(j, (0,))[0] for g in gone | {attach} for j in m.nb[g] + [g]):
                    continue                                # not the group that carries (or sits by) the charge
                S = S0 - gone
                if m.side(attach, -1, cut=[(attach, g) for g in gone] +
                          [(x, y) for x in S0 for y in m.nb[x] if y not in S0]) != S:
                    continue
                g = self._loss(S, f._imarks, attach, radical, name, _lh, f.score * s * 0.85,
                               f'{f.how}, then −{name}', f, step)
                if g is not None:
                    out.append(g)
        for c, v in f._imarks.items():                      # an acylium ion loses CO
            o = m.carbonyl_o(c)
            if v[0] <= 0 or o is None or o not in S0:
                continue
            rest = [j for j in m.heavy(c) if j in S0 and j != o]
            if len(rest) != 1 or m.el[rest[0]] != 'C':     # R–C≡O+ → R+ needs a carbon to take the charge
                continue
            S = S0 - {c, o}
            base = {k: w for k, w in f._imarks.items() if k != c}
            hc = m.hcap(S)
            if hc is None:
                continue
            target = _sub(f.counts, {'C': 1, 'O': 1})
            dh = target.get('H', 0) - hc - sum(w[2] for k, w in base.items() if k in S)
            marks = self.place(S, base, rest[0], v[0], 0, dh)
            g = self.make(S, marks, f.odd, f'{f.how}, then −CO', f.score * 0.75, step, lost_name='CO', parent=f)
            if g is not None:
                out.append(g)
        for x, v in f._imarks.items():                      # a carboxylate ion loses CO2
            c = next((j for j in m.heavy(x) if j in S0), None) if v[0] < 0 and m.el[x] == 'O' else None
            o = m.carbonyl_o(c) if c is not None else None
            if o is None or o not in S0:
                continue
            rest = [j for j in m.heavy(c) if j in S0 and j not in (o, x)]
            if len(rest) != 1:
                continue
            S = S0 - {c, o, x}
            base = {k: w for k, w in f._imarks.items() if k != x}
            hc = m.hcap(S)
            if hc is None:
                continue
            target = _sub(f.counts, {'C': 1, 'O': 2})
            dh = target.get('H', 0) - hc - sum(w[2] for k, w in base.items() if k in S)
            marks = self.place(S, base, rest[0], v[0], 0, dh)
            g = self.make(S, marks, f.odd, f'{f.how}, then −CO2', f.score * 0.8, step, lost_name='CO2', parent=f)
            if g is not None:
                out.append(g)
        if f.z > 0 and not f.odd:                           # an onium ion loses an alkene, its H moving to N / O
            for y in S0:
                if m.el[y] not in ('N', 'O') or m.aro[y]:
                    continue
                for c in m.heavy(y):
                    if c not in S0 or not m.saturated(c) or not m.single(c, y):
                        continue
                    R = m.side(c, y)
                    if R is None or not R <= S0 or any(k in f._imarks for k in R):
                        continue
                    if not any(m.h(j) for j in m.heavy(c) if j in R):
                        continue
                    S = S0 - R
                    hc = m.hcap(S)
                    if hc is None:
                        continue
                    kept = _clean(m.counts(S, hc + sum(w[2] for k, w in f._imarks.items() if k in S)))
                    lost = _clean(_sub(f.counts, kept))
                    g = self.make(S, {k: w for k, w in f._imarks.items() if k in S}, f.odd,
                                  f'{f.how}, then −{chem.formula_text(lost)} (alkene; H to {m.el[y]})',
                                  f.score * 0.6, step, parent=f)
                    if g is not None:
                        out.append(g)
        if f.counts in ({'C': 7, 'H': 7}, {'C': 6, 'H': 5}) and f.z > 0 and not f.odd:
            g = self.formula_only({'C': 2, 'H': 2}, False, f'{f.how}, then −C2H2', f.score * 0.5, f.atoms, f, step)
            if g is not None:
                out.append(g)
        return out

    def merged(self):
        best = {}
        for f in self.out:
            key = (f.formula, f.z, f.odd)
            old = best.get(key)
            if old is None:
                best[key] = f
                continue
            if f.score > old.score:
                f.also = [old.how] + [h for h in old.also if h != f.how]
                best[key] = f
            elif f.how != old.how and f.how not in old.also:
                old.also.append(f.how)
        out = list(best.values())
        for f in out:                                       # agreeing rules make an ion more likely
            f.score = min(1.0, f.score + 0.05 * min(len(f.also), 3))
        out.sort(key=lambda f: -f.mz)
        return out[:500]


def predict(doc, atom_ids=None, mode='ESI+', steps=2):
    """(precursor Fragment, [Fragment] highest m/z first) of the structure(s) with these atoms; (None, []) when
    the structure cannot be read. A structure drawn with a charge is used as it is (mode is ignored)."""
    m = _Mol(doc, atom_ids)
    if m.mol is None or not m.n:
        return None, []
    e = _Engine(m, mode)
    if e.pre.get('H', 0) < 0:
        return None, []
    pre = Fragment(e.mz(e.pre, e.z), e.z, e.odd, e.pre, ion_text(e.pre, e.z, e.odd), '', f'precursor ion {e.label}',
                   1.0, atoms=frozenset(m.owner[m.order[i]] for i in range(m.n)), step=0)
    pre.label = e.label
    pre.drawn_ion = e.charged
    frags = e.run(steps)
    everything = pre.atoms
    for f in frags:
        f.lost_atoms = frozenset(everything - f.atoms) if f.S is not None else frozenset()
    return pre, frags


# ------------------------------------------------------------------------------------------- drawing a fragment
def draw_doc(doc, frag):
    """(Doc holding the fragment ion, '') with the drawing's atom ids and places, its charge, unpaired electron and
    H as worked out; (None, why) when it cannot be drawn."""
    from model import Doc
    if frag.S is None:
        return None, 'This fragment is worked out from the formula only, so it has no structure to draw.'
    full, grown = doc.expanded()
    owners = {}
    for a, more in grown.items():
        for g in more:
            owners[g] = a
    piece = set()
    for fid in frag.S:
        o = owners.get(fid, fid)
        if o not in doc.atoms:
            return None, 'The fragment could not be matched to the drawing.'
        piece.add(o)
    for o in piece:
        group = {o} | grown.get(o, set())
        if not group <= frag.S:
            return None, f'This fragment splits the {doc.atoms[o].label} label: expand it first, then draw it.'
    for fid in frag.marks:
        o = owners.get(fid, fid)
        if fid != o or doc.atoms[o].label:
            return None, f'The charge of this fragment sits inside the {doc.atoms[o].label} label: expand it first.'
    d = Doc()
    d.style = dict(doc.style)
    d.insert_dict(doc.to_dict([('atom', a) for a in piece]), keep_ids=True)
    target = {}
    for a in d.atoms:
        if d.atoms[a].label:
            continue
        cut = sum(doc.bonds[b].order for _n, b in doc.neighbors(a) if _n not in piece)
        target[a] = doc.implicit_h(a) + cut + frag.marks.get(a, (0, 0, 0))[2]
    for a, (dq, rad, _dh) in frag.marks.items():
        d.atoms[a].charge += dq
        d.atoms[a].radical = min(2, d.atoms[a].radical + rad)
    for a in target:
        d.atoms[a].hcount = None

    def natural(a):
        return d.implicit_h(a)
    deficit = {a: natural(a) - t for a, t in target.items() if t < natural(a)}
    for a in sorted(deficit, key=lambda k: -deficit[k]):        # too few H for the bonds: make a double bond
        guard = 0
        while deficit.get(a, 0) > 0 and guard < 6:
            guard += 1
            nbrs = [(n, b) for n, b in d.neighbors(a) if not d.atoms[n].label and d.bonds[b].order < 3]
            pair = next(((n, b) for n, b in nbrs if deficit.get(n, 0) > 0), None)
            if pair is None and deficit[a] >= 2:
                spare = next(((n, b) for n, b in nbrs if target.get(n, 0) >= 1 and n not in deficit), None)
                if spare is not None:                       # one of a's missing H moves to the neighbour
                    target[spare[0]] -= 1
                    target[a] += 1
                    deficit[spare[0]] = 1
                    deficit[a] -= 1
                    pair = spare
            if pair is None:
                break
            n, b = pair
            d.bonds[b].order += 1
            deficit[a] -= 1
            deficit[n] -= 1
    for a in [a for a in deficit if deficit[a] > 0]:      # the two ends of a new diene: a–b=c–e becomes a=b–c=e
        if deficit.get(a, 0) <= 0:
            continue
        for b, ab in d.neighbors(a):
            if deficit.get(b, 0) > 0 or d.bonds[ab].order >= 3:
                continue
            for c, bc in d.neighbors(b):
                if c == a:
                    continue
                if d.bonds[bc].order != 2:
                    _flip_ring(d, bc)
                if d.bonds[bc].order != 2:
                    continue
                e = next(((e, ce) for e, ce in d.neighbors(c) if e not in (a, b) and deficit.get(e, 0) > 0
                          and d.bonds[ce].order < 3), None)
                if e is None:
                    continue
                d.bonds[bc].order = 1
                d.bonds[ab].order += 1
                d.bonds[e[1]].order += 1
                deficit[a] -= 1
                deficit[e[0]] -= 1
                break
            if deficit.get(a, 0) <= 0:
                break
    for a, t in target.items():
        if natural(a) != t:
            d.atoms[a].hcount = max(0, t)
    if any(d.valence_problem(a) for a in d.atoms):
        return None, ('Drawing this fragment needs bonds moved round (a quinoid or rearranged structure); it is '
                      'listed, but not drawn automatically.')
    return d, ''


def _flip_ring(d, bond):
    """Swaps single and double bonds round a six-membered ring of alternating bonds holding bond (the other Kekulé
    structure of a benzene ring), so that bond becomes double."""
    a0, b0 = d.bonds[bond].a, d.bonds[bond].b
    for path in _rings6(d, a0, b0):
        ids = [d.bond_between(path[k], path[(k + 1) % 6]) for k in range(6)]
        orders = [d.bonds[i].order for i in ids]
        if sorted(orders) == [1, 1, 1, 2, 2, 2] and all(orders[k] != orders[(k + 1) % 6] for k in range(6)):
            for i in ids:
                d.bonds[i].order = 3 - d.bonds[i].order
            return


def _rings6(d, a, b):
    """Six-membered rings through bond a–b, as atom lists starting a, b."""
    out = []

    def walk(path):
        if len(path) == 6:
            if d.bond_between(path[-1], a) is not None:
                out.append(list(path))
            return
        for n, _bd in d.neighbors(path[-1]):
            if n not in path:
                walk(path + [n])
    walk([a, b])
    return out
