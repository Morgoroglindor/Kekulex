"""
model.py - the drawing: atoms and bonds (any number of separate structures), reaction arrows and text, with the
geometry the editor needs (where the next bond goes, rings, chains, the side of a double bond, implicit hydrogens) and
the MOL conversion RDKit reads.

    doc = Doc()
    a = doc.add_atom(0, 0); b = doc.add_atom(1, 0, 'O'); doc.add_bond(a, b, 2)
    doc.to_molblock()            # V2000 MOL text (y up, 1.5 A bonds)
    Doc.from_dict(doc.to_dict()) # the .sed file is this dict as JSON

Coordinates are in bond lengths (one bond = 1.0) with y pointing down, as on the screen. Every item (atom, bond, arrow,
text) has an integer id from one counter, so a selection can hold any of them as (kind, id).
"""
import copy
import math
from dataclasses import asdict, dataclass

MERGE = 0.3                                                 # a new atom this close to an existing one is that atom
MOL_BOND = 1.5                                              # bond length written into MOL files (A, RDKit's own)

# Usual valences (the lowest that fits the bonds is used) and what a charge does to them
VALENCE = {'H': (1,), 'B': (3,), 'C': (4,), 'N': (3, 5), 'O': (2,), 'F': (1,), 'Si': (4,), 'P': (3, 5), 'S': (2, 4, 6),
           'Cl': (1, 3, 5, 7), 'Se': (2, 4, 6), 'Br': (1, 3, 5), 'I': (1, 3, 5, 7), 'As': (3, 5), 'Te': (2, 4, 6)}
ACCEPTS_CHARGE_AS_BOND = {'N', 'O', 'P', 'S', 'Se', 'As', 'F', 'Cl', 'Br', 'I'}   # N+ / O+ gain a bond, C+ / B- lose one


@dataclass
class Atom:
    x: float
    y: float
    el: str = 'C'
    charge: int = 0


@dataclass
class Bond:
    a: int
    b: int
    order: int = 1                                          # 1, 2, 3
    stereo: str = ''                                        # '' | 'wedge' | 'hash' (the narrow end is atom a)
    side: int | None = None                                 # double bond: None = automatic, 0 = centred, +1 / -1 = side


@dataclass
class Arrow:
    x1: float
    y1: float
    x2: float
    y2: float
    kind: str = 'forward'                                   # 'forward' | 'equilibrium'


@dataclass
class Text:
    x: float                                                # top-left corner
    y: float
    text: str


KINDS = {'atom': Atom, 'bond': Bond, 'arrow': Arrow, 'text': Text}


def angle_of(dx, dy):
    return math.atan2(dy, dx)


def snap_angle(a, step_deg=15):
    s = math.radians(step_deg)
    return round(a / s) * s


class Doc:
    def __init__(self):
        self.atoms = {}
        self.bonds = {}
        self.arrows = {}
        self.texts = {}
        self.style = {}                                     # overrides of render.DEFAULT_STYLE
        self._next = 1

    # ------------------------------------------------------------------------------------------- items
    def _id(self):
        self._next += 1
        return self._next - 1

    def table(self, kind):
        return {'atom': self.atoms, 'bond': self.bonds, 'arrow': self.arrows, 'text': self.texts}[kind]

    def add_atom(self, x, y, el='C', charge=0):
        i = self._id()
        self.atoms[i] = Atom(x, y, el, charge)
        return i

    def atom_at_or_new(self, x, y, el='C', exclude=()):
        """The existing atom within MERGE of (x, y), else a new one there."""
        hit = self.find_atom(x, y, MERGE, exclude)
        return hit if hit is not None else self.add_atom(x, y, el)

    def add_bond(self, a, b, order=1, stereo=''):
        """Bond a-b (the existing one if the two are already bonded; None if a == b)."""
        if a == b:
            return None
        old = self.bond_between(a, b)
        if old is not None:
            return old
        i = self._id()
        self.bonds[i] = Bond(a, b, order, stereo)
        return i

    def add_arrow(self, x1, y1, x2, y2, kind='forward'):
        i = self._id()
        self.arrows[i] = Arrow(x1, y1, x2, y2, kind)
        return i

    def add_text(self, x, y, text):
        i = self._id()
        self.texts[i] = Text(x, y, text)
        return i

    def bond_between(self, a, b):
        for i, bd in self.bonds.items():
            if (bd.a == a and bd.b == b) or (bd.a == b and bd.b == a):
                return i
        return None

    def neighbors(self, a):
        """[(neighbour atom, bond id)]"""
        out = []
        for i, bd in self.bonds.items():
            if bd.a == a:
                out.append((bd.b, i))
            elif bd.b == a:
                out.append((bd.a, i))
        return out

    def remove(self, kind, i):
        if kind == 'atom':
            self.remove_atom(i)
        elif kind == 'bond':
            self.remove_bond(i)
        else:
            self.table(kind).pop(i, None)

    def remove_atom(self, a):
        for _n, b in self.neighbors(a):
            self.remove_bond(b, prune=False)
        self.atoms.pop(a, None)

    def remove_bond(self, i, prune=True):
        """Remove a bond; a carbon left without bonds goes with it (it would show up as a lone CH4)."""
        bd = self.bonds.pop(i, None)
        if bd and prune:
            for a in (bd.a, bd.b):
                at = self.atoms.get(a)
                if at and at.el == 'C' and not at.charge and not self.neighbors(a):
                    del self.atoms[a]

    def merge_atom(self, src, dst):
        """Atom src is dropped onto dst: its bonds move to dst (a bond to dst itself is dropped) and src goes."""
        for n, b in self.neighbors(src):
            bd = self.bonds[b]
            if n == dst or self.bond_between(n, dst) is not None:
                del self.bonds[b]
                continue
            if bd.a == src:
                bd.a = dst
            else:
                bd.b = dst
        del self.atoms[src]

    def is_empty(self):
        return not (self.atoms or self.arrows or self.texts)

    # ------------------------------------------------------------------------------------------- hit testing
    def find_atom(self, x, y, r, exclude=()):
        best, bd = None, r
        for i, a in self.atoms.items():
            if i in exclude:
                continue
            d = math.hypot(a.x - x, a.y - y)
            if d <= bd:
                best, bd = i, d
        return best

    def find_bond(self, x, y, r):
        best, bd = None, r
        for i, b in self.bonds.items():
            p, q = self.atoms[b.a], self.atoms[b.b]
            d = seg_dist(x, y, p.x, p.y, q.x, q.y)
            if d <= bd:
                best, bd = i, d
        return best

    def find_arrow(self, x, y, r):
        best, bd = None, r
        for i, a in self.arrows.items():
            d = seg_dist(x, y, a.x1, a.y1, a.x2, a.y2)
            if d <= bd:
                best, bd = i, d
        return best

    # ------------------------------------------------------------------------------------------- chemistry
    def bond_sum(self, a):
        return sum(self.bonds[b].order for _n, b in self.neighbors(a))

    def _valences(self, a):
        at = self.atoms[a]
        base = VALENCE.get(at.el)
        if base is None:
            return None
        if at.charge and at.el in ACCEPTS_CHARGE_AS_BOND:
            return tuple(v + at.charge for v in base if v + at.charge >= 0)
        return tuple(v - abs(at.charge) for v in base if v - abs(at.charge) >= 0)

    def implicit_h(self, a):
        """Hydrogens the atom carries without them being drawn (0 for elements without a usual valence)."""
        vals = self._valences(a)
        if not vals:
            return 0
        s = self.bond_sum(a)
        for v in vals:
            if v >= s:
                return v - s
        return 0

    def valence_error(self, a):
        vals = self._valences(a)
        return bool(vals) and self.bond_sum(a) > max(vals)

    def fragments(self, atom_ids=None):
        """Connected structures, as sorted lists of atom ids (only atoms in atom_ids, when given)."""
        todo = set(self.atoms if atom_ids is None else (a for a in atom_ids if a in self.atoms))
        out = []
        while todo:
            seed = min(todo)
            comp, stack = {seed}, [seed]
            while stack:
                for n, _b in self.neighbors(stack.pop()):
                    if n in todo and n not in comp:
                        comp.add(n); stack.append(n)
            todo -= comp
            out.append(sorted(comp))
        out.sort(key=lambda c: (self.atoms[c[0]].x, self.atoms[c[0]].y))
        return out

    # ------------------------------------------------------------------------------------------- geometry
    def bond_angles(self, a, skip=None):
        at = self.atoms[a]
        return [angle_of(self.atoms[n].x - at.x, self.atoms[n].y - at.y) for n, b in self.neighbors(a) if b != skip]

    def next_bond_angle(self, a):
        """Direction of a new bond from atom a: up-right for a lone atom, zigzag (away from the crowd) after one bond,
        straight on after a triple bond or two double bonds, else the middle of the widest gap between the bonds."""
        nb = self.neighbors(a)
        at = self.atoms[a]
        if not nb:
            return math.radians(-30)
        angles = self.bond_angles(a)
        if len(nb) == 1:
            base = angles[0]
            if any(self.bonds[b].order == 3 for _n, b in nb):
                return base + math.pi
            cands = [base + math.radians(120), base - math.radians(120)]

            def room(ang):
                x, y = at.x + math.cos(ang), at.y + math.sin(ang)
                others = [math.hypot(o.x - x, o.y - y) for i, o in self.atoms.items() if i != a and i != nb[0][0]]
                return min(others, default=9.0)
            r = [room(c) for c in cands]
            if abs(r[0] - r[1]) > 1e-6:
                return cands[0] if r[0] > r[1] else cands[1]
            return min(cands, key=lambda c: abs(math.sin(c)))   # a free choice: the more horizontal one (a zigzag)
        angles.sort()
        gaps = [((angles[(k + 1) % len(angles)] - angles[k]) % (2 * math.pi) or 2 * math.pi, angles[k])
                for k in range(len(angles))]
        gap, start = max(gaps)
        return start + gap / 2

    def free_direction(self, a):
        """Unit vector pointing away from the bonds of atom a (where its hydrogens and labels can go)."""
        at = self.atoms[a]
        sx = sy = 0.0
        for n, _b in self.neighbors(a):
            o = self.atoms[n]
            d = math.hypot(o.x - at.x, o.y - at.y) or 1.0
            sx += (o.x - at.x) / d; sy += (o.y - at.y) / d
        d = math.hypot(sx, sy)
        if d < 1e-6:
            return (1.0, 0.0) if not self.neighbors(a) else (0.0, -1.0)
        return (-sx / d, -sy / d)

    def ring_of_bond(self, i, limit=8):
        """Atoms of the smallest ring through bond i (in order around the ring), or None."""
        bd = self.bonds[i]
        prev = {bd.a: None}
        frontier = [bd.a]
        for _depth in range(limit - 1):
            nxt = []
            for u in frontier:
                for n, b in self.neighbors(u):
                    if b == i or n in prev:
                        continue
                    prev[n] = u
                    if n == bd.b:
                        path = [n]
                        while prev[path[-1]] is not None:
                            path.append(prev[path[-1]])
                        return path
                    nxt.append(n)
            frontier = nxt
        return None

    def double_side(self, i):
        """Where a double bond's second line goes: +1 / -1 = that side of a->b (left normal = +1), 0 = centred."""
        bd = self.bonds[i]
        if bd.side is not None:
            return bd.side
        p, q = self.atoms[bd.a], self.atoms[bd.b]
        nx, ny = -(q.y - p.y), q.x - p.x                   # normal of a->b

        def side_of(x, y):
            v = (x - p.x) * nx + (y - p.y) * ny
            return 0 if abs(v) < 1e-6 else (1 if v > 0 else -1)
        ring = self.ring_of_bond(i)
        if ring:
            cx = sum(self.atoms[k].x for k in ring) / len(ring)
            cy = sum(self.atoms[k].y for k in ring) / len(ring)
            s = side_of(cx, cy)
            if s:
                return s
        na = [n for n, _b in self.neighbors(bd.a) if n != bd.b]
        nb = [n for n, _b in self.neighbors(bd.b) if n != bd.a]
        if not na and not nb:
            return 0
        if (not na and len(nb) != 1) or (not nb and len(na) != 1):
            return 0                                        # C=O of a carbonyl, =CH2 on a branch point: centred
        total = sum(side_of(self.atoms[n].x, self.atoms[n].y) for n in na + nb)
        if total == 0:
            return 0
        return 1 if total > 0 else -1

    def ring_points(self, n, center, start_angle):
        R = 1 / (2 * math.sin(math.pi / n))
        cx, cy = center
        return [(cx + R * math.cos(start_angle + 2 * math.pi * k / n), cy + R * math.sin(start_angle + 2 * math.pi * k / n))
                for k in range(n)]

    def ring_on_bond(self, i, n):
        """Points of an n-ring fused on bond i, on the side away from the bond's other neighbours."""
        bd = self.bonds[i]
        p, q = self.atoms[bd.a], self.atoms[bd.b]
        mx, my = (p.x + q.x) / 2, (p.y + q.y) / 2
        L = math.hypot(q.x - p.x, q.y - p.y) or 1.0
        nx, ny = -(q.y - p.y) / L, (q.x - p.x) / L
        score = 0.0
        for a in (bd.a, bd.b):
            for o, _b in self.neighbors(a):
                if o in (bd.a, bd.b):
                    continue
                score += (self.atoms[o].x - mx) * nx + (self.atoms[o].y - my) * ny
        s = -1 if score > 0 else 1
        apothem = 1 / (2 * math.tan(math.pi / n))
        cx, cy = mx + s * nx * apothem, my + s * ny * apothem
        return self.ring_points(n, (cx, cy), angle_of(p.x - cx, p.y - cy))

    def ring_on_atom(self, a, n):
        at = self.atoms[a]
        ang = self.next_bond_angle(a) if self.neighbors(a) else math.radians(90)
        R = 1 / (2 * math.sin(math.pi / n))
        cx, cy = at.x + R * math.cos(ang), at.y + R * math.sin(ang)
        return self.ring_points(n, (cx, cy), ang + math.pi)

    def add_ring(self, pts, aromatic=False):
        """Atoms and bonds of a ring through pts (existing atoms are reused); aromatic: alternate double bonds where
        they fit with the double bonds already there. Returns the ring's atom ids."""
        ids = [self.atom_at_or_new(x, y) for x, y in pts]
        n = len(ids)
        edges = [(ids[k], ids[(k + 1) % n]) for k in range(n)]
        bonds = [self.add_bond(u, v) for u, v in edges]
        if aromatic and n % 2 == 0:
            def has_double(atom, skip):
                return any(self.bonds[b].order >= 2 for _o, b in self.neighbors(atom) if b not in skip)

            def conflicts(par):
                c = 0
                for k in range(par, n, 2):
                    u, v = edges[k]
                    if self.bonds[bonds[k]].order == 1 and (has_double(u, bonds) or has_double(v, bonds)):
                        c += 1
                for k in range(1 - par, n, 2):
                    if self.bonds[bonds[k]].order >= 2:
                        c += 1
                return c
            par = min((0, 1), key=conflicts)
            for k in range(par, n, 2):
                u, v = edges[k]
                b = self.bonds[bonds[k]]
                if b.order == 1 and not has_double(u, bonds) and not has_double(v, bonds):
                    b.order = 2
        return ids

    def chain_points(self, x0, y0, x1, y1):
        """Zigzag chain from (x0, y0) toward (x1, y1) (direction snapped to 30 degrees): the points after the first."""
        ang = snap_angle(angle_of(x1 - x0, y1 - y0), 30)
        dist = math.hypot(x1 - x0, y1 - y0)
        step = math.sqrt(3) / 2
        n = max(1, round(dist / step))
        ux, uy = math.cos(ang), math.sin(ang)
        vx, vy = -uy, ux
        sgn = -1 if vy > 0 else 1                           # the first zigzag step goes up on the screen
        if abs(uy) > 0.99:                                  # (to the right for a vertical chain)
            sgn = 1 if vx > 0 else -1
        return [(x0 + k * step * ux + (k % 2) * 0.5 * vx * sgn, y0 + k * step * uy + (k % 2) * 0.5 * vy * sgn)
                for k in range(1, n + 1)]

    # ------------------------------------------------------------------------------------------- transforms
    def atoms_of(self, sel):
        """Atom ids moved along with a selection [(kind, id)] (the atoms of selected bonds too)."""
        out = set()
        for kind, i in sel:
            if kind == 'atom':
                out.add(i)
            elif kind == 'bond' and i in self.bonds:
                out.update((self.bonds[i].a, self.bonds[i].b))
        return out

    def points_of(self, sel):
        pts = [(self.atoms[a].x, self.atoms[a].y) for a in self.atoms_of(sel) if a in self.atoms]
        for kind, i in sel:
            if kind == 'arrow' and i in self.arrows:
                ar = self.arrows[i]; pts += [(ar.x1, ar.y1), (ar.x2, ar.y2)]
            elif kind == 'text' and i in self.texts:
                pts.append((self.texts[i].x, self.texts[i].y))
        return pts

    def transform(self, sel, fn):
        """Apply fn(x, y) -> (x, y) to everything in sel."""
        for a in self.atoms_of(sel):
            at = self.atoms[a]; at.x, at.y = fn(at.x, at.y)
        for kind, i in sel:
            if kind == 'arrow' and i in self.arrows:
                ar = self.arrows[i]
                ar.x1, ar.y1 = fn(ar.x1, ar.y1); ar.x2, ar.y2 = fn(ar.x2, ar.y2)
            elif kind == 'text' and i in self.texts:
                t = self.texts[i]; t.x, t.y = fn(t.x, t.y)

    def all_items(self):
        return ([('atom', i) for i in self.atoms] + [('bond', i) for i in self.bonds] +
                [('arrow', i) for i in self.arrows] + [('text', i) for i in self.texts])

    # ------------------------------------------------------------------------------------------- files
    def to_dict(self, sel=None):
        """The drawing (or the items in sel, with the bonds between selected atoms) as plain data."""
        if sel is None:
            keep = None
        else:
            keep = set(sel)
            atoms = self.atoms_of(sel)
            keep |= {('atom', a) for a in atoms}
            keep |= {('bond', i) for i, b in self.bonds.items() if b.a in atoms and b.b in atoms}

        def items(kind):
            return [dict(id=i, **asdict(v)) for i, v in self.table(kind).items() if keep is None or (kind, i) in keep]
        return dict(format='structure-editor', version=1, style=dict(self.style), atoms=items('atom'),
                    bonds=items('bond'), arrows=items('arrow'), texts=items('text'))

    @classmethod
    def from_dict(cls, d):
        doc = cls()
        doc.insert_dict(d, keep_ids=True)
        doc.style = dict(d.get('style') or {})
        return doc

    def insert_dict(self, d, dx=0.0, dy=0.0, keep_ids=False):
        """Add the items of a to_dict() result (shifted by dx, dy); returns the new items as [(kind, id)]."""
        new, idmap = [], {}
        for kind in ('atom', 'bond', 'arrow', 'text'):
            for item in d.get(kind + 's', []):
                item = dict(item)
                old = item.pop('id')
                if kind == 'bond':
                    if item['a'] not in idmap or item['b'] not in idmap:
                        continue
                    item['a'], item['b'] = idmap[item['a']], idmap[item['b']]
                for k in ('x', 'x1', 'x2'):
                    if k in item:
                        item[k] += dx
                for k in ('y', 'y1', 'y2'):
                    if k in item:
                        item[k] += dy
                if keep_ids:
                    i = old
                    self._next = max(self._next, i + 1)
                else:
                    i = self._id()
                if kind == 'atom':
                    idmap[old] = i
                self.table(kind)[i] = KINDS[kind](**item)
                new.append((kind, i))
        return new

    def copy(self):
        return copy.deepcopy(self)

    def to_molblock(self, atom_ids=None, title=''):
        """V2000 MOL block of the structures (or of the atoms in atom_ids) with wedges and charges; y up, 1.5 A bonds.
        Returns (text, [atom id in file order])."""
        order = [a for a in self.atoms if atom_ids is None or a in atom_ids]
        index = {a: k + 1 for k, a in enumerate(order)}
        bonds = [b for b in self.bonds.values() if b.a in index and b.b in index]
        lines = [title, '  StructureEditor 2D', '',
                 f'{len(order):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000']
        for a in order:
            at = self.atoms[a]
            lines.append(f'{at.x * MOL_BOND:10.4f}{-at.y * MOL_BOND:10.4f}{0:10.4f} {at.el:<3} 0  0  0  0  0  0  0  0  0  0  0  0')
        for b in bonds:
            st = {'wedge': 1, 'hash': 6}.get(b.stereo, 0)
            lines.append(f'{index[b.a]:3d}{index[b.b]:3d}{b.order:3d}{st:3d}')
        charged = [(index[a], self.atoms[a].charge) for a in order if self.atoms[a].charge]
        for k in range(0, len(charged), 8):
            chunk = charged[k:k + 8]
            lines.append(f'M  CHG{len(chunk):3d}' + ''.join(f'{i:4d}{c:4d}' for i, c in chunk))
        lines.append('M  END')
        return '\n'.join(lines) + '\n', order

    def insert_rdkit(self, mol, center=(0.0, 0.0)):
        """Add an RDKit molecule (laid out in 2D when it has no coordinates), centred on center; returns the new items."""
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        mol = Chem.Mol(mol)
        if mol.GetNumConformers() == 0 or mol.GetConformer().Is3D():
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(mol)
        try:
            Chem.ReapplyMolBlockWedging(mol)
        except (AttributeError, RuntimeError):
            pass
        if not any(b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) for b in mol.GetBonds()):
            Chem.WedgeMolBonds(mol, mol.GetConformer())
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            pass
        conf = mol.GetConformer()
        pos = [conf.GetAtomPosition(k) for k in range(mol.GetNumAtoms())]
        lens = sorted(math.hypot(pos[b.GetBeginAtomIdx()].x - pos[b.GetEndAtomIdx()].x,
                                 pos[b.GetBeginAtomIdx()].y - pos[b.GetEndAtomIdx()].y) for b in mol.GetBonds())
        L = lens[len(lens) // 2] if lens else MOL_BOND
        L = L or MOL_BOND
        cx = sum(p.x for p in pos) / max(len(pos), 1)
        cy = sum(p.y for p in pos) / max(len(pos), 1)
        ids, new = [], []
        for k, atom in enumerate(mol.GetAtoms()):
            i = self.add_atom(center[0] + (pos[k].x - cx) / L, center[1] - (pos[k].y - cy) / L, atom.GetSymbol(),
                              atom.GetFormalCharge())
            ids.append(i); new.append(('atom', i))
        for b in mol.GetBonds():
            bt = b.GetBondType()
            order = 3 if bt == Chem.BondType.TRIPLE else 2 if bt == Chem.BondType.DOUBLE else 1
            stereo = {Chem.BondDir.BEGINWEDGE: 'wedge', Chem.BondDir.BEGINDASH: 'hash'}.get(b.GetBondDir(), '')
            i = self._id()
            self.bonds[i] = Bond(ids[b.GetBeginAtomIdx()], ids[b.GetEndAtomIdx()], order, stereo)
            new.append(('bond', i))
        return new


def seg_dist(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / L2))
    return math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))
