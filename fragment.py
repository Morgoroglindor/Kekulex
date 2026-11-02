"""
fragment.py - ChemDraw's fragmentation tools: drag a line across bonds of a structure.

    crossed_bonds(doc, x1, y1, x2, y2)     # bond ids the line crosses
    mass_fragments(doc, x1, y1, ...)       # Mass Fragmentation: formula and exact mass of each side, as if the
                                           # bonds broke homolytically (each atom keeps its hydrogens: radicals)
    disconnect(doc, line, retro=True)      # Retrosynthesis (product -> ⇒ -> pieces) / Synthesis (pieces -> product):
                                           # the pieces are copies with the cut bonds removed (filled up with H)

Not a prediction: the user picks the bonds.
"""
import math

import chem
import model


def _cross(ax, ay, bx, by, cx, cy, dx, dy):
    """Whether segments a-b and c-d cross."""
    def orient(px, py, qx, qy, rx, ry):
        return (qx - px) * (ry - py) - (qy - py) * (rx - px)
    d1, d2 = orient(cx, cy, dx, dy, ax, ay), orient(cx, cy, dx, dy, bx, by)
    d3, d4 = orient(ax, ay, bx, by, cx, cy), orient(ax, ay, bx, by, dx, dy)
    return d1 * d2 < 0 and d3 * d4 < 0


def crossed_bonds(doc, x1, y1, x2, y2):
    out = []
    for i, b in doc.bonds.items():
        pa, pb = doc.atoms[b.a], doc.atoms[b.b]
        if _cross(pa.x, pa.y, pb.x, pb.y, x1, y1, x2, y2):
            out.append(i)
    return out


def _side(x1, y1, x2, y2, px, py):
    return (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)


def pieces(doc, bonds, line):
    """The molecules the cut bonds belong to, split without those bonds: ([atom ids of side A], [of side B]) - each
    side may hold several pieces; a ring bond alone gives everything on one side."""
    cut = set(bonds)
    touched = set()
    for i in cut:
        touched |= {doc.bonds[i].a, doc.bonds[i].b}
    mols = [set(f) for f in doc.fragments() if set(f) & touched]
    atoms = set().union(*mols) if mols else set()
    adj = {a: set() for a in atoms}
    for i, b in doc.bonds.items():
        if i not in cut and b.a in atoms and b.b in atoms:
            adj[b.a].add(b.b); adj[b.b].add(b.a)
    comps, seen = [], set()
    for a in atoms:
        if a in seen:
            continue
        stack, comp = [a], set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u); stack.extend(adj[u] - comp)
        seen |= comp
        comps.append(comp)
    x1, y1, x2, y2 = line
    sides = ([], [])
    for comp in comps:
        cx = sum(doc.atoms[a].x for a in comp) / len(comp)
        cy = sum(doc.atoms[a].y for a in comp) / len(comp)
        sides[0 if _side(x1, y1, x2, y2, cx, cy) > 0 else 1].append(comp)
    return comps, sides


def _counts(doc, atoms):
    """{element: n} of these atoms with the hydrogens they have in the intact structure (nicknames in full)."""
    full, grown = doc.expanded()
    ids = chem.expand_ids(grown, atoms)
    counts = {}
    for a in ids:
        at = full.atoms[a]
        counts[at.el] = counts.get(at.el, 0) + 1
        h = full.implicit_h(a)
        if h:
            counts['H'] = counts.get('H', 0) + h
    return counts


def mass_fragments(doc, line):
    """[(formula, exact mass, (x, y) where to put the label)] for the sides of the cut (one entry when only a ring
    bond is crossed); [] when the line crosses no bond."""
    bonds = crossed_bonds(doc, *line)
    if not bonds:
        return []
    comps, sides = pieces(doc, bonds, line)
    groups = [s for s in sides if s] if len(comps) > 1 else [[set().union(*comps)]]
    bottom = max(doc.atoms[a].y for a in set().union(*comps))
    out = []
    for side in groups:                                     # each label under its own fragment
        atoms = set().union(*side)
        counts = _counts(doc, atoms)
        cx, _cy = _centroid(doc, atoms)
        out.append((chem.formula_text(counts), chem._mono(counts), (cx, bottom + 0.75)))
    if len(out) == 2 and abs(out[0][2][0] - out[1][2][0]) < 2.2:   # side by side under the structure
        (fa, ma, (xa, ya)), (fb, mb, (xb, yb)) = out
        mid = (xa + xb) / 2
        left, right = (0, 1) if xa <= xb else (1, 0)
        out[left] = (out[left][0], out[left][1], (mid - 1.2, ya))
        out[right] = (out[right][0], out[right][1], (mid + 1.2, yb))
    return out


def _centroid(doc, atoms):
    return (sum(doc.atoms[a].x for a in atoms) / len(atoms), sum(doc.atoms[a].y for a in atoms) / len(atoms))


def disconnect(doc, line, retro=True):
    """Draw the reaction for a disconnection: retro - the structure, a retrosynthetic arrow, then its pieces with
    '+' between them; synthesis - the pieces with '+', a reaction arrow, then the structure. Returns the new items
    ([] when the line crosses no bond)."""
    bonds = crossed_bonds(doc, *line)
    if not bonds:
        return []
    comps, _sides = pieces(doc, bonds, line)
    product = set().union(*comps)
    cut_pairs = {frozenset((doc.bonds[i].a, doc.bonds[i].b)) for i in bonds}
    xs = [doc.atoms[a].x for a in product]; ys = [doc.atoms[a].y for a in product]
    cy = (min(ys) + max(ys)) / 2
    gap = 1.0
    new = []
    sub_dicts = []
    for comp in comps:
        d = doc.to_dict([('atom', a) for a in comp])
        d['bonds'] = [b for b in d['bonds'] if frozenset((b['a'], b['b'])) not in cut_pairs]
        d['fills'] = []
        cxs = [a['x'] for a in d['atoms']]; cys = [a['y'] for a in d['atoms']]
        sub_dicts.append((d, min(cxs), max(cxs), (min(cys) + max(cys)) / 2))
    widths = [mx - mn for _d, mn, mx, _c in sub_dicts]
    plus_w = 1.0
    if retro:
        x = max(xs) + gap
        new.append(('arrow', doc.add_arrow(x, cy, x + 2.5, cy, 'retro')))
        x += 2.5 + gap
        for k, (d, mn, mx, mid) in enumerate(sub_dicts):
            if k:
                new.append(('text', doc.add_text(x, cy - 0.35, '+', align='center')))
                x += plus_w
            new += doc.insert_dict(d, x - mn, cy - mid)
            x += widths[k] + gap
    else:
        total = sum(widths) + gap * len(widths) + plus_w * (len(widths) - 1) + 2.5 + gap
        x = min(xs) - total
        for k, (d, mn, mx, mid) in enumerate(sub_dicts):
            if k:
                new.append(('text', doc.add_text(x, cy - 0.35, '+', align='center')))
                x += plus_w
            new += doc.insert_dict(d, x - mn, cy - mid)
            x += widths[k] + gap
        new.append(('arrow', doc.add_arrow(x, cy, x + 2.5, cy, 'forward')))
    return new
