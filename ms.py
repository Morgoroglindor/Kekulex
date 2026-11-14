"""Mass spectrometry helpers: breaking a bond into fragment ions (as in EI / CID fragmentation schemes) and matching a
measured m/z to the structures drawn or to elemental compositions.

Fragments are worked out by electron bookkeeping on the two pieces left when a bond breaks:
  - a radical cation (an atom carrying both • and +, as in M•+): the bond's electrons go one to each side and pair
    with the radical site, so one piece keeps the charge (an even-electron cation) and the other leaves as a radical -
    "A+ + B•" or "A• + B+" (alpha cleavage of acetone M•+: CH3C≡O+ m/z 43 + CH3•);
  - an even-electron cation (+ without •, as in [M+H]+): inductive cleavage - the charged site leaves neutral and the
    other piece takes the + ("A+ + B");
  - a neutral molecule: homolysis "A• + B•" or heterolysis "A+ + B−" / "A− + B+".
A bond in a ring leaves one piece, so it is not offered.
"""
import chem
from model import Doc


def sides(doc, b):
    """The atoms on each side of bond b ([set, set]) or None when it is in a ring."""
    bd = doc.bonds[b]
    seen, todo = {bd.a}, [bd.a]
    while todo:
        x = todo.pop()
        for n, bb in doc.neighbors(x):
            if bb != b and n not in seen:
                seen.add(n); todo.append(n)
    if bd.b in seen:
        return None
    other, todo = {bd.b}, [bd.b]
    while todo:
        x = todo.pop()
        for n, bb in doc.neighbors(x):
            if bb != b and n not in other:
                other.add(n); todo.append(n)
    return [seen, other]


def options(doc, b):
    """[(spec, description)] for breaking bond b: spec is (mark on bd.a's side, mark on bd.b's side, cleared atoms),
    a mark being '+', '-', '•' or ''. Empty when the bond is in a ring."""
    parts = sides(doc, b)
    if parts is None:
        return []
    atoms = parts[0] | parts[1]
    radical_sites = [a for a in atoms if doc.atoms[a].radical and doc.atoms[a].charge]
    charge = sum(doc.atoms[a].charge for a in atoms)
    unpaired = sum(doc.atoms[a].radical for a in atoms)
    out = []
    if unpaired % 2 and charge:                             # an odd-electron ion such as M•+
        sign = '+' if charge > 0 else '-'
        site = radical_sites[:1]
        out.append(((sign, '•', site), f'charge stays left: A{sign} + B•'))
        out.append((('•', sign, site), f'charge stays right: A• + B{sign}'))
        return out
    if charge:                                              # even-electron ion, e.g. [M+H]+: inductive cleavage
        sign = '+' if charge > 0 else '-'
        for k, (mine, theirs) in enumerate(((parts[0], parts[1]), (parts[1], parts[0]))):
            sites = [a for a in theirs if doc.atoms[a].charge and (doc.atoms[a].charge > 0) == (charge > 0)]
            if sites:
                spec = (sign, '', sites[:1]) if k == 0 else ('', sign, sites[:1])
                out.append((spec, 'inductive: ' + ('A' + sign + ' + B (neutral)' if k == 0 else 'A (neutral) + B' + sign)))
        out.append((('•', '•', []), 'homolysis: A• + B•'))
        return out
    out.append((('•', '•', []), 'homolysis: A• + B•'))
    out.append((('+', '-', []), 'heterolysis: A+ + B−'))
    out.append((('-', '+', []), 'heterolysis: A− + B+'))
    return out


def pieces(doc, b, spec):
    """The two fragments of breaking bond b as new Docs (atom ids kept): [(Doc, cut atom)] for bd.a's and bd.b's side."""
    parts = sides(doc, b)
    bd = doc.bonds[b]
    mark_a, mark_b, cleared = spec
    out = []
    for side, cut, mark in ((parts[0], bd.a, mark_a), (parts[1], bd.b, mark_b)):
        d = Doc()
        d.style = dict(doc.style)
        d.insert_dict(doc.to_dict([('atom', a) for a in side]), keep_ids=True)
        for a in cleared:                                   # the old ion site gives up its charge and radical
            if a in d.atoms:
                d.atoms[a].charge = 0
                d.atoms[a].radical = 0
        at = d.atoms[cut]
        if mark == '•':
            at.radical = min(at.radical + 1, 2)
        elif mark == '+':
            at.charge += 1
        elif mark == '-':
            at.charge -= 1
        out.append((d, cut))
    return out


def describe(doc, atom_ids=None, sep='  '):
    """'C2H3O+  m/z 43.0178' for an ion, 'CH3•  15.0235 (neutral)' for a neutral piece (sep between the two)."""
    info = chem.analyse(doc, atom_ids)
    if not info.atoms:
        return ''
    dots = '•' * (info.radicals % 3)
    if info.charge:
        z = abs(info.charge)
        mz = (info.exact - info.charge * chem.ELECTRON) / z
        return f'{info.formula}{dots}{chem.charge_text(info.charge)}{sep}m/z {mz:.4f}'
    return f'{info.formula}{dots}{sep}{info.exact:.4f} (neutral)'


# ------------------------------------------------------------------------------------------- matching a measured m/z
def ppm(measured, calc):
    return (measured - calc) / calc * 1e6


def drawn_matches(doc, units, mz, tol_ppm):
    """[(where, ion label, formula, calculated m/z, ppm, atom ids)] of the drawn structures' ions (adducts included,
    charged structures as themselves) within tol_ppm of mz. units: [(name, atom ids)]."""
    out = []
    for name, ids in units:
        info = chem.analyse(doc, ids)
        if not info.atoms:
            continue
        for ion in info.ions:
            d = ppm(mz, ion.mz)
            if abs(d) <= tol_ppm:
                out.append((name, ion.label, chem.formula_text(ion.counts), ion.mz, d, ids))
    out.sort(key=lambda r: abs(r[4]))
    return out


DEFAULT_LIMITS = {'C': 80, 'H': 160, 'N': 10, 'O': 25, 'S': 3}
MAX_LIMITS = {'C': 120, 'H': 240, 'N': 20, 'O': 40, 'S': 6, 'P': 4, 'F': 20, 'Cl': 6, 'Br': 4, 'I': 3, 'Si': 6, 'B': 4,
              'Na': 1, 'K': 1}


def limits_from_text(text):
    """'C H N O S Cl' or 'C0-40 N0-4 O ...' -> {element: most}; unknown symbols are skipped."""
    import re
    out = {}
    for tok in text.replace(',', ' ').split():
        m = re.fullmatch(r'([A-Z][a-z]?)(?:(\d+)\s*-\s*)?(\d+)?', tok)
        if not m or m.group(1) not in MAX_LIMITS:
            continue
        out[m.group(1)] = int(m.group(3)) if m.group(3) else DEFAULT_LIMITS.get(m.group(1), MAX_LIMITS[m.group(1)])
    if 'H' not in out and out:
        out['H'] = MAX_LIMITS['H']
    return out


def compositions(mz, z=1, tol_ppm=5.0, limits=None, max_results=50):
    """Elemental compositions of an ion of charge z whose m/z is within tol_ppm of mz:
    [(formula, calculated m/z, ppm, RDBE, 'even' / 'odd' electron)], best first. RDBE (rings plus double bonds) of
    the ion must be at least -0.5; H follows from the other elements."""
    limits = dict(limits or DEFAULT_LIMITS)
    mass = {e: chem.PT.GetMostCommonIsotopeMass(e) for e in limits}
    target = mz * abs(z) + z * chem.ELECTRON                 # the neutral atoms' mass
    tol = target * tol_ppm * 1e-6
    heavy = [e for e in limits if e != 'H']
    heavy.sort(key=lambda e: -mass[e])
    hmax = limits.get('H', 0)
    out = []

    def valence(e):
        return {'C': 4, 'Si': 4, 'N': 3, 'P': 3, 'B': 3, 'O': 2, 'S': 2}.get(e, 1)

    def rec(k, m, counts):
        if m > target + tol:
            return
        if k == len(heavy):
            rest = target - m
            h = round(rest / mass['H']) if 'H' in mass else 0
            if h < 0 or h > hmax:
                return
            total = m + h * mass.get('H', 0)
            if abs(total - target) > tol:
                return
            c = dict(counts)
            if h:
                c['H'] = h
            if not c:
                return
            # RDBE = 1 + sum(n (v - 2)) / 2 over all atoms (H and halogens: v = 1)
            rdbe = 1 + sum(n * (valence(e) - 2) for e, n in c.items()) / 2
            if rdbe < -0.5 or (rdbe > 0.5 * sum(c.values()) + 1):
                return
            nc = c.get('C', 0)
            if nc and not (0.2 <= h / nc <= 3.1 and c.get('N', 0) / nc <= 1.3 and c.get('O', 0) / nc <= 1.2
                           and c.get('S', 0) / nc <= 0.8 and c.get('P', 0) / nc <= 0.3):
                return                                      # outside the usual ratios of real compounds
            if c.get('Na') and c.get('K'):
                return
            electrons = 'odd' if abs(rdbe - round(rdbe)) < 1e-9 else 'even'
            calc = (total - z * chem.ELECTRON) / abs(z)
            out.append((chem.formula_text(c), calc, ppm(mz, calc), rdbe, electrons))
            return
        e = heavy[k]
        for n in range(0, limits[e] + 1):
            mm = m + n * mass[e]
            if mm > target + tol:
                break
            if n:
                counts[e] = n
            else:
                counts.pop(e, None)
            rec(k + 1, mm, counts)
        counts.pop(e, None)
    rec(0, 0.0, {})
    out.sort(key=lambda r: (abs(r[2]), -r[3]))
    return out[:max_results]


def rdbe_of(counts):
    v = {'C': 4, 'Si': 4, 'N': 3, 'P': 3, 'B': 3, 'O': 2, 'S': 2}
    return 1 + sum(n * (v.get(e, 1) - 2) for e, n in counts.items()) / 2

