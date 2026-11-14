"""
chem.py - what the drawn structure is: formula, masses, elemental analysis, the m/z and isotope pattern of the usual
LC-MS adducts, SMILES, InChIKey and R/S, E/Z labels (RDKit).

    info = chem.analyse(doc, atom_ids)      # Info (atom_ids None: everything drawn); nicknames count in full
    chem.isotope_pattern(counts, z)         # [(m/z, % of the highest)] for a formula {element: n}
    chem.mol_from_text(text)                # RDKit molecule from SMILES / MOL text (None if neither)

Masses are monoisotopic (most abundant isotope of each element) and average (standard atomic weights). The m/z of an
ion is (neutral atoms + adduct atoms - z electrons) / |z|, so a structure drawn with a charge (a quaternary ammonium,
a zwitterion drawn with its charges) gives the right ions as well.
"""
from dataclasses import dataclass, field

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')                              # RDKit's own warnings would go to the console

ELECTRON = 0.000548579909
PT = Chem.GetPeriodicTable()
H = PT.GetMostCommonIsotopeMass('H')


def _mono(formula):
    """Monoisotopic mass of a small neutral formula given as {element: count}."""
    return sum(key_mass(el) * n for el, n in formula.items())


# Formula keys: an element ('C'), or a labelled isotope: '13C', '15N', and 'D' / 'T' for hydrogen's
def iso_key(el, iso):
    if not iso:
        return el
    return {2: 'D', 3: 'T'}.get(iso, f'{iso}{el}') if el == 'H' else f'{iso}{el}'


def split_key(key):
    """(element, mass number or 0) of a formula key: 'C' -> ('C', 0), '13C' -> ('C', 13), 'D' -> ('H', 2)."""
    if key in ('D', 'T'):
        return 'H', 2 if key == 'D' else 3
    k = 0
    while k < len(key) and key[k].isdigit():
        k += 1
    return (key[k:], int(key[:k])) if k else (key, 0)


def key_mass(key, average=False):
    """Monoisotopic (or average) mass of one atom of a formula key; a labelled isotope has its own exact mass."""
    el, iso = split_key(key)
    if iso:
        return PT.GetMassForIsotope(el, iso)
    return PT.GetAtomicWeight(el) if average else PT.GetMostCommonIsotopeMass(el)


_SUP = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')


def key_text(key, html=False):
    """A formula key as written: C, D, ¹³C (or <sup>13</sup>C)."""
    el, iso = split_key(key)
    if not iso or key in ('D', 'T'):
        return key
    return f'<sup>{iso}</sup>{el}' if html else str(iso).translate(_SUP) + el


# (label, molecules, atoms gained, atoms lost, charge gained). The label's 'M' takes the count.
ADDUCTS_POS = (('M+H', 1, {'H': 1}, {}, 1), ('M+NH4', 1, {'N': 1, 'H': 4}, {}, 1), ('M+Na', 1, {'Na': 1}, {}, 1),
               ('M+K', 1, {'K': 1}, {}, 1), ('M+H-H2O', 1, {'H': 1}, {'H': 2, 'O': 1}, 1), ('M+2H', 1, {'H': 2}, {}, 2),
               ('2M+H', 2, {'H': 1}, {}, 1), ('2M+Na', 2, {'Na': 1}, {}, 1))
ADDUCTS_NEG = (('M-H', 1, {}, {'H': 1}, -1), ('M+Cl', 1, {'Cl': 1}, {}, -1), ('M+HCOO', 1, {'C': 1, 'H': 1, 'O': 2}, {}, -1),
               ('M+CH3COO', 1, {'C': 2, 'H': 3, 'O': 2}, {}, -1), ('M-H-H2O', 1, {}, {'H': 3, 'O': 1}, -1),
               ('M-2H', 1, {}, {'H': 2}, -2), ('2M-H', 2, {}, {'H': 1}, -1))


@dataclass
class Ion:
    label: str
    z: int
    mz: float
    counts: dict                                            # its formula, for the isotope pattern


@dataclass
class Info:
    atoms: int = 0
    formula: str = ''                                       # Hill order, e.g. C9H9NO5 (charge not included)
    formula_html: str = ''                                  # with <sub> numbers and the charge as <sup>
    counts: dict = field(default_factory=dict)              # {element: number}, hydrogens included
    charge: int = 0
    exact: float = 0.0                                      # monoisotopic mass of the neutral atoms (no electrons)
    average: float = 0.0
    analysis: str = ''                                      # elemental analysis: 'C, 51.19; H, 4.30; ...'
    ions: list = field(default_factory=list)                # [Ion]
    smiles: str = ''
    inchikey: str = ''
    errors: list = field(default_factory=list)              # human-readable problems (valence ...)
    fragments: int = 0
    repeat: str = ''                                        # polymer brackets with a letter (n): the masses are for 1
    radicals: int = 0                                       # unpaired electrons drawn (•): odd -> [M]•+ ions


def charge_text(z):
    if not z:
        return ''
    return ('' if abs(z) == 1 else str(abs(z))) + ('+' if z > 0 else '−')


def ion_label(base, z, radical=False):
    """[M+H]+, or [M]•+ for an odd-electron ion (a radical cation, as in EI mass spectra)."""
    return f'[{base}]{"•" if radical else ""}{charge_text(z)}'


def hill(counts):
    """Hill order: C, H, then alphabetical (alphabetical throughout when there is no carbon). A labelled isotope
    follows its element (C, ¹³C, H, D)."""
    def order(key):
        el, iso = split_key(key)
        return el, iso
    els = sorted(counts, key=order)
    if any(split_key(e)[0] == 'C' for e in counts):
        first = [e for e in els if split_key(e)[0] == 'C'] + [e for e in els if split_key(e)[0] == 'H']
        els = first + [e for e in els if e not in first]
    return els


def formula_text(counts):
    return ''.join(key_text(e) + (str(counts[e]) if counts[e] > 1 else '') for e in hill(counts) if counts[e] > 0)


def formula_html(counts):
    return ''.join(key_text(e, True) + (f'<sub>{counts[e]}</sub>' if counts[e] > 1 else '') for e in hill(counts)
                   if counts[e] > 0)


def to_mol(doc, atom_ids=None):
    """(RDKit molecule or None, [atom id in molecule order], error text). Nicknames go in as their first atom; use
    doc.expanded() first for the real chemistry."""
    block, order = doc.to_molblock(atom_ids)
    if not order:
        return None, order, ''
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if mol is None:
        return None, order, 'The structure could not be read.'
    for k, a in enumerate(order):                           # radicals and typed H counts (CH2••, CH3•, CH3+ ...)
        at = doc.atoms[a]
        if not at.label and (at.radical or at.hcount is not None):
            atom = mol.GetAtomWithIdx(k)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(doc.implicit_h(a))
            atom.SetNumRadicalElectrons(at.radical)
        if not at.label and at.iso:                             # 13C, D: [13CH3], [2H] in the SMILES
            mol.GetAtomWithIdx(k).SetIsotope(at.iso)
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignChiralTypesFromBondDirs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception as e:                                   # valence and kekulisation problems
        return None, order, str(e).split('\n')[0]
    return mol, order, ''


def expand_ids(grown, ids):
    """The atom ids of an expanded() copy that stand for the drawing's atoms ids."""
    out = set(ids)
    for a in ids:
        out |= grown.get(a, set())
    return out


def analyse(doc, atom_ids=None):
    info = Info()
    ids = [a for a in doc.atoms if atom_ids is None or a in atom_ids]
    info.atoms = len(ids)
    if not ids:
        return info
    info.fragments = len(doc.fragments(ids))
    for a in ids:
        problem = doc.valence_problem(a)
        if problem:
            info.errors.append(problem)
    full, grown = doc.expanded()
    repeats = doc.repeat_units(set(ids))
    ids = [a for a in full.atoms if a in expand_ids(grown, ids)]
    counts, mono, avg, charge, unknown = _tally(full, ids)
    info.radicals = sum(full.atoms[a].radical for a in ids)      # unpaired electrons (•)
    if unknown:
        info.errors.append('Unknown element: ' + ', '.join(sorted(set(unknown))))
    symbolic = []
    for rep, unit in repeats:                               # polymer brackets: the unit inside counts rep times
        unit = set(expand_ids(grown, unit)) & set(ids)
        uc, um, ua, uq, _u = _tally(full, unit)
        if rep.isdigit() and int(rep) >= 1:
            k = int(rep) - 1
            for e, c in uc.items():
                counts[e] = counts.get(e, 0) + k * c
            mono, avg, charge = mono + k * um, avg + k * ua, charge + k * uq
        else:
            symbolic.append((rep, uc))
    _masses(info, counts, charge, mono, avg)
    if symbolic:                                            # (C2H4)n: masses and ions are for n = 1
        ends = dict(counts)
        for _rep, uc in symbolic:
            for e, c in uc.items():
                ends[e] = ends.get(e, 0) - c
        ends = {e: c for e, c in ends.items() if c > 0}
        units = ''.join(f'({formula_text(uc)}){rep}' for rep, uc in symbolic)
        units_html = ''.join('(' + formula_html(uc) + f')<sub><i>{rep}</i></sub>' for rep, uc in symbolic)
        ends_html = formula_html(ends)
        info.formula = (formula_text(ends) if ends else '') + units
        info.formula_html = ends_html + units_html + (f'<sup>{charge_text(charge)}</sup>' if charge else '')
        info.repeat = ', '.join(rep for rep, _uc in symbolic)
    mol, _order, err = to_mol(full, set(ids))
    if mol is not None:
        _identifiers(info, mol)
    elif err and not info.errors:
        info.errors.append(err)
    return info


def partial_charges(doc):
    """{atom id: Gasteiger partial charge, its H included} of every atom drawn (a nickname: its group's sum). Atoms
    whose structure RDKit cannot read are left out."""
    from rdkit.Chem import AllChem
    full, grown = doc.expanded()
    owner = {a: a for a in doc.atoms}
    for a, more in grown.items():
        for g in more:
            owner[g] = a
    out = {}
    for frag in full.fragments():
        mol, order, _err = to_mol(full, set(frag))
        if mol is None:
            continue
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except Exception:
            continue
        for k, a in enumerate(order):
            atom = mol.GetAtomWithIdx(k)
            q = atom.GetDoubleProp('_GasteigerCharge')
            if atom.HasProp('_GasteigerHCharge'):
                q += atom.GetDoubleProp('_GasteigerHCharge')    # all its implicit H together
            if q != q or a not in owner:                    # NaN: an element Gasteiger has no parameters for
                continue
            out[owner[a]] = out.get(owner[a], 0.0) + q
    return out


def _tally(full, ids):
    """(element counts, monoisotopic mass, average mass, charge, unknown elements) of these atoms of the expanded
    drawing full, with their implicit hydrogens."""
    counts = {}
    mono = avg = 0.0
    charge = 0
    unknown = []
    for a in ids:
        at = full.atoms[a]
        key = iso_key(at.el, at.iso) if at.iso and not at.label else at.el
        try:
            m = key_mass(key)
            w = key_mass(key, average=True)
        except RuntimeError:
            unknown.append(at.el); continue
        if not m:                                           # an isotope the element does not have: the usual one
            key, m, w = at.el, key_mass(at.el), key_mass(at.el, True)
        hs = full.implicit_h(a)
        counts[key] = counts.get(key, 0) + 1
        if hs:
            counts['H'] = counts.get('H', 0) + hs
        mono += m + hs * H
        avg += w + hs * PT.GetAtomicWeight('H')
        charge += at.charge
    return counts, mono, avg, charge, unknown


def analyse_mol(mol):
    """analyse() for an RDKit molecule built elsewhere (a biopolymer sequence): formula, masses, adducts, SMILES."""
    info = Info()
    counts, mono, avg, charge = {}, 0.0, 0.0, 0
    for at in mol.GetAtoms():
        el, hs = at.GetSymbol(), at.GetTotalNumHs()
        if el == 'H':
            hs, el = 1, None
        if el:
            counts[el] = counts.get(el, 0) + 1
            mono += PT.GetMostCommonIsotopeMass(el)
            avg += PT.GetAtomicWeight(el)
            info.atoms += 1
        if hs:
            counts['H'] = counts.get('H', 0) + hs
            mono += hs * H
            avg += hs * PT.GetAtomicWeight('H')
        charge += at.GetFormalCharge()
    info.fragments = len(Chem.GetMolFrags(mol))
    _masses(info, counts, charge, mono, avg)
    _identifiers(info, mol)
    return info


def _identifiers(info, mol):
    try:
        info.smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
    except Exception:
        info.smiles = Chem.MolToSmiles(mol)
    try:
        info.inchikey = Chem.MolToInchiKey(mol) or ''
    except Exception:
        info.inchikey = ''


def _masses(info, counts, charge, mono, avg):
    """Formula, masses, elemental analysis and adduct ions of info from the element counts."""
    info.counts, info.charge, info.exact, info.average = counts, charge, mono, avg
    els = hill(counts)
    info.formula = formula_text(counts)
    info.formula_html = formula_html(counts)
    if charge or info.radicals:
        info.formula_html += f'<sup>{"•" * info.radicals}{charge_text(charge)}</sup>'
    if avg:
        by_el = {}                                          # by element (a labelled isotope counts with it)
        for e in els:
            el = split_key(e)[0]
            by_el[el] = by_el.get(el, 0.0) + key_mass(e, True) * counts[e]
        info.analysis = '; '.join(f'{e}, {100 * w / avg:.2f}' for e, w in by_el.items())
    if charge:                                              # a charged structure is an ion already
        info.ions.append(Ion(ion_label('M', charge, info.radicals % 2), charge, (mono - charge * ELECTRON) / abs(charge),
                             dict(counts)))
    for label, n, gain, loss, dz in ADDUCTS_POS + ADDUCTS_NEG:
        z = n * charge + dz
        if z == 0 or (charge and n > 1):
            continue
        if charge and (z > 0) != (charge > 0):
            continue
        ion = {e: n * c for e, c in counts.items()}
        for e, c in gain.items():
            ion[e] = ion.get(e, 0) + c
        for e, c in loss.items():
            ion[e] = ion.get(e, 0) - c
        if any(c < 0 for c in ion.values()):
            continue
        info.ions.append(Ion(ion_label(label, z, n * info.radicals % 2), z, (n * mono + _mono(gain) - _mono(loss) - z * ELECTRON) / abs(z), ion))


# ------------------------------------------------------------------------------------------------ isotope pattern
_isotopes = {}


def isotopes(el):
    """[(mass, abundance 0-1)] of the stable isotopes of an element (a labelled isotope key such as 13C or D: that one
    alone)."""
    if split_key(el)[1]:
        return [(key_mass(el), 1.0)]
    if el not in _isotopes:
        z = PT.GetAtomicNumber(el)
        base = PT.GetMostCommonIsotope(z)
        out = []
        for a in range(max(1, base - 8), base + 9):
            try:
                ab = PT.GetAbundanceForIsotope(z, a)
            except RuntimeError:
                ab = 0.0
            if ab > 0:
                out.append((PT.GetMassForIsotope(z, a), ab / 100.0))
        _isotopes[el] = out or [(PT.GetMostCommonIsotopeMass(el), 1.0)]
    return _isotopes[el]


def isotope_pattern(counts, z=0, min_percent=0.05, peaks=8):
    """The isotope peaks of a formula {element: n}, z the charge (0: neutral masses): [(m/z, % of the highest)],
    peaks one nominal mass apart (M, M+1, M+2 ...), each at its abundance-weighted mass."""
    dist = {0: (1.0, 0.0)}                                  # nominal offset -> (probability, probability x mass)
    for el, n in counts.items():
        iso = isotopes(el)
        if n <= 0:
            continue
        base_mass = iso[max(range(len(iso)), key=lambda k: iso[k][1])][0]
        for _ in range(n):
            new = {}
            for off, (p, pm) in dist.items():
                for mass, ab in iso:
                    k = off + round(mass - base_mass)
                    if k < -2 or k > peaks + 2:
                        continue
                    q, qm = new.get(k, (0.0, 0.0))
                    new[k] = (q + p * ab, qm + pm * ab + p * ab * mass)
            top = max(v[0] for v in new.values())
            dist = {k: v for k, v in new.items() if v[0] > top * 1e-7}
    top = max(v[0] for v in dist.values())
    out = []
    for k in sorted(dist):
        p, pm = dist[k]
        rel = 100 * p / top
        if rel < min_percent or k > peaks:
            continue
        mass = pm / p
        out.append(((mass - z * ELECTRON) / abs(z) if z else mass, rel))
    return out


# ------------------------------------------------------------------------------------------------ stereo labels
def cip_labels(doc):
    """({atom id: 'R' | 'S'}, {bond id: 'E' | 'Z'}) for the stereocentres and double bonds the drawing defines."""
    full, _grown = doc.expanded()
    atoms, bonds = {}, {}
    for frag in full.fragments():
        mol, order, _err = to_mol(full, set(frag))
        if mol is None:
            continue
        try:
            from rdkit.Chem import rdCIPLabeler
            rdCIPLabeler.AssignCIPLabels(mol)
        except Exception:
            continue
        for k, atom in enumerate(mol.GetAtoms()):
            if atom.HasProp('_CIPCode') and order[k] in doc.atoms:
                atoms[order[k]] = atom.GetProp('_CIPCode')
        for b in mol.GetBonds():
            if b.HasProp('_CIPCode'):
                i = doc.bond_between(order[b.GetBeginAtomIdx()], order[b.GetEndAtomIdx()])
                if i is not None:
                    bonds[i] = b.GetProp('_CIPCode')
    return atoms, bonds


# ------------------------------------------------------------------------------------------------ reading
def mol_from_text(text):
    """An RDKit molecule from SMILES (one line) or a MOL block; None if the text is neither."""
    text = text or ''
    if 'M  END' in text or 'V2000' in text or 'V3000' in text:
        return Chem.MolFromMolBlock(text, removeHs=False)  # as it is: its first line (the title) may be empty
    text = text.strip('\r\n ')
    if not text:
        return None
    first = text.split()[0] if len(text.splitlines()) == 1 else ''
    if not first or len(first) > 2000:
        return None
    return Chem.MolFromSmiles(first)


def mols_from_file(path):
    """The molecules in a .mol / .sdf / .smi / .cdx file (an error if none can be read)."""
    low = path.lower()
    if low.endswith('.sdf') or low.endswith('.sd'):
        mols = [m for m in Chem.SDMolSupplier(path, removeHs=False) if m is not None]
    elif low.endswith('.smi') or low.endswith('.smiles') or low.endswith('.txt'):
        with open(path, encoding='utf-8', errors='replace') as f:
            mols = [Chem.MolFromSmiles(line.split()[0]) for line in f if line.strip() and not line.startswith('#')]
        mols = [m for m in mols if m is not None]
    elif low.endswith('.cdx'):                               # ChemDraw's binary format: the structures only
        params = Chem.CDXMLParserParams()
        params.format = Chem.CDXMLFormat.CDX
        mols = [m for m in Chem.MolsFromCDXMLFile(path, params) if m is not None]
    else:
        m = Chem.MolFromMolFile(path, removeHs=False)
        mols = [m] if m is not None else []
    if not mols:
        raise ValueError('No structure could be read from this file.')
    return mols


# ------------------------------------------------------------------------------------------------ properties
def properties(doc, atom_ids=None):
    """[(name, value text, tooltip)] of calculated properties (RDKit; atom_ids None: everything drawn): Crippen logP
    and molar refractivity, topological polar surface area, H-bond donors / acceptors, rotatable bonds, rings, heavy
    atoms, fraction sp3 carbons, and whether Lipinski's rule of five holds."""
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    full, grown = doc.expanded()
    ids = expand_ids(grown, [a for a in doc.atoms if atom_ids is None or a in atom_ids])
    mol, _order, _err = to_mol(full, ids)
    if mol is None:
        return []
    logp, mr = Crippen.MolLogP(mol), Crippen.MolMR(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    violations = sum((mw > 500, logp > 5, hbd > 5, hba > 10))
    return [
        ('LogP', f'{logp:.2f}', 'Octanol / water partition coefficient (Crippen method)'),
        ('MR', f'{mr:.2f} cm³/mol', 'Molar refractivity (Crippen method)'),
        ('tPSA', f'{rdMolDescriptors.CalcTPSA(mol):.2f} Å²', 'Topological polar surface area'),
        ('H-bond donors', str(hbd), 'Lipinski donors (OH and NH)'),
        ('H-bond acceptors', str(hba), 'Lipinski acceptors (N and O)'),
        ('Rotatable bonds', str(rdMolDescriptors.CalcNumRotatableBonds(mol)), ''),
        ('Rings (aromatic)', f'{rdMolDescriptors.CalcNumRings(mol)} ({rdMolDescriptors.CalcNumAromaticRings(mol)})', ''),
        ('Heavy atoms', str(mol.GetNumHeavyAtoms()), 'Atoms other than hydrogen'),
        ('Fraction sp³ C', f'{rdMolDescriptors.CalcFractionCSP3(mol):.2f}', ''),
        ('Stereocentres', str(len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))), ''),
        ('Rule of five', 'passes' if violations <= 1 else f'{violations} violations',
         'Lipinski: MW ≤ 500, logP ≤ 5, ≤ 5 donors, ≤ 10 acceptors (one violation allowed)'),
    ]


ANALYSIS_KEYS = ('formula', 'exact', 'mz', 'mw', 'analysis')       # Show analysis / Paste all, in this order


def live_caption(doc, atom_ids, keys):
    """The text of a pasted analysis caption for these atoms: one line per key - an ANALYSIS_KEYS entry, 'name'
    (the IUPAC name) or 'prop:<name>' (a chem.properties row). None when the atoms make no structure."""
    info = analyse(doc, atom_ids)
    if not info.atoms:
        return None
    props = None
    lines = []
    per = f' ({info.repeat} = 1)' if info.repeat else ''
    for key in keys:
        if key == 'formula':
            lines.append(f'Chemical Formula: {info.formula}{charge_text(info.charge)}')
        elif key == 'exact':
            lines.append(f'Exact Mass: {info.exact:.2f}{per}')
        elif key == 'mw':
            lines.append(f'Molecular Weight: {info.average:.2f}{per}')
        elif key == 'mz':
            lines.append('m/z: ' + ', '.join(f'{m:.2f} ({p:.1f}%)' for m, p in isotope_pattern(info.counts, 0)[:4]))
        elif key == 'analysis':
            lines.append(f'Elemental Analysis: {info.analysis}')
        elif key == 'name':
            name = iupac_name(doc, atom_ids) if not info.errors else ''
            lines.append(name or '(no name found)')
        elif key == 'ms':                                   # a fragment: formula and m/z (or neutral mass)
            import ms
            lines.append(ms.describe(doc, atom_ids, chr(10)))
        elif key.startswith('prop:'):
            if props is None:
                props = {n: v for n, v, _tip in properties(doc, atom_ids)}
            if key[5:] in props:
                lines.append(f'{key[5:]}: {props[key[5:]]}')
    return '\n'.join(lines)


# ------------------------------------------------------------------------------------------------ names
def nest_enclosures(name):
    """IUPAC nesting of enclosing marks, innermost out: ( ), [ ], { }, ( ) ... A name that already has [ ] or { }
    is left alone. '2-(4-(2-methylpropyl)phenyl)propanoic acid' -> '2-[4-(2-methylpropyl)phenyl]propanoic acid'."""
    if '[' in name or '{' in name or name.count('(') != name.count(')'):
        return name
    pairs, stack = {}, []
    for k, ch in enumerate(name):
        if ch == '(':
            stack.append(k)
        elif ch == ')':
            if not stack:
                return name
            pairs[stack.pop()] = k
    height = {}                                             # how many levels of enclosures sit inside each pair

    def measure(open_k):
        if open_k not in height:
            inner = [o for o in pairs if open_k < o < pairs[open_k]]
            height[open_k] = 1 + max((measure(o) for o in inner), default=-1)
        return height[open_k]
    for o in pairs:
        measure(o)
    out = list(name)
    marks = ('()', '[]', '{}')
    for o, c in pairs.items():
        m = marks[height[o] % 3]
        out[o], out[c] = m[0], m[1]
    return ''.join(out)


def iupac_name(doc, atom_ids=None):
    """A systematic (IUPAC) name for the structure, made offline by Openclatura (rule-based, IUPAC 2013 rules), or
    '' when it cannot name it. Nicknames count in full; stereo comes from the wedges and double-bond geometry."""
    full, grown = doc.expanded()
    ids = expand_ids(grown, [a for a in doc.atoms if atom_ids is None or a in atom_ids])
    mol, _order, _err = to_mol(full, ids)
    if mol is None:
        return ''
    try:
        import openclatura
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        name = openclatura.name_smiles(smiles)
    except Exception:
        import logging
        logging.getLogger('kekulex').warning('no name for %s', Chem.MolToSmiles(mol), exc_info=True)
        return ''
    import logging
    logging.getLogger('kekulex').debug('named %s: %s', Chem.MolToSmiles(mol), name)
    return nest_enclosures(name.strip()) if name else ''
