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
    return sum(PT.GetMostCommonIsotopeMass(el) * n for el, n in formula.items())


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


def charge_text(z):
    if not z:
        return ''
    return ('' if abs(z) == 1 else str(abs(z))) + ('+' if z > 0 else '−')


def ion_label(base, z):
    return f'[{base}]{charge_text(z)}'


def hill(counts):
    """Hill order: C, H, then alphabetical (alphabetical throughout when there is no carbon)."""
    els = sorted(counts)
    if 'C' in counts:
        els = ['C'] + (['H'] if 'H' in counts else []) + [e for e in els if e not in ('C', 'H')]
    return els


def formula_text(counts):
    return ''.join(e + (str(counts[e]) if counts[e] > 1 else '') for e in hill(counts) if counts[e] > 0)


def to_mol(doc, atom_ids=None):
    """(RDKit molecule or None, [atom id in molecule order], error text). Nicknames go in as their first atom; use
    doc.expanded() first for the real chemistry."""
    block, order = doc.to_molblock(atom_ids)
    if not order:
        return None, order, ''
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if mol is None:
        return None, order, 'The structure could not be read.'
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
        if doc.valence_error(a):
            at = doc.atoms[a]
            info.errors.append(f'{at.el} with {doc.bond_sum(a)} bonds{" and charge " + charge_text(at.charge) if at.charge else ""} is over its valence.')
    full, grown = doc.expanded()
    ids = [a for a in full.atoms if a in expand_ids(grown, ids)]
    counts = {}
    mono = avg = 0.0
    charge = 0
    unknown = []
    for a in ids:
        at = full.atoms[a]
        try:
            m = PT.GetMostCommonIsotopeMass(at.el)
            w = PT.GetAtomicWeight(at.el)
        except RuntimeError:
            unknown.append(at.el); continue
        hs = full.implicit_h(a)
        counts[at.el] = counts.get(at.el, 0) + 1
        if hs:
            counts['H'] = counts.get('H', 0) + hs
        mono += m + hs * H
        avg += w + hs * PT.GetAtomicWeight('H')
        charge += at.charge
    if unknown:
        info.errors.append('Unknown element: ' + ', '.join(sorted(set(unknown))))
    info.counts, info.charge, info.exact, info.average = counts, charge, mono, avg
    els = hill(counts)
    info.formula = formula_text(counts)
    info.formula_html = ''.join(e + (f'<sub>{counts[e]}</sub>' if counts[e] > 1 else '') for e in els)
    if charge:
        info.formula_html += f'<sup>{charge_text(charge)}</sup>'
    if avg:
        info.analysis = '; '.join(f'{e}, {100 * PT.GetAtomicWeight(e) * counts[e] / avg:.2f}' for e in els)
    if charge:                                              # a charged structure is an ion already
        info.ions.append(Ion(ion_label('M', charge), charge, (mono - charge * ELECTRON) / abs(charge), dict(counts)))
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
        info.ions.append(Ion(ion_label(label, z), z, (n * mono + _mono(gain) - _mono(loss) - z * ELECTRON) / abs(z), ion))
    mol, _order, err = to_mol(full, set(ids))
    if mol is not None:
        try:
            info.smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        except Exception:
            info.smiles = Chem.MolToSmiles(mol)
        try:
            info.inchikey = Chem.MolToInchiKey(mol) or ''
        except Exception:
            info.inchikey = ''
    elif err and not info.errors:
        info.errors.append(err)
    return info


# ------------------------------------------------------------------------------------------------ isotope pattern
_isotopes = {}


def isotopes(el):
    """[(mass, abundance 0-1)] of the stable isotopes of an element."""
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
