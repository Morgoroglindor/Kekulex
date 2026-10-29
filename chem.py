"""
chem.py - what the drawn structure is: formula, masses, the m/z of the usual LC-MS adducts, SMILES and InChIKey (RDKit).

    info = chem.analyse(doc, atom_ids)      # Info (atom_ids None: everything drawn)
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


# (label, molecules, neutral atoms gained (+) or lost (-), charge gained). The label's 'M' takes the count.
ADDUCTS_POS = (('M+H', 1, _mono({'H': 1}), 1), ('M+NH4', 1, _mono({'N': 1, 'H': 4}), 1),
               ('M+Na', 1, _mono({'Na': 1}), 1), ('M+K', 1, _mono({'K': 1}), 1),
               ('M+H-H2O', 1, _mono({'H': 1}) - _mono({'H': 2, 'O': 1}), 1),
               ('M+2H', 1, _mono({'H': 2}), 2), ('2M+H', 2, _mono({'H': 1}), 1), ('2M+Na', 2, _mono({'Na': 1}), 1))
ADDUCTS_NEG = (('M-H', 1, -_mono({'H': 1}), -1), ('M+Cl', 1, _mono({'Cl': 1}), -1),
               ('M+HCOO', 1, _mono({'C': 1, 'H': 1, 'O': 2}), -1), ('M+CH3COO', 1, _mono({'C': 2, 'H': 3, 'O': 2}), -1),
               ('M-H-H2O', 1, -_mono({'H': 3, 'O': 1}), -1), ('M-2H', 1, -_mono({'H': 2}), -2),
               ('2M-H', 2, -_mono({'H': 1}), -1))


@dataclass
class Info:
    atoms: int = 0
    formula: str = ''                                       # Hill order, e.g. C9H9NO5 (charge not included)
    formula_html: str = ''                                  # with <sub> numbers and the charge as <sup>
    charge: int = 0
    exact: float = 0.0                                      # monoisotopic mass of the neutral atoms (no electrons)
    average: float = 0.0
    ions: list = field(default_factory=list)                # [(label, z, m/z)]
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


def to_mol(doc, atom_ids=None):
    """(RDKit molecule or None, [atom id in molecule order], error text)."""
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
    counts = {}
    mono = avg = 0.0
    charge = 0
    unknown = []
    for a in ids:
        at = doc.atoms[a]
        try:
            m = PT.GetMostCommonIsotopeMass(at.el)
            w = PT.GetAtomicWeight(at.el)
        except RuntimeError:
            unknown.append(at.el); continue
        hs = doc.implicit_h(a)
        counts[at.el] = counts.get(at.el, 0) + 1
        if hs:
            counts['H'] = counts.get('H', 0) + hs
        mono += m + hs * H
        avg += w + hs * PT.GetAtomicWeight('H')
        charge += at.charge
    if unknown:
        info.errors.append('Unknown element: ' + ', '.join(sorted(set(unknown))))
    info.charge, info.exact, info.average = charge, mono, avg
    els = hill(counts)
    info.formula = ''.join(e + (str(counts[e]) if counts[e] > 1 else '') for e in els)
    info.formula_html = ''.join(e + (f'<sub>{counts[e]}</sub>' if counts[e] > 1 else '') for e in els)
    if charge:
        info.formula_html += f'<sup>{charge_text(charge)}</sup>'
    if charge:                                              # a charged structure is an ion already
        info.ions.append((ion_label('M', charge), charge, (mono - charge * ELECTRON) / abs(charge)))
    for label, n, dm, dz in ADDUCTS_POS + ADDUCTS_NEG:
        z = n * charge + dz
        if z == 0 or (charge and n > 1):
            continue
        if charge and (z > 0) != (charge > 0):
            continue
        info.ions.append((ion_label(label, z), z, (n * mono + dm - z * ELECTRON) / abs(z)))
    mol, _order, err = to_mol(doc, ids)
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
    """The molecules in a .mol / .sdf / .smi file (an error if none can be read)."""
    low = path.lower()
    if low.endswith('.sdf') or low.endswith('.sd'):
        mols = [m for m in Chem.SDMolSupplier(path, removeHs=False) if m is not None]
    elif low.endswith('.smi') or low.endswith('.smiles') or low.endswith('.txt'):
        with open(path, encoding='utf-8', errors='replace') as f:
            mols = [Chem.MolFromSmiles(line.split()[0]) for line in f if line.strip() and not line.startswith('#')]
        mols = [m for m in mols if m is not None]
    else:
        m = Chem.MolFromMolFile(path, removeHs=False)
        mols = [m] if m is not None else []
    if not mols:
        raise ValueError('No structure could be read from this file.')
    return mols
