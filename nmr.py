"""
nmr.py - 1H and 13C NMR shift prediction with CSP5 (MIT; a graph neural network trained on measured spectra, solvent
aware, with an uncertainty for each shift).

    if nmr.available():
        result = nmr.predict(doc, atom_ids, '1H', 'cdcl3')   # Result: signals (shift, +/-, atoms, H count, multiplicity)
        nmr.journal_text(result)                               # '1H NMR (CDCl3, predicted) δ 8.22 (d, 2H), ...'

CSP5 is not built for Windows by its authors: the editor uses its source with the optional peak-matching helpers left
out (tools/patch_csp5.py); prediction does not need them. The first prediction loads the model (several seconds).
Multiplicities are first-order, with expected coupling constants (see couplings()): vicinal sp3 7 Hz (Karplus on an
RDKit 3D conformer for ring bonds: axial-axial ~10, axial-equatorial ~3), aromatic ortho 8 / meta 2 (5-ring 3.5, next
to a ring N 5), alkene cis 10-11 / trans 16, sp2-sp2 single bond 11, aldehyde 2.5; OH / NH / SH as broad singlets.
"""
import importlib.util
from dataclasses import dataclass, field

from rdkit import Chem

import chem

SOLVENTS = (('CDCl₃', 'cdcl3'), ('DMSO-d₆', 'dmso-d6'), ('CD₃OD', 'cd3od'), ('D₂O', 'd2o'), ('Acetone-d₆', 'acetone-d6'),
            ('CD₂Cl₂', 'cd2cl2'), ('C₆D₆', 'c6d6'), ('CD₃CN', 'cd3cn'), ('THF-d₈', 'thf-d8'),
            ('Toluene-d₈', 'toluene-d8'), ('DMF-d₇', 'dmf-d7'))
SOLVENT_NAMES = {key: name for name, key in SOLVENTS}
MULTIPLETS = {0: 's', 1: 'd', 2: 't', 3: 'q', 4: 'quint', 5: 'sext', 6: 'sept'}
J_FREE = 7.0                                                # vicinal H-C-C-H across a freely turning bond, Hz


def available():
    return importlib.util.find_spec('csp5') is not None and importlib.util.find_spec('torch') is not None


@dataclass
class Signal:
    shift: float                                            # ppm (average over equivalent atoms)
    std: float                                              # uncertainty, ppm (0 when not known)
    atoms: list                                             # the drawing's atom ids (carbons, or the H-bearing atoms)
    n_h: int = 0                                            # 1H: hydrogens in the signal
    couplings: list = field(default_factory=list)           # 1H: [(J Hz, number of hydrogens)], largest J first
    multiplicity: str = ''                                  # 1H: s, d, t, q, dd, dt, td, ddd ..., m, br s

    def j_text(self):
        """'8.4, 2.1' (the coupling constants of a resolved multiplet), '' for s / br s / m."""
        if self.multiplicity in ('s', 'br s', 'm', ''):
            return ''
        return ', '.join(f'{j:.1f}' for j, _n in self.couplings)


@dataclass
class Result:
    nucleus: str
    solvent: str
    signals: list = field(default_factory=list)             # highest shift first
    warnings: list = field(default_factory=list)


def _mol_for(doc, atom_ids):
    """(RDKit molecule with explicit H, {its atom index: drawing atom id}) for the drawn structure; nicknames count
    in full and their atoms point at the nickname's atom."""
    full, grown = doc.expanded()
    owner = {}
    for a, new in grown.items():
        for b in new:
            owner[b] = a
    ids = chem.expand_ids(grown, [a for a in doc.atoms if atom_ids is None or a in atom_ids])
    mol, order, err = chem.to_mol(full, ids)
    if mol is None:
        raise ValueError(err or 'The structure could not be read.')
    return mol, {k: owner.get(a, a) for k, a in enumerate(order)}


def predict(doc, atom_ids, nucleus='1H', solvent='cdcl3'):
    """Predicted signals of the structure (atom_ids None: everything drawn)."""
    import csp5
    mol, to_drawing = _mol_for(doc, atom_ids)
    if Chem.MolToSmiles(mol).count('.'):
        raise ValueError('Select one structure for NMR prediction (the drawing holds several).')
    smiles = Chem.MolToSmiles(mol)
    model = 'CSP5q-13C' if nucleus == '13C' else 'CSP5q-1H'
    r = csp5.predict_smiles([smiles], nucleus=nucleus, model_name=model, solvent=solvent, random_seed=1)
    if r.failures or not r.molecule_records:
        raise ValueError('CSP5 could not predict this structure.')
    rec = r.molecule_records[0].mol
    match = rec.GetSubstructMatch(mol)                      # drawn atom k -> CSP5's atom match[k]
    if not match:
        raise ValueError('The prediction could not be matched to the drawing.')
    rec_to_drawing = {rec_idx: to_drawing[k] for k, rec_idx in enumerate(match)}
    ranks = list(Chem.CanonicalRankAtoms(rec, breakTies=False))
    groups = {}                                             # symmetry class -> [(shift, std, drawing atom, rec atom, H)]
    for _i, row in r.predictions.iterrows():
        idx = int(row['atom_index'])
        heavy = idx
        if nucleus == '1H':
            nbrs = rec.GetAtomWithIdx(idx).GetNeighbors()
            if not nbrs:
                continue
            heavy = nbrs[0].GetIdx()
        if heavy not in rec_to_drawing:
            continue
        std = float(row.get('shift_std_ppm', 0.0) or 0.0)
        std = 0.0 if std != std else std                    # NaN: no uncertainty given
        groups.setdefault(ranks[idx], []).append((float(row['shift_ppm']), std, rec_to_drawing[heavy], heavy, idx))
    result = Result(nucleus, solvent, warnings=list(r.warnings or []))
    if nucleus == '1H':                                     # couplings on the drawn molecule: it keeps the wedges
        hmol = Chem.AddHs(mol)
        hmatch = rec.GetSubstructMatch(hmol)
        if len(hmatch) == hmol.GetNumAtoms():
            rec_to_h = {rec_idx: k for k, rec_idx in enumerate(hmatch)}
            jmol, jranks = hmol, list(Chem.CanonicalRankAtoms(hmol, breakTies=False))
        else:
            rec_to_h, jmol, jranks = None, rec, ranks
        conf = _conformer(jmol)
    for members in groups.values():
        shift = sum(m[0] for m in members) / len(members)
        std = max(m[1] for m in members)
        atoms = sorted({m[2] for m in members})
        sig = Signal(round(shift, 2 if nucleus == '1H' else 1), round(std, 2), atoms)
        if nucleus == '1H':
            sig.n_h = len(members)
            heavy = members[0][3]
            if rec.GetAtomWithIdx(heavy).GetSymbol() in ('O', 'N', 'S'):
                sig.multiplicity = 'br s'                  # exchangeable: usually no resolved coupling
            else:
                hs = [rec_to_h.get(m[4], m[4]) if rec_to_h is not None else m[4] for m in members]
                sig.couplings = couplings(jmol, hs, jranks, conf)
                sig.multiplicity = multiplet_name(sig.couplings)
        result.signals.append(sig)
    result.signals.sort(key=lambda s: -s.shift)
    return result


EXCHANGING_SOLVENTS = ('d2o', 'cd3od')


def exchanged(result, signal):
    """True for an OH / NH / SH signal in D2O or CD3OD: the hydrogen exchanges with deuterium and is not seen."""
    return result.nucleus == '1H' and signal.multiplicity == 'br s' and result.solvent in EXCHANGING_SOLVENTS


def journal_text(result, mhz=400):
    """The predicted spectrum the way experimental sections give it."""
    solvent = SOLVENT_NAMES.get(result.solvent, result.solvent)
    if result.nucleus == '1H':
        parts = [f'{s.shift:.2f} ({s.multiplicity}, ' + (f'J = {s.j_text()} Hz, ' if s.j_text() else '') + f'{s.n_h}H)'
                 for s in result.signals if not exchanged(result, s)]
        return f'¹H NMR ({mhz} MHz, {solvent}, predicted) δ ' + ', '.join(parts) + '.'
    shifts = []
    for s in result.signals:
        shifts.append(f'{s.shift:.1f}' + (f' ({len(s.atoms)}C)' if len(s.atoms) > 1 else ''))
    return f'¹³C NMR ({mhz // 4} MHz, {solvent}, predicted) δ ' + ', '.join(shifts) + '.'


def spectrum(result, mhz=400, points=6000):
    """(ppm array, intensity array) of a simulated spectrum: Lorentzian lines; 1H multiplets split by their
    couplings with binomial heights, 13C one line per carbon (equivalent carbons add up)."""
    import numpy as np
    if result.nucleus == '1H':
        lo, hi = min([0.0] + [s.shift - 0.5 for s in result.signals]), max([10.0] + [s.shift + 0.5 for s in result.signals])
        width = 1.2 / mhz                                   # 1.2 Hz lines
    else:
        lo, hi = min([0.0] + [s.shift - 10 for s in result.signals]), max([200.0] + [s.shift + 10 for s in result.signals])
        width = 0.15
    x = np.linspace(lo, hi, points if result.nucleus == '13C' else points * 4)
    y = np.zeros_like(x)
    for s in result.signals:
        if exchanged(result, s):
            continue
        if result.nucleus == '1H':
            lw = width * (6 if s.multiplicity == 'br s' else 1)
            for pos, w in lines(s, mhz):
                y += s.n_h * w * (lw / 2) ** 2 / ((x - pos) ** 2 + (lw / 2) ** 2)
        else:
            y += len(s.atoms) * (width / 2) ** 2 / ((x - s.shift) ** 2 + (width / 2) ** 2)
    return x, y


def lines(signal, mhz=400):
    """[(ppm, relative height)] of a first-order multiplet: split by each coupling in turn (binomial), heights sum to 1."""
    from math import comb
    out = [(signal.shift, 1.0)]
    if signal.multiplicity == 'br s':
        return out
    for j, n in signal.couplings:
        step = j / mhz
        out = [(pos + (k - n / 2) * step, w * comb(n, k) / 2 ** n) for pos, w in out for k in range(n + 1)]
    return out


# ------------------------------------------------------------------------------------------------- coupling constants
def _conformer(mol):
    """The lowest-energy (MMFF) of a few 3D conformers of mol (explicit H), for the ring dihedrals and the alkene
    geometry, or None."""
    from rdkit.Chem import AllChem
    m = Chem.Mol(mol)
    try:
        ids = list(AllChem.EmbedMultipleConfs(m, numConfs=12, randomSeed=0xF00D, numThreads=0))
        if not ids:
            return None
        energies = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=1000, numThreads=0)
        best = min(range(len(ids)), key=lambda k: energies[k][1] if energies[k][0] == 0 else energies[k][1] + 1e6)
        return m.GetConformer(ids[best])
    except Exception:
        return None


def _karplus(phi_deg):
    """Vicinal 3J(H,H) in Hz from the H-C-C-H dihedral (Karplus, sp3-sp3): 180 deg ~10, 60 deg ~3, 90 deg ~1.4."""
    from math import cos, radians
    c = cos(radians(phi_deg))
    return 7.76 * c * c - 1.10 * c + 1.40


def _ring_size(mol, a, b):
    """Size of the smallest ring holding the bond a-b (0 when none)."""
    idx = mol.GetBondBetweenAtoms(a, b).GetIdx()
    sizes = [len(r) for r in mol.GetRingInfo().BondRings() if idx in r]
    return min(sizes) if sizes else 0


def _is_aldehyde_c(atom):
    return atom.GetSymbol() == 'C' and any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(atom).GetSymbol()
                                           == 'O' for b in atom.GetBonds())


def _vicinal(mol, h, c1, c2, h2, conf):
    """Expected 3J between h (on c1) and h2 (on c2), Hz."""
    from rdkit.Chem import rdMolTransforms
    a1, a2 = mol.GetAtomWithIdx(c1), mol.GetAtomWithIdx(c2)
    bond = mol.GetBondBetweenAtoms(c1, c2)
    sp2 = Chem.HybridizationType.SP2
    if _is_aldehyde_c(a1) or _is_aldehyde_c(a2):
        return 2.5                                          # CHO to the next CH
    if bond.GetIsAromatic():
        if _ring_size(mol, c1, c2) == 5:                    # furan 1.8 / 3.4, pyrrole 2.6 / 3.5, thiophene 5 / 3.5
            hetero = {nb.GetSymbol() for a in (a1, a2) for nb in a.GetNeighbors() if nb.GetIsAromatic()} - {'C'}
            return 5.0 if 'S' in hetero else 1.8 if 'O' in hetero else 2.6 if hetero else 3.5
        next_to_n = any(nb.GetSymbol() == 'N' and nb.GetIsAromatic() for a in (a1, a2) for nb in a.GetNeighbors())
        return 5.0 if next_to_n else 8.0                    # pyridine H2-H3 ~5, benzene ortho ~8
    if bond.GetBondType() == Chem.BondType.DOUBLE:
        size = _ring_size(mol, c1, c2)
        if size and size <= 8:
            return 10.0                                     # cis in a ring
        if conf is not None:
            phi = abs(rdMolTransforms.GetDihedralDeg(conf, h, c1, c2, h2))
            return 11.0 if phi < 90 else 16.0               # cis / trans
        return 16.0
    if a1.GetHybridization() == sp2 and a2.GetHybridization() == sp2:
        return 11.0                                         # H-C=C-C=C-H across the single bond (s-trans)
    if conf is not None and _ring_size(mol, c1, c2) and sp2 not in (a1.GetHybridization(), a2.GetHybridization()):
        return round(_karplus(rdMolTransforms.GetDihedralDeg(conf, h, c1, c2, h2)), 1)
    return J_FREE


def _same_ring(mol, size, *atoms):
    """True when the atoms lie in one ring of this size."""
    ids = {a.GetIdx() for a in atoms}
    return any(len(r) == size and ids <= set(r) for r in mol.GetRingInfo().AtomRings())


def couplings(mol, h, ranks, conf=None):
    """[(J Hz, number of hydrogens)] of hydrogen h (molecule with explicit H; or a list of the equivalent hydrogens of
    one signal, averaged - a ring CH2), largest first: vicinal 3J to the hydrogens on the neighbouring carbons and 4J
    meta (2 Hz) on aromatic rings. Hydrogens equivalent to h, those on O / N / S (exchanging) and couplings under 1 Hz
    are left out; groups within 1 Hz of each other count as one (a CH2 between CH3 and CH2 is a sextet, as usually
    reported)."""
    hs = list(h) if isinstance(h, (list, tuple)) else [h]
    per_h = [_partners(mol, x, ranks, conf) for x in hs]
    partners = {}                                           # symmetry class -> [J] over all the hydrogens
    for p in per_h:
        for rank, js in p.items():
            partners.setdefault(rank, []).extend(js)
    groups = sorted(((sum(js) / len(js), max(1, round(len(js) / len(hs)))) for js in partners.values()
                     if sum(js) / len(js) >= 1.0), reverse=True)
    merged = []
    for j, n in groups:
        if merged and merged[-1][0] - j < 1.0:
            pj, pn = merged[-1]
            merged[-1] = ((pj * pn + j * n) / (pn + n), pn + n)
        else:
            merged.append((j, n))
    return [(round(j, 1), n) for j, n in merged]


def _partners(mol, h, ranks, conf):
    """{symmetry class of the coupled hydrogens: [J Hz, one per hydrogen]} of hydrogen h."""
    c1 = mol.GetAtomWithIdx(h).GetNeighbors()[0]
    partners = {}
    for c2 in c1.GetNeighbors():
        if c2.GetSymbol() in ('H', 'O', 'N', 'S'):
            continue
        for h2 in c2.GetNeighbors():
            if h2.GetSymbol() != 'H' or ranks[h2.GetIdx()] == ranks[h]:
                continue
            partners.setdefault(ranks[h2.GetIdx()], []).append(_vicinal(mol, h, c1.GetIdx(), c2.GetIdx(),
                                                                        h2.GetIdx(), conf))
        if c1.GetIsAromatic() and c2.GetIsAromatic():       # meta: c1 - c2 - c3 in one aromatic ring
            for c3 in c2.GetNeighbors():
                if c3.GetIdx() == c1.GetIdx() or not c3.GetIsAromatic() or not _same_ring(mol, 6, c1, c2, c3):
                    continue
                for h3 in c3.GetNeighbors():
                    if h3.GetSymbol() == 'H' and ranks[h3.GetIdx()] != ranks[h]:
                        partners.setdefault(ranks[h3.GetIdx()], []).append(2.0)
    return partners


def multiplet_name(groups):
    """First-order name of a multiplet from its couplings: s, d, t, q ... (one group), dd, dt, td, ddd ... (up to three
    groups of at most three hydrogens each), else m."""
    if not groups:
        return 's'
    if len(groups) == 1:
        return MULTIPLETS.get(groups[0][1], 'm')
    if len(groups) > 3 or any(n > 3 for _j, n in groups):
        return 'm'
    return ''.join(MULTIPLETS[n] for _j, n in groups)
