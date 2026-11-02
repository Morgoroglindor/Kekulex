"""
md_export.py - one conformer of a 3D model as input for molecular dynamics programs: PDB, MOL2 (with Gasteiger charges;
what antechamber / ACPYPE, CGenFF and LigParGen read to make a force-field topology), SDF, XYZ, GROMACS .gro, AMBER
restart (.inpcrd / .rst7) and a LAMMPS data file. These hold coordinates and connectivity: force-field parameters must
come from those tools, as no general force field can be derived here.
"""
from collections import Counter

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

BOX_PADDING = 10.0                                          # Å of vacuum on each side in .gro and LAMMPS boxes
RESIDUE = 'LIG'
FORMATS = (('PDB, with bonds (*.pdb)', 'pdb'), ('MOL2, with Gasteiger charges: antechamber, ACPYPE, CGenFF (*.mol2)', 'mol2'),
           ('GROMACS coordinates (*.gro)', 'gro'), ('AMBER coordinates (*.inpcrd *.rst7)', 'inpcrd'),
           ('LAMMPS data, atom style full (*.data *.lmp)', 'lammps'), ('SDF, this conformer (*.sdf)', 'sdf'),
           ('SDF, all conformers (*.sdf)', 'sdf_all'), ('XYZ (*.xyz)', 'xyz'))
NOTE = ('The file holds coordinates and bonds (MOL2 and LAMMPS also Gasteiger charges). For a simulation, get the force '
        'field parameters from a tool such as ACPYPE / antechamber (GAFF), CGenFF or LigParGen, which read the MOL2 file.')


def atom_names(mol):
    """Unique names as MD programs want them: C1, C2, O1, H1 ..."""
    seen = Counter()
    out = []
    for a in mol.GetAtoms():
        seen[a.GetSymbol()] += 1
        out.append(f'{a.GetSymbol()}{seen[a.GetSymbol()]}'[:4])
    return out


def charges(mol):
    """Gasteiger charges summing to the molecule's formal charge; (charges, False) as zeros if they cannot be had."""
    m = Chem.Mol(mol)
    try:
        AllChem.ComputeGasteigerCharges(m, throwOnParamFailure=True)
        q = np.array([a.GetDoubleProp('_GasteigerCharge') for a in m.GetAtoms()])
    except (RuntimeError, ValueError, KeyError):
        return np.zeros(mol.GetNumAtoms()), False
    if not np.isfinite(q).all():
        return np.zeros(mol.GetNumAtoms()), False
    return q + (Chem.GetFormalCharge(mol) - q.sum()) / len(q), True


def sybyl_type(a):
    """The SYBYL atom type MOL2 readers expect."""
    sym, hyb = a.GetSymbol(), a.GetHybridization()
    sp, sp2, sp3 = Chem.HybridizationType.SP, Chem.HybridizationType.SP2, Chem.HybridizationType.SP3
    orders = [b.GetBondTypeAsDouble() for b in a.GetBonds()]
    heavy_o = [n for n in a.GetNeighbors() if n.GetSymbol() == 'O' and n.GetDegree() == 1]
    if sym == 'C':
        if a.GetIsAromatic():
            return 'C.ar'
        if hyb == sp2:
            cation = a.GetDegree() == 3 and sum(n.GetSymbol() == 'N' and n.GetDegree() == 3 for n in a.GetNeighbors()) == 3
            return 'C.cat' if cation else 'C.2'
        return 'C.1' if hyb == sp else 'C.3'
    if sym == 'N':
        if a.GetIsAromatic():
            return 'N.ar'
        if a.GetFormalCharge() > 0 and a.GetDegree() == 4:
            return 'N.4'
        if hyb == sp:
            return 'N.1'
        if any(_is_carbonyl_c(n) for n in a.GetNeighbors()) and 3.0 not in orders and 2.0 not in orders:
            return 'N.am'
        if 2.0 in orders:
            return 'N.2'
        return 'N.pl3' if hyb == sp2 else 'N.3'
    if sym == 'O':
        nbrs = a.GetNeighbors()
        if a.GetDegree() == 1 and nbrs[0].GetSymbol() in ('C', 'P'):
            terminal = [n for n in nbrs[0].GetNeighbors() if n.GetSymbol() == 'O' and n.GetDegree() == 1]
            if len(terminal) >= 2:
                return 'O.co2'
        return 'O.2' if hyb in (sp2, sp) and a.GetDegree() == 1 and not a.GetFormalCharge() else 'O.3'
    if sym == 'S':
        if len(heavy_o) == 1 and a.GetDegree() == 3:
            return 'S.O'
        if len(heavy_o) >= 2:
            return 'S.O2'
        return 'S.2' if 2.0 in orders or a.GetIsAromatic() else 'S.3'
    if sym == 'P':
        return 'P.3'
    return sym if hyb != sp3 or sym in ('H', 'F', 'Cl', 'Br', 'I', 'Na', 'K', 'Li', 'Ca', 'Mg', 'Zn') else sym


def _is_carbonyl_c(a):
    return a.GetSymbol() == 'C' and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(a).GetSymbol() in ('O', 'S')
                                        for b in a.GetBonds())


def _pdb(mol, xyz):
    names = atom_names(mol)
    out = []
    for i, a in enumerate(mol.GetAtoms()):
        name = names[i] if len(a.GetSymbol()) == 2 or len(names[i]) == 4 else ' ' + names[i]
        q = a.GetFormalCharge()
        charge = f'{abs(q)}{"+" if q > 0 else "-"}' if q else '  '
        out.append(f'HETATM{i + 1:5d} {name:<4s} {RESIDUE} A   1    {xyz[i][0]:8.3f}{xyz[i][1]:8.3f}{xyz[i][2]:8.3f}'
                   f'  1.00  0.00          {a.GetSymbol().upper():>2s}{charge}')
    for a in mol.GetAtoms():
        nbrs = sorted(b.GetOtherAtomIdx(a.GetIdx()) + 1 for b in a.GetBonds()
                      for _ in range(max(1, int(b.GetBondTypeAsDouble()))))
        for k in range(0, len(nbrs), 4):
            out.append(f'CONECT{a.GetIdx() + 1:5d}' + ''.join(f'{n:5d}' for n in nbrs[k:k + 4]))
    return '\n'.join(out + ['END', ''])


def _mol2(mol, xyz, name):
    names = atom_names(mol)
    q, _ok = charges(mol)
    out = ['@<TRIPOS>MOLECULE', name or RESIDUE, f'{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0', 'SMALL', 'GASTEIGER', '',
           '@<TRIPOS>ATOM']
    for i, a in enumerate(mol.GetAtoms()):
        out.append(f'{i + 1:7d} {names[i]:<6s}{xyz[i][0]:10.4f}{xyz[i][1]:10.4f}{xyz[i][2]:10.4f} {sybyl_type(a):<6s}'
                   f'   1 {RESIDUE} {q[i]:9.4f}')
    out.append('@<TRIPOS>BOND')
    for k, b in enumerate(mol.GetBonds()):
        kind = {Chem.BondType.SINGLE: '1', Chem.BondType.DOUBLE: '2', Chem.BondType.TRIPLE: '3',
                Chem.BondType.AROMATIC: 'ar'}.get(b.GetBondType(), '1')
        out.append(f'{k + 1:6d} {b.GetBeginAtomIdx() + 1:5d} {b.GetEndAtomIdx() + 1:5d} {kind}')
    out += ['@<TRIPOS>SUBSTRUCTURE', f'     1 {RESIDUE}         1 TEMP              0 ****  ****    0 ROOT', '']
    return '\n'.join(out)


def _gro(mol, xyz, name):
    names = atom_names(mol)
    p = (xyz - xyz.min(axis=0) + BOX_PADDING) / 10.0        # nm, in the box
    box = (xyz.max(axis=0) - xyz.min(axis=0) + 2 * BOX_PADDING) / 10.0
    out = [name or RESIDUE, f'{mol.GetNumAtoms():5d}']
    for i in range(mol.GetNumAtoms()):
        out.append(f'{1:5d}{RESIDUE:<5s}{names[i]:>5s}{(i + 1) % 100000:5d}{p[i][0]:8.3f}{p[i][1]:8.3f}{p[i][2]:8.3f}')
    out.append(f'{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}')
    return '\n'.join(out + [''])


def _inpcrd(mol, xyz, name):
    flat = [f'{v:12.7f}' for v in xyz.flatten()]
    lines = [''.join(flat[k:k + 6]) for k in range(0, len(flat), 6)]
    return '\n'.join([(name or RESIDUE)[:80], f'{mol.GetNumAtoms():6d}'] + lines + [''])


def angles_and_dihedrals(mol):
    angles, dihedrals = [], []
    for a in mol.GetAtoms():
        n = sorted(x.GetIdx() for x in a.GetNeighbors())
        angles += [(n[i], a.GetIdx(), n[j]) for i in range(len(n)) for j in range(i + 1, len(n))]
    for b in mol.GetBonds():
        j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        for x in mol.GetAtomWithIdx(j).GetNeighbors():
            for y in mol.GetAtomWithIdx(k).GetNeighbors():
                i, l = x.GetIdx(), y.GetIdx()
                if i not in (k, l) and l != j:
                    dihedrals.append((i, j, k, l))
    return angles, dihedrals


def _lammps(mol, xyz, name):
    q, _ok = charges(mol)
    typed = [sybyl_type(a) for a in mol.GetAtoms()]
    kinds = sorted(set(typed))
    pt = Chem.GetPeriodicTable()
    masses = {t: pt.GetAtomicWeight(a.GetSymbol()) for t, a in zip(typed, mol.GetAtoms())}
    tid = {t: k + 1 for k, t in enumerate(kinds)}
    angles, dihedrals = angles_and_dihedrals(mol)

    def table(tuples):
        """Each term's type number: by the sorted-direction atom types it joins."""
        keys = [min(tuple(typed[i] for i in t), tuple(typed[i] for i in reversed(t))) for t in tuples]
        ids = {k: n + 1 for n, k in enumerate(sorted(set(keys)))}
        return [ids[k] for k in keys], ids
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    bt, bids = table(bonds)
    at, aids = table(angles)
    dt, dids = table(dihedrals)
    lo, hi = xyz.min(axis=0) - BOX_PADDING, xyz.max(axis=0) + BOX_PADDING
    out = [f'{name or RESIDUE}: atom types {" ".join(f"{tid[t]}={t}" for t in kinds)}; coefficients are for you to set', '',
           f'{mol.GetNumAtoms()} atoms', f'{len(bonds)} bonds', f'{len(angles)} angles', f'{len(dihedrals)} dihedrals', '0 impropers', '',
           f'{len(kinds)} atom types', f'{len(bids)} bond types', f'{len(aids)} angle types', f'{len(dids)} dihedral types', '',
           f'{lo[0]:.4f} {hi[0]:.4f} xlo xhi', f'{lo[1]:.4f} {hi[1]:.4f} ylo yhi', f'{lo[2]:.4f} {hi[2]:.4f} zlo zhi', '',
           'Masses', '']
    out += [f'{tid[t]} {masses[t]:.4f}  # {t}' for t in kinds]
    out += ['', 'Atoms  # full', '']
    out += [f'{i + 1} 1 {tid[typed[i]]} {q[i]:.4f} {xyz[i][0]:.4f} {xyz[i][1]:.4f} {xyz[i][2]:.4f}' for i in range(len(typed))]
    for title, rows, types in (('Bonds', bonds, bt), ('Angles', angles, at), ('Dihedrals', dihedrals, dt)):
        if rows:
            out += ['', title, '']
            out += [f'{n + 1} {types[n]} ' + ' '.join(str(i + 1) for i in r) for n, r in enumerate(rows)]
    return '\n'.join(out + [''])


def text(model, conf_id, fmt):
    """The file's content for one conformer (sdf_all: for all, lowest energy first)."""
    mol = model.mol
    name = model.name.split()[-1] if model.name else ''
    if fmt == 'sdf_all':
        return ''.join(_sdf(model, c) for c, _e in model.confs)
    xyz = model.positions(conf_id)
    if fmt == 'pdb':
        return _pdb(mol, xyz)
    if fmt == 'mol2':
        return _mol2(mol, xyz, name)
    if fmt == 'gro':
        return _gro(mol, xyz, name)
    if fmt == 'inpcrd':
        return _inpcrd(mol, xyz, name)
    if fmt == 'lammps':
        return _lammps(mol, xyz, name)
    if fmt == 'sdf':
        return _sdf(model, conf_id)
    if fmt == 'xyz':
        return Chem.MolToXYZBlock(mol, confId=conf_id)
    raise ValueError(f'Unknown format {fmt}')


def _sdf(model, conf_id):
    energy = dict(model.confs).get(conf_id)
    tag = '' if energy is None else f'>  <Energy_{model.method}_kcal_mol>\n{energy:.4f}\n\n'
    return Chem.MolToMolBlock(model.mol, confId=conf_id) + tag + '$$$$\n'


def write(model, conf_id, fmt, path):
    with open(path, 'w', encoding='ascii', errors='replace', newline='\n') as f:
        f.write(text(model, conf_id, fmt))
