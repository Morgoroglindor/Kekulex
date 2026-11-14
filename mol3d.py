"""
mol3d.py - the chemistry behind the 3D tab: a drawn structure to 3D conformers (RDKit ETKDG + MMFF94, UFF when MMFF has
no parameters) and the measurements on them. No Qt here.
"""
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms

ENERGY_SAME = 0.02                                          # kcal/mol: closer than this counts as the same conformer
LARGE = 150                                                 # heavy atoms: fewer conformers and a warning
EMBED_TIMEOUT = 60                                          # seconds for the search of one molecule


@dataclass
class Model3D:
    mol: object                                             # RDKit molecule with explicit H and the conformers
    source: object = None                                   # the molecule as sent, to search conformers again
    name: str = ''
    method: str = ''                                        # 'MMFF94', 'UFF' or '' (not minimised)
    confs: list = field(default_factory=list)               # [(conformer id, energy in kcal/mol or None)], lowest first
    notes: list = field(default_factory=list)               # warnings to show with the model

    def positions(self, conf_id):
        return self.mol.GetConformer(conf_id).GetPositions()


FORCE_FIELDS = ('MMFF94', 'MMFF94s', 'UFF')


def _force_field(mol, wanted=None):
    """The force field to use: wanted if it has parameters for mol, else the best that does ('' for none)."""
    if wanted == 'UFF' and AllChem.UFFHasAllMoleculeParams(mol):
        return 'UFF'
    if AllChem.MMFFHasAllMoleculeParams(mol):
        return wanted if wanted in ('MMFF94', 'MMFF94s') else 'MMFF94'
    if AllChem.UFFHasAllMoleculeParams(mol):
        return 'UFF'
    return ''


def _optimise(mol, method, max_iters=2000):
    """[(converged, energy)] per conformer, in the order of mol's conformers."""
    if method in ('MMFF94', 'MMFF94s'):
        return AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=max_iters, numThreads=0, mmffVariant=method)
    if method == 'UFF':
        return AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=max_iters, numThreads=0)
    return [(1, None)] * mol.GetNumConformers()


def _unspecified_centres(mol):
    return sum(1 for _a, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
               if tag == '?')


def build(mol, name='', n_confs=20, seed=0xF00D, force_field=None):
    """A Model3D of the molecule mol (keeps its stereo). Raises ValueError with a readable text when there is no 3D
    geometry to be had."""
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError('There is no structure to build.')
    source = mol
    heavy = mol.GetNumHeavyAtoms()
    if mol.GetNumConformers():                              # a drawing: double-bond geometry from its 2D coordinates
        mol = Chem.Mol(mol)
        Chem.DetectBondStereoChemistry(mol, mol.GetConformer())
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    notes = []
    unknown = _unspecified_centres(mol)
    if unknown:
        notes.append(f'{unknown} stereocentre{"s" if unknown > 1 else ""} not specified in the drawing: '
                     'the 3D model shows one arbitrary choice.')
    m = Chem.AddHs(mol)
    m.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 0
    params.timeout = EMBED_TIMEOUT
    count = n_confs if heavy <= 50 else 5 if heavy <= LARGE else 1
    ids = list(AllChem.EmbedMultipleConfs(m, count, params))
    if not ids:                                             # strained or large: start from random coordinates
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(m, count, params))
    if not ids:
        raise ValueError('RDKit could not find a 3D geometry for this structure (too strained, or too large).')
    method = _force_field(m, force_field)
    if method == 'UFF' and force_field != 'UFF':
        notes.append('MMFF94 has no parameters for this structure: minimised with UFF, a rougher force field.')
    elif not method:
        notes.append('No force field has parameters for this structure (metals?): the geometry is not minimised '
                     'and bond lengths may be wrong.')
    if heavy > LARGE:
        notes.append('A large molecule: this is one rough shape, not a realistic fold.')
    results = _optimise(m, method)
    if method and any(r[0] != 0 for r in results):
        notes.append('Some conformers did not fully converge: use Minimize.')
    ranked = sorted(zip(ids, results), key=lambda c: (c[1][1] is None, c[1][1] or 0))
    confs = []
    for cid, (_conv, e) in ranked:
        if e is None or not confs or all(abs(e - other) > ENERGY_SAME for _c, other in confs):
            confs.append((cid, e))
    try:
        rdMolAlign.AlignMolConformers(m, confIds=[c for c, _e in confs])    # switching conformers keeps the view steady
    except RuntimeError:
        pass
    return Model3D(m, source, name, method, confs, notes)


def minimize(model, conf_id, max_iters=5000):
    """Minimises one conformer further; returns (its new energy or None, converged). The energies stay sorted."""
    if not model.method:
        return None, False
    ff = _field(model, conf_id)
    ff.Initialize()
    converged = ff.Minimize(maxIts=max_iters) == 0
    e = ff.CalcEnergy()
    model.confs = sorted(((c, e if c == conf_id else old) for c, old in model.confs), key=lambda c: c[1])
    return e, converged


def _field(model, conf_id):
    if model.method == 'UFF':
        return AllChem.UFFGetMoleculeForceField(model.mol, confId=conf_id)
    props = AllChem.MMFFGetMoleculeProperties(model.mol, mmffVariant=model.method)
    return AllChem.MMFFGetMoleculeForceField(model.mol, props, confId=conf_id)


def minimize_all(model, max_iters=5000):
    """Minimises every conformer further; returns how many did not converge."""
    return sum(not minimize(model, c, max_iters)[1] for c, _e in list(model.confs))


def rms_gradient(model, conf_id):
    """How far the conformer is from a minimum: RMS force, kcal/mol/Å (None without a force field)."""
    if not model.method:
        return None
    ff = _field(model, conf_id)
    ff.Initialize()
    g = np.array(ff.CalcGrad())
    return float(np.sqrt((g * g).sum() / max(len(g), 1)))


def measure(model, conf_id, atoms):
    """(what, value, unit) for 2, 3 or 4 atoms: distance, angle or dihedral; None otherwise."""
    conf = model.mol.GetConformer(conf_id)
    if len(atoms) == 2:
        return 'Distance', rdMolTransforms.GetBondLength(conf, *atoms), 'Å'
    if len(atoms) == 3:
        return 'Angle', rdMolTransforms.GetAngleDeg(conf, *atoms), '°'
    if len(atoms) == 4:
        return 'Dihedral', rdMolTransforms.GetDihedralDeg(conf, *atoms), '°'
    return None


def atom_label(mol, i):
    return f'{mol.GetAtomWithIdx(i).GetSymbol()}{i + 1}'


# ----------------------------------------------------------------------------------------------- surfaces and shape
PROBE = 1.4                                                 # Å: a water molecule, for the solvent-accessible surface
SURFACES = (('sas', 'Solvent accessible'), ('vdw', 'Van der Waals'))


def vdw_radii(mol):
    table = Chem.GetPeriodicTable()
    return np.array([table.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])


def _sphere(n):
    """n points spread evenly over the unit sphere (golden spiral)."""
    k = np.arange(n) + 0.5
    z = 1 - 2 * k / n
    r = np.sqrt(1 - z * z)
    phi = np.pi * (1 + 5 ** 0.5) * k
    return np.stack((r * np.cos(phi), r * np.sin(phi), z), axis=1)


def surface_dots(xyz, radii, n=240):
    """Shrake-Rupley: (points, owning atom, area in Å²) of the outer surface of the union of the spheres."""
    unit = _sphere(n)
    points, owners, area = [], [], 0.0
    reach = radii.max() * 2
    for i, (c, r) in enumerate(zip(xyz, radii)):
        near = np.where((np.linalg.norm(xyz - c, axis=1) < r + reach) & (np.arange(len(xyz)) != i))[0]
        pts = c + unit * r
        free = np.ones(n, dtype=bool)
        for j in near:
            free &= np.linalg.norm(pts - xyz[j], axis=1) >= radii[j]
        points.append(pts[free])
        owners.append(np.full(int(free.sum()), i))
        area += 4 * np.pi * r * r * free.sum() / n
    return np.concatenate(points), np.concatenate(owners), float(area)


def surface_radii(model, kind):
    radii = vdw_radii(model.mol)
    return radii + PROBE if kind == 'sas' else radii


def atom_values(model, scheme):
    """A number per atom for colouring a surface: Gasteiger partial charge, or Crippen's logP contribution
    (hydrophobicity). None for schemes that take the atom's own colour."""
    mol = model.mol
    if scheme == 'charge':
        m = Chem.Mol(mol)
        AllChem.ComputeGasteigerCharges(m)
        vals = [m.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(m.GetNumAtoms())]
        return np.nan_to_num(np.array(vals, dtype=float))
    if scheme == 'hydrophobicity':
        from rdkit.Chem import rdMolDescriptors
        return np.array([c[0] for c in rdMolDescriptors._CalcCrippenContribs(mol)], dtype=float)
    return None


def _donor_hydrogens(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1 and any(n.GetAtomicNum() in (7, 8) for n in a.GetNeighbors())]


def _acceptors(mol):
    """O always; N with a lone pair (not N+ with four bonds, not an amide-like or aromatic NH)."""
    out = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 8:
            out.append(a.GetIdx())
        elif a.GetAtomicNum() == 7 and a.GetFormalCharge() <= 0 and a.GetDegree() < 4 and not (a.GetIsAromatic() and a.GetTotalNumHs()):
            out.append(a.GetIdx())
    return out


def hydrogen_bonds(model, conf_id, max_dist=2.5, min_angle=120.0):
    """[(hydrogen, acceptor)] in the conformer: H...A at most max_dist Å, D-H...A at least min_angle degrees, and at
    least four bonds apart (no 1,2-diol style near-contacts)."""
    mol = model.mol
    xyz = model.positions(conf_id)
    paths = Chem.GetDistanceMatrix(mol)
    found = []
    for h in _donor_hydrogens(mol):
        d = mol.GetAtomWithIdx(h).GetNeighbors()[0].GetIdx()
        for a in _acceptors(mol):
            if a == d or paths[h][a] < 4:
                continue
            ha = xyz[a] - xyz[h]
            dist = float(np.linalg.norm(ha))
            if dist > max_dist:
                continue
            hd = xyz[d] - xyz[h]
            cos = float(hd @ ha) / (np.linalg.norm(hd) * dist)
            if np.degrees(np.arccos(np.clip(cos, -1, 1))) >= min_angle:
                found.append((h, a))
    return found


def shape_properties(model, conf_id):
    """[(name, text)] of the conformer's surfaces, volume and shape, and its hydrogen-bond counts."""
    from rdkit.Chem import Descriptors3D, rdMolDescriptors
    mol = model.mol
    xyz = model.positions(conf_id)
    _p, _o, sas = surface_dots(xyz, surface_radii(model, 'sas'))
    _p, _o, vdw_area = surface_dots(xyz, surface_radii(model, 'vdw'))
    volume = AllChem.ComputeMolVolume(mol, confId=conf_id, gridSpacing=0.25)
    ovality = vdw_area / (4 * np.pi * (3 * volume / (4 * np.pi)) ** (2 / 3)) if volume else float('nan')
    pmi = [getattr(Descriptors3D, f'PMI{k}')(mol, confId=conf_id) for k in (1, 2, 3)]
    rows = [('Solvent-accessible area', f'{sas:.0f} Å² (probe {PROBE} Å)'),
            ('Van der Waals area', f'{vdw_area:.0f} Å²'),
            ('Van der Waals volume', f'{volume:.0f} Å³'),
            ('Ovality', f'{ovality:.3f} (1 = a sphere)'),
            ('Principal moments', ', '.join(f'{v:.0f}' for v in pmi) + ' amu·Å²'),
            ('Radius of gyration', f'{Descriptors3D.RadiusOfGyration(mol, confId=conf_id):.2f} Å'),
            ('H-bond donors / acceptors', f'{rdMolDescriptors.CalcNumHBD(mol)} / {rdMolDescriptors.CalcNumHBA(mol)}'),
            ('H-bonds inside this conformer', str(len(hydrogen_bonds(model, conf_id))))]
    return rows
