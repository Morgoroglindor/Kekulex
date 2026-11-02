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
