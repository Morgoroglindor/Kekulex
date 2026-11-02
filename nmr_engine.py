"""
nmr_engine.py - runs the CSP5 NMR models (Lin et al., MIT licence, CSP5 0.3.1) with numpy alone, so the program does not
carry PyTorch. The network (PaiNN), the 3D geometry, the graph and the quantile calibration follow CSP5's release code
step by step; tools/export_nmr_models.py wrote its weights to nmr_models/ and tests/nmr_reference.json holds what real
CSP5 predicts, which tests/test_nmr_engine.py compares with.
"""
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Lipinski, rdMolDescriptors

Q90_ZSCORE = 1.2815515655446004
SHIFT_WINDOW = {'13C': (-50.0, 300.0), '1H': (-5.0, 20.0)}
ATOMIC_NUMBER = {'13C': 6, '1H': 1}
MODELS = {'13C': 'CSP5q-13C', '1H': 'CSP5q-1H'}
MIN_CONTACT = 0.90                                          # Å: closer atoms mean a broken geometry
EMBED_TRIES = 20
EMBED_TIMEOUT = 5                                           # s per attempt
FF_ITERS = 200
_cache = {}


@dataclass
class Prediction:
    mol: object                                             # the molecule with H and its conformer
    atom_index: list = field(default_factory=list)
    shift_ppm: list = field(default_factory=list)
    shift_std_ppm: list = field(default_factory=list)


def models_dir():
    base = getattr(sys, '_MEIPASS', None) or os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, 'nmr_models')


def available():
    return all(os.path.exists(os.path.join(models_dir(), name + '.npz')) for name in MODELS.values())


class _Model:
    def __init__(self, name):
        with np.load(os.path.join(models_dir(), name + '.npz')) as z:
            self.w = {k: z[k] for k in z.files}
        with open(os.path.join(models_dir(), name + '.json'), encoding='utf-8') as f:
            self.meta = json.load(f)
        self.tokens = {int(k): v for k, v in self.meta['tokens'].items()}
        self.depth = 1 + max(int(k.split('.')[1]) for k in self.w if k.startswith('convs.'))

    def linear(self, name, x):
        out = x @ self.w[name + '.weight'].T
        bias = self.w.get(name + '.bias')
        return out if bias is None else out + bias

    def forward(self, tokens, pos, edges, targets, solvent_id):
        """The network's 99 raw quantiles of each target atom."""
        src, dst = edges[:, 0], edges[:, 1]
        n = len(tokens)
        z = self.w['embedding.weight'][tokens]
        v = np.full((n, 3, z.shape[1]), np.finfo(np.float32).eps, dtype=np.float32)
        rij = pos[src] - pos[dst]
        dist = np.linalg.norm(rij, axis=-1, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            inv = np.where(dist > 0, 1.0 / dist, np.zeros_like(dist))
        rij = rij * inv
        rbf = self._bessel(dist)
        for k in range(self.depth):
            ds, dv = self._conv(k, z, v, rbf, rij, src, dst, n)
            z, v = z + ds, v + dv
            ds, dv = self._update(k, z, v)
            z, v = z + ds, v + dv
        pooled = z[targets]
        pred = self.linear('mlp.2', _silu(self.linear('mlp.0', pooled)))
        pred = pred + self.w['solvent_bias.weight'][solvent_id]
        emb = np.broadcast_to(self.w['solvent_embedding.weight'][solvent_id], (len(targets), self.w['solvent_embedding.weight'].shape[1]))
        hidden = _silu(self.linear('solvent_adapter.0', np.concatenate([pooled, emb], axis=-1)))
        return pred + self.linear('solvent_adapter.3', hidden)

    def _bessel(self, dist):
        cutoff, p = self.meta['cutoff'], self.meta['envelope_exponent'] + 1
        a, b, c = -(p + 1) * (p + 2) / 2, p * (p + 2), -p * (p + 1) / 2
        x = dist * np.float32(1.0 / cutoff)
        with np.errstate(divide='ignore', invalid='ignore'):
            env = 1.0 / x + a * x ** (p - 1) + b * x ** p + c * x ** (p + 1)
        env = np.where(x < 1, env, np.zeros_like(x))
        return env * np.sin(self.w['bessel_basis.frequencies'] * x)

    def _conv(self, k, z, v, rbf, rij, src, dst, n):
        s = self.linear(f'convs.{k}.dense_phi', _silu(self.linear(f'convs.{k}.dense_s', z)))[dst]
        sw = s * self.linear(f'convs.{k}.dense_w', rbf)
        sw1, sw2, sw3 = np.split(sw, 3, axis=-1)
        ds = np.zeros_like(z)
        np.add.at(ds, src, sw1)
        dv_edge = v[dst] * sw2[:, None, :] + rij[:, :, None] * sw3[:, None, :]
        dv = np.zeros_like(v)
        np.add.at(dv, src, dv_edge)
        return ds, dv

    def _update(self, k, z, v):
        v_v = self.linear(f'updates.{k}.lin_v', v)
        v_u = self.linear(f'updates.{k}.lin_u', v)
        v_prod = (v_u * v_v).sum(axis=1)
        v_norm = np.linalg.norm(v_v, axis=1)
        a = self.linear(f'updates.{k}.dense_a',
                        _silu(self.linear(f'updates.{k}.dense1', np.concatenate([z, v_norm], axis=-1))))
        a_vv, a_sv, a_ss = np.split(a, 3, axis=-1)
        return v_prod * a_sv + a_ss, v_u * a_vv[:, None, :]


def _silu(x):
    with np.errstate(over='ignore'):                        # a very negative x gives -0, as it should
        return x / (1.0 + np.exp(-x))


def _model(nucleus):
    name = MODELS[nucleus]
    if name not in _cache:
        _cache[name] = _Model(name)
    return _cache[name]


# ----------------------------------------------------------------------------------------------- geometry
def _geometry_key(mol):
    m = Chem.RemoveHs(Chem.Mol(mol))
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)


def _flat_polycyclic(mol):
    """Large fused aromatics make ETKDG hang: 2D coordinates and a force field do instead."""
    heavy = mol.GetNumHeavyAtoms()
    if not heavy or rdMolDescriptors.CalcNumRings(mol) < 10 or Lipinski.NumRotatableBonds(mol):
        return False
    return sum(a.GetIsAromatic() for a in mol.GetAtoms()) / heavy >= 0.9


def _embed(smiles, seed):
    """The molecule with H and CSP5's one geometry: embedded in canonical atom order, MMFF (else UFF) minimised."""
    start = Chem.MolFromSmiles(smiles)
    if start is None:
        raise ValueError('The structure could not be read.')
    mol = Chem.AddHs(start)
    mol.RemoveAllConformers()
    key = _geometry_key(mol)
    canon = Chem.AddHs(Chem.MolFromSmiles(key))
    match = mol.GetSubstructMatch(canon, useChirality=True)
    if len(match) != canon.GetNumAtoms() or canon.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError('No 3D geometry could be made for this structure.')
    done = False
    if _flat_polycyclic(Chem.RemoveHs(Chem.Mol(canon))):
        import random
        AllChem.Compute2DCoords(canon)
        base = canon.GetConformer()
        rng = random.Random(seed)
        coords = [base.GetAtomPosition(i) for i in range(canon.GetNumAtoms())]
        canon.RemoveAllConformers()
        conf = Chem.Conformer(canon.GetNumAtoms())
        for i, p in enumerate(coords):
            conf.SetAtomPosition(i, (float(p.x), float(p.y), (rng.random() - 0.5) * 0.1))
        canon.AddConformer(conf, assignId=True)
        done = True
    else:
        params = AllChem.ETKDGv3()
        params.timeout = EMBED_TIMEOUT
        for attempt in range(EMBED_TRIES):
            canon.RemoveAllConformers()
            params.randomSeed = int((seed + attempt) & 0x7FFFFFFF)
            params.useRandomCoords = attempt > 0
            if AllChem.EmbedMolecule(canon, params=params) == 0:
                done = True
                break
    if not done:
        raise ValueError('No 3D geometry could be made for this structure.')
    status = -1
    if AllChem.MMFFHasAllMoleculeParams(canon):
        status = AllChem.MMFFOptimizeMolecule(canon, maxIters=FF_ITERS)
    if status not in (0, 1) and AllChem.UFFHasAllMoleculeParams(canon):
        status = AllChem.UFFOptimizeMolecule(canon, maxIters=FF_ITERS)
    if status not in (0, 1):
        raise ValueError('The 3D geometry of this structure could not be minimised.')
    source = canon.GetConformer()
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for c, i in enumerate(match):
        conf.SetAtomPosition(int(i), source.GetAtomPosition(c))
    mol.AddConformer(conf, assignId=True)
    return mol


def _edges(pos, n_neighbors, cutoff):
    n = len(pos)
    if n == 1:
        return np.array([[0, 0]], dtype=np.int64)
    d = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
    out = []
    for s in range(n):
        nb = []
        for t in np.argsort(d[s])[1:]:
            if d[s, t] >= cutoff:
                break
            nb.append(int(t))
            if len(nb) >= n_neighbors:
                break
        out += [(s, t) for t in (nb or [s])]
    return np.array(out, dtype=np.int64)


def _check_geometry(pos):
    if not np.isfinite(pos).all():
        raise ValueError('The 3D geometry of this structure is broken.')
    d = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(d, np.inf)
    if d.min() < MIN_CONTACT:
        raise ValueError('The 3D geometry of this structure has atoms on top of each other.')


# ----------------------------------------------------------------------------------------------- quantiles
def _calibrate(q, cal):
    q = np.sort(q.astype(np.float64), axis=1)
    if not cal:
        return q
    median = q[:, 49:50]
    q = median + float(cal.get('median_scale', 1.0)) * (q - median)
    if cal.get('level_map') is not None:
        pos = (np.clip(np.asarray(cal['level_map'], dtype=np.float64), 0.01, 0.99) - 0.01) * 100.0
        lower = np.clip(np.floor(pos).astype(np.int64), 0, 97)
        weight = np.clip(pos - lower, 0.0, 1.0)
        q = q[:, lower] * (1.0 - weight) + q[:, lower + 1] * weight
    tail = cal.get('tail')
    if tail:
        q05, q50, q95 = q[:, 4].copy(), q[:, 49].copy(), q[:, 94].copy()
        for i, c in enumerate(tail['lower_multipliers']):
            q[:, i] = np.minimum(q[:, i], q05 - float(c) * (q50 - q05))
        for i, c in zip(range(95, 99), tail['upper_multipliers']):
            q[:, i] = np.maximum(q[:, i], q95 + float(c) * (q95 - q50))
    return q


# ----------------------------------------------------------------------------------------------- prediction
def predict(smiles, nucleus='1H', solvent='cdcl3', seed=1):
    """The predicted shifts of one molecule (SMILES). Raises ValueError with a readable text."""
    RDLogger.DisableLog('rdApp.*')
    model = _model(nucleus)
    meta = model.meta
    if meta['solvents'].count(solvent) == 0:
        solvent = 'unknown'
    mol = _embed(smiles, seed)
    atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
    targets = np.array([i for i, z in enumerate(atoms) if z == ATOMIC_NUMBER[nucleus]], dtype=np.int64)
    if not targets.size:
        raise ValueError(f'This structure has no {nucleus} atoms.')
    missing = sorted({z for z in atoms if z not in model.tokens})
    if missing:
        raise ValueError('CSP5 cannot use element number ' + ', '.join(map(str, missing)) + '.')
    pos = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float32)
    _check_geometry(pos.astype(np.float64))
    tokens = np.array([model.tokens[z] for z in atoms], dtype=np.int64)
    edges = _edges(pos, meta['n_neighbors'], meta['cutoff'])
    raw = model.forward(tokens, pos, edges, targets, meta['solvents'].index(solvent))
    q = _calibrate(raw, meta.get('calibration'))
    shift = np.round(q[:, 49], 2)
    low, high = SHIFT_WINDOW[nucleus]
    if not np.isfinite(q).all() or (shift < low).any() or (shift > high).any():
        raise ValueError('CSP5 could not predict this structure (shifts out of range).')
    std = np.round(np.maximum(np.abs(q[:, 89] - q[:, 9]) / (2.0 * Q90_ZSCORE), 1.0e-4), 4)
    return Prediction(mol, [int(i) for i in targets], [float(x) for x in shift], [float(x) for x in std])
