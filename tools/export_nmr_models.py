"""
export_nmr_models.py - writes CSP5's 1H and 13C weights to nmr_models/ for nmr_engine.py, and what real CSP5 predicts to
tests/nmr_reference.json. Only needed when CSP5 changes; it needs PyTorch and CSP5 (requirements-nmr.txt), the program
itself does not.
"""
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
REFERENCE = ['CCO', 'CC(=O)O', 'c1ccccc1', 'Cc1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
             'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'CC(=O)Nc1ccc(O)cc1', 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
             'C[C@H](N)C(=O)O', 'N[C@@H](Cc1ccccc1)C(=O)O', 'C/C=C/C', 'C/C=C\\C', 'CCN(CC)CC', 'CS(C)=O',
             'ClCCl', 'FC(F)(F)c1ccccc1', 'Brc1ccccc1', 'O=C1CCCCC1', 'C1CCCCC1', 'c1ccncc1', 'c1ccsc1', 'c1ccoc1',
             'c1ccc2ccccc2c1', 'CC(C)(C)OC(=O)N', 'COc1ccccc1', 'O=[N+]([O-])c1ccccc1', 'CC#N', 'C#CC', 'C=CC=O',
             'NCCc1c[nH]c2ccccc12', 'CC(C)=CCCC(C)=CCO', 'CC1=CCC(CC1)C(C)=C', 'O=C(O)CCC(=O)O', 'CCCCCCCCCC(=O)O',
             'CP(=O)(O)O', 'C[Si](C)(C)C', 'CC(=O)[O-].[Na+]'.split('.')[0],
             'O=C([O-])[C@@H]1Cc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc2[N+]1=CC=C1C[C@@H](C(=O)O)NC(C(=O)[O-])=C1',
             'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O', 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',
             'CC(=O)OC1=CC=CC=C1', 'CCOC(=O)C', 'OCCO', 'c1ccc(cc1)-c1ccccc1', 'NC(=O)N', 'CC(O)=O']
SOLVENTS = ['cdcl3', 'dmso-d6', 'cd3od', 'd2o', 'acetone-d6', 'cd2cl2', 'c6d6', 'cd3cn', 'thf-d8', 'toluene-d8', 'dmf-d7']


def export():
    import torch
    import csp5
    from csp5 import api
    from csp5._runtime import cascade_core as core
    from csp5.model_registry import get_model_spec_by_id
    import nmr_engine
    out = ROOT / 'nmr_models'
    out.mkdir(exist_ok=True)
    package = Path(csp5.__file__).parent
    for name in nmr_engine.MODELS.values():
        spec = next(s for s in (get_model_spec_by_id(m) for m in _ids()) if s.model_name == name)
        artifacts = api._load_artifacts_cached(spec.model_id, 'cpu')
        state = torch.load(package / 'models' / name / 'best_model.pt', map_location='cpu', weights_only=False)
        state = core._unwrap_state_dict(state)
        np.savez_compressed(out / (name + '.npz'), **{k: v.cpu().numpy() for k, v in state.items()})
        meta = json.loads((package / 'models' / name / 'model_metadata.json').read_text(encoding='utf-8'))
        assert meta['quantile_parameterization'] == 'direct' and meta['output_scale'] == 1.0 and meta['output_bias'] == 0.0
        tokens = {}
        for z in range(1, 119):
            atom = Chem_atom(z)
            tokens[z] = int(artifacts.preprocessor.atom_tokenizer(artifacts.preprocessor.atom_features(atom)))
        config = artifacts.model.config
        (out / (name + '.json')).write_text(json.dumps({
            'tokens': tokens, 'solvents': meta['solvent_vocab'], 'calibration': meta.get('calibration'),
            'n_neighbors': spec.graph_n_neighbors, 'cutoff': config.cutoff,
            'envelope_exponent': config.envelope_exponent, 'graph_cutoff': spec.graph_cutoff_angstrom}), encoding='utf-8')
    shutil.copy(next(Path(csp5.__file__).parent.parent.glob('csp5-*.dist-info/licenses/LICENSE')), out / 'CSP5_LICENSE.txt')


def Chem_atom(z):
    from rdkit import Chem
    return Chem.Atom(z)


def _ids():
    from csp5 import model_registry
    return list(model_registry._MODEL_SPECS_BY_ID)


def reference(path=ROOT / 'tests' / 'nmr_reference.json', smiles=REFERENCE):
    import csp5
    rows = []
    for k, smi in enumerate(smiles):
        for nucleus in ('1H', '13C'):
            solvent = SOLVENTS[k % len(SOLVENTS)]
            model = 'CSP5q-13C' if nucleus == '13C' else 'CSP5q-1H'
            r = csp5.predict_smiles([smi], nucleus=nucleus, model_name=model, solvent=solvent, random_seed=1)
            if r.failures or not r.molecule_records:
                rows.append(dict(smiles=smi, nucleus=nucleus, solvent=solvent, failed=True))
                continue
            p = r.predictions
            rows.append(dict(smiles=smi, nucleus=nucleus, solvent=solvent, atom_index=[int(x) for x in p['atom_index']],
                             shift=[float(x) for x in p['shift_ppm']], std=[float(x) for x in p['shift_std_ppm']]))
    Path(path).write_text(json.dumps(rows), encoding='utf-8')
    return rows


if __name__ == '__main__':
    export()
    print(len(reference()), 'reference predictions')
