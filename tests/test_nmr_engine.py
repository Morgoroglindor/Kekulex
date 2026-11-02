"""The numpy engine against what real CSP5 predicted (tests/nmr_reference.json, from tools/export_nmr_models.py)."""
import json
import os

import numpy as np
import pytest
from rdkit import Chem

import nmr_engine

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'nmr_reference.json'), encoding='utf-8') as f:
    REFERENCE = json.load(f)


def test_reference_is_broad():
    assert len(REFERENCE) >= 90
    assert len({r['solvent'] for r in REFERENCE}) >= 8


@pytest.mark.parametrize('row', [r for r in REFERENCE if not r.get('failed')],
                         ids=lambda r: f"{r['nucleus']}-{r['smiles'][:20]}")
def test_matches_csp5(row):
    p = nmr_engine.predict(row['smiles'], row['nucleus'], row['solvent'])
    assert p.atom_index == row['atom_index']
    assert np.abs(np.array(p.shift_ppm) - row['shift']).max() <= 0.011       # both are rounded to 0.01 ppm
    assert np.abs(np.array(p.shift_std_ppm) - row['std']).max() <= 0.011


def test_no_pytorch_needed():
    source = open(os.path.join(os.path.dirname(HERE), 'nmr_engine.py'), encoding='utf-8').read()
    assert 'import torch' not in source and 'import pandas' not in source


def test_readable_errors():
    with pytest.raises(ValueError):
        nmr_engine.predict('not a molecule')
    with pytest.raises(ValueError, match='no 13C'):
        nmr_engine.predict('O', '13C')


def test_unknown_solvent_is_still_predicted():
    assert nmr_engine.predict('CCO', '1H', 'no such solvent').shift_ppm
