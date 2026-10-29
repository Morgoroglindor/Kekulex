"""Formula, masses, adduct m/z, SMILES round trips (model <-> RDKit)."""
import math

import pytest

import chem
from model import Doc


def doc_of(smiles):
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text(smiles))
    return doc


def canon(smiles):
    from rdkit import Chem
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def ion(info, label):
    return {i.label: i.mz for i in info.ions}[label]


def test_caffeine_masses():
    info = chem.analyse(doc_of('Cn1cnc2c1c(=O)n(C)c(=O)n2C'))
    assert info.formula == 'C8H10N4O2'
    assert info.exact == pytest.approx(194.080376, abs=2e-6)
    assert info.average == pytest.approx(194.19, abs=0.01)
    assert ion(info, '[M+H]+') == pytest.approx(195.087652, abs=2e-6)
    assert ion(info, '[M+Na]+') == pytest.approx(217.069596, abs=2e-6)
    assert ion(info, '[M-H]−') == pytest.approx(193.073100, abs=2e-6)
    assert ion(info, '[M+2H]2+') == pytest.approx(98.047464, abs=2e-6)


def test_betalamic_acid():
    info = chem.analyse(doc_of('OC(=O)C1CC(=CC=O)C=C(N1)C(=O)O'))
    assert info.formula == 'C9H9NO5'
    assert ion(info, '[M+H]+') == pytest.approx(212.05535, abs=1e-5)


def test_charged_structure_is_its_own_ion():
    info = chem.analyse(doc_of('C[N+](C)(C)C'))           # tetramethylammonium
    assert info.charge == 1 and info.formula == 'C4H12N'
    assert ion(info, '[M]+') == pytest.approx(74.096426, abs=2e-6)
    assert ion(info, '[M+H]2+') == pytest.approx((74.096426 + 1.007276) / 2, abs=2e-6)
    assert not any('−' in i.label for i in info.ions)


@pytest.mark.parametrize('smiles', ['C[C@@H](N)C(=O)O', 'C/C=C/C', 'C/C=C\\C', 'c1ccc2ccccc2c1', 'O=C([O-])CC[NH3+]',
                                    'N#Cc1ccncc1', 'C[C@H]1CC[C@@H](O)CC1', 'O=S(=O)(O)c1ccccc1', 'C1CC1Br'])
def test_smiles_round_trip(smiles):
    info = chem.analyse(doc_of(smiles))
    assert not info.errors
    assert info.smiles == canon(smiles)


def test_drawn_atoms_get_hydrogens():
    doc = Doc()
    c = doc.add_atom(0, 0); o = doc.add_atom(1, 0, 'O'); n = doc.add_atom(-1, 0, 'N')
    doc.add_bond(c, o); doc.add_bond(c, n)
    assert doc.implicit_h(c) == 2 and doc.implicit_h(o) == 1 and doc.implicit_h(n) == 2
    assert chem.analyse(doc).formula == 'CH5NO'
    doc.atoms[n].charge = 1
    assert doc.implicit_h(n) == 3


def test_valence_error_reported():
    doc = Doc()
    c = doc.add_atom(0, 0)
    for k in range(5):
        doc.add_bond(c, doc.add_atom(math.cos(k), math.sin(k)))
    assert doc.valence_error(c)
    assert chem.analyse(doc).errors


def test_mol_from_text():
    assert chem.mol_from_text('CCO') is not None
    assert chem.mol_from_text('not a smiles at all') is None
    block, _ = doc_of('CCO').to_molblock()
    assert chem.mol_from_text(block).GetNumAtoms() == 3


def test_hill_order():
    assert chem.analyse(doc_of('ClCCBr')).formula == 'C2H4BrCl'
    assert chem.analyse(doc_of('O')).formula == 'H2O'
