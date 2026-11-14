"""Radicals, carbenes, radical ions and typed hydrogen counts (MS fragments, reactive intermediates)."""
import pytest
from PyQt6.QtTest import QTest

import cdxml
import chem
import model
import render
from model import Doc


@pytest.fixture
def win(qapp):
    import kekulex as se
    w = se.MainWindow()
    w.resize(1300, 850)
    w.show()
    qapp.processEvents()
    yield w
    w.dirty = False
    w.close()


def _doc(smiles):
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles), (3, 3)); return d


def _el(d, el):
    return next(a for a, at in d.atoms.items() if at.el == el)


@pytest.mark.parametrize('smiles, formula, radicals', [
    ('[CH3]', 'CH3', [1]),                                  # methyl radical
    ('[CH2]', 'CH2', [2]),                                  # carbene
    ('C[CH2]', 'C2H5', [0, 1]),                             # ethyl radical
    ('C[O]', 'CH3O', [0, 1]),                               # methoxyl
    ('CC#[O+]', 'C2H3O', [0, 0, 0]),                        # acylium (even-electron MS fragment)
    ('C[CH2+]', 'C2H5', [0, 0]),                            # ethyl cation, six electrons
])
def test_radicals_and_cations_from_smiles(smiles, formula, radicals):
    d = _doc(smiles)
    info = chem.analyse(d)
    assert info.formula == formula and not info.errors
    assert sorted(at.radical for at in d.atoms.values()) == sorted(radicals)


def test_radical_takes_one_hydrogen_away():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.atoms[a].radical = 1
    assert d.implicit_h(a) == 2 and not d.valence_problem(a)
    assert chem.analyse(d).formula == 'C2H5'
    d.atoms[a].radical = 2                                  # a carbene carbon
    assert chem.analyse(d).formula == 'C2H4'


def test_radical_cation_is_an_odd_electron_ion():
    d = _doc('CC(C)=O')                                     # acetone M•+ (EI)
    o = _el(d, 'O')
    d.atoms[o].charge, d.atoms[o].radical = 1, 1
    info = chem.analyse(d)
    assert info.formula == 'C3H6O' and not info.errors
    assert info.ions[0].label == '[M]•+'
    assert info.ions[0].mz == pytest.approx(58.0413, abs=1e-3)   # 58.0419 - one electron
    assert '•+' in info.formula_html


def test_rdkit_sees_the_radical():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.atoms[a].radical = 1
    full, _grown = d.expanded()
    mol, _order, err = chem.to_mol(full)
    assert not err
    assert sum(at.GetNumRadicalElectrons() for at in mol.GetAtoms()) == 1
    assert sum(at.GetTotalNumHs() for at in mol.GetAtoms()) == 5


def test_molblock_writes_rad_lines():
    d = _doc('C[CH2]')
    block, _order = d.to_molblock()
    assert 'M  RAD  1' in block


@pytest.mark.parametrize('typed, el, charge, hs, radical', [
    ('CH2', 'C', 0, 2, 0), ('CH3•', 'C', 0, 3, 1), ('C••', 'C', 0, None, 2), ('O•+', 'O', 1, None, 1),
    ('ch2.', 'C', 0, 2, 1), ('N•', 'N', 0, None, 1), ('CH3+', 'C', 1, 3, 0)])
def test_typed_labels(typed, el, charge, hs, radical):
    assert model.parse_label(typed) == ('element', el, charge)
    assert model.label_extras(typed) == (hs, radical)


def test_typing_ch2_on_a_carbon_keeps_two_hydrogens(win):
    c = win.canvas
    c.set_doc(_doc('CCC'))
    a = next(a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 1)
    c.edit_label(a); c.editor[0].setText('CH2'); c.commit_editor()
    at = c.doc.atoms[a]
    assert at.hcount == 2 and at.radical == 0                # a hypovalent carbon, as drawn
    assert chem.analyse(c.doc).formula == 'C3H7'
    c.edit_label(a); assert c.editor[0].text() == 'CH2'; c.editor[0].setText('CH2•'); c.commit_editor()
    assert at.radical == 1                                   # 2 H is the radical's own number: not fixed
    assert at.hcount is None
    assert chem.analyse(c.doc).formula == 'C3H7'
    c.edit_label(a); c.editor[0].setText('CH3'); c.commit_editor()
    assert at.radical == 0 and at.hcount is None and chem.analyse(c.doc).formula == 'C3H8'


def test_dot_key_cycles_the_radical(win):
    c = win.canvas
    c.set_doc(_doc('CC'))
    a = next(iter(c.doc.atoms))
    c.hover = ('atom', a)
    for want in (1, 2, 0):
        QTest.keyClick(c, '.')
        assert c.doc.atoms[a].radical == want


def test_radical_dot_is_drawn_beside_the_atom():
    d = _doc('C[CH2]')
    a = next(a for a, at in d.atoms.items() if at.radical)
    lab = render.Renderer(d).labels[a]
    texts = [p[0] for p in lab['pieces']]
    assert texts[0] == 'C' and 'H' in texts and '•' in texts
    o = _doc('CC(C)=O')
    ox = _el(o, 'O')
    o.atoms[ox].charge, o.atoms[ox].radical = 1, 1
    assert '•+' in [p[0] for p in render.Renderer(o).labels[ox]['pieces']]


def test_radicals_round_trip_through_files(tmp_path):
    d = _doc('C[CH2]')
    a = next(a for a, at in d.atoms.items() if at.radical)
    back = Doc.from_dict(d.to_dict())
    assert back.atoms[a].radical == 1
    path = str(tmp_path / 'r.cdxml')
    cdxml.write(d, path)
    assert 'Radical="Doublet"' in open(path, encoding='utf-8').read()
    back, _notes = cdxml.read(path)
    assert sorted(at.radical for at in back.atoms.values()) == [0, 1]
    assert chem.analyse(back).formula == 'C2H5'
