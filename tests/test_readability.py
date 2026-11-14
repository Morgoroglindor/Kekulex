"""Nothing drawn sits on anything else: labels, H, charges, radical dots, lone pairs, atom numbers, R/S, map
numbers and captions keep clear of each other and of the bonds (overlaps.find); plus pictures on the drawing."""
import json
import random

import pytest
from PyQt6 import QtCore, QtGui, QtWidgets

import chem
import model
import ms
import overlaps
import reactions
import render
from model import Doc

SMILES = ['CCO', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
          'C[C@H](N)C(=O)O', 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O', 'FC(F)(F)c1ccc(Br)cc1Cl', 'CS(C)(=O)=O',
          'C[NH3+]', 'O=[N+]([O-])c1ccccc1', 'P(=O)(O)(O)O', 'OP(=O)(O)OP(=O)(O)O',
          'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'NC(=O)c1cccnc1',
          'O=C(O)CC(O)(CC(=O)O)C(=O)O', 'CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O', 'O=S(=O)(O)c1ccc(N)cc1',
          'CN(C)CCOC(c1ccccc1)c1ccccc1', 'OC(=O)C(F)(F)F', 'ClC(Cl)(Cl)Cl', 'OB(O)c1ccccc1', 'NC(N)=[NH2+]',
          'C[N+](C)(C)CC([O-])=O', 'NCCCC[C@H](N)C(=O)O', 'NC(=N)NCCC[C@H](N)C(=O)O', 'c1ccc2[nH]ccc2c1',
          'OC(=O)c1cc(O)c(O)c(O)c1', 'COc1cc(C=O)ccc1O', 'CC(=O)c1ccc(cc1)S(=O)(=O)N', 'C[S+](C)C', 'CC[O-]',
          '[CH2-]C=C', 'C[CH2+]', 'OS(=O)(=O)O', 'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O']


def _doc(smiles):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text(smiles), (5, 5))
    return d


def _variants(smi):
    rng = random.Random(smi)
    yield 'plain', _doc(smi)
    d = _doc(smi); d.style = dict(show_stereo=True, show_atom_numbers=True); yield 'numbers', d
    d = _doc(smi)
    for k, at in enumerate(d.atoms.values(), 1):
        at.map = k
    yield 'map', d
    d = _doc(smi)
    for at in d.atoms.values():
        if at.el not in ('C', 'H'):
            at.lone_pairs = -1
    yield 'lone pairs', d
    d = _doc(smi)
    for at in d.atoms.values():
        if rng.random() < 0.25 and at.el in ('C', 'N', 'O'):
            at.charge, at.radical = rng.choice([(1, 0), (0, 1), (1, 1), (-1, 0), (0, 2)])
    yield 'ions', d


@pytest.mark.parametrize('smiles', SMILES)
def test_nothing_drawn_on_anything_else(qapp, smiles):
    bad = []
    for name, d in _variants(smiles):
        bad += [f'{name}: {i.kind} {i.text}' for i in overlaps.find(d)]
    assert not bad, '\n'.join(bad)


def test_the_checker_sees_an_overlap(qapp):
    d = _doc('CCO')
    a, b = list(d.atoms)[:2]
    d.atoms[a].el = d.atoms[b].el = 'N'                     # two labels pressed together
    d.atoms[b].x, d.atoms[b].y = d.atoms[a].x + 0.25, d.atoms[a].y
    assert overlaps.find(d)


def test_two_oh_side_by_side_put_their_h_apart(qapp):
    d = _doc('OP(=O)(O)OP(=O)(O)O')
    assert not overlaps.find(d)


def test_layout_spreads_crowded_end_atoms(qapp):
    import layout2d
    import numpy as np
    m = chem.mol_from_text('OP(=O)(O)OP(=O)(O)O')
    layout2d.compute(m)
    conf = m.GetConformer()
    xy = conf.GetPositions()[:, :2]
    bond = float(np.median([layout2d._length(conf, b) for b in m.GetBonds()]))
    need = layout2d._crowd_pairs(m)
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2) / bond
    assert (d[need > 0] >= need[need > 0] - 1e-6).all()    # every pair as far apart as its labels need


def test_labels_with_nicknames_and_charges_turn_round():
    assert model.label_order('OTBS', True) == (['TBS', 'O'], 1)
    assert model.label_order('NHBoc', True) == (['Boc', 'H', 'N'], 2)
    assert model.label_order('NMe3+', True) == (['Me3', 'N', '+'], 1)
    assert model.label_order('OTBS', False)[1] == 0


# ------------------------------------------------------------------------------------------- documents Kekulex makes
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


def test_analysis_captions_of_neighbours_do_not_meet(win):
    c = win.canvas
    d = Doc()
    for k, s in enumerate(('CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                           'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O')):
        d.insert_rdkit(chem.mol_from_text(s), (4 + 7 * k, 4))
    c.set_doc(d)
    c.add_analysis()
    assert not overlaps.find(c.doc)


def test_repeated_cleavages_stack_without_overlap(win):
    c = win.canvas
    d = _doc('CCCC(C)=O')
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    d.atoms[o].charge, d.atoms[o].radical = 1, 1
    c.set_doc(d)
    co = next(a for a, at in d.atoms.items() if at.el == 'C' and len(d.neighbors(a)) == 3)
    for b in [bb for n, bb in d.neighbors(co) if d.atoms[n].el == 'C']:
        for k in range(len(ms.options(c.doc, b))):
            c.cleave_bond(b, k)
    assert not overlaps.find(c.doc)


def test_resonance_forms_are_textbook_and_spaced(win):
    from rdkit import Chem
    d = _doc('[O-][N+](=O)c1ccccc1')
    forms = [Chem.MolToSmiles(m) for m in reactions.resonance_forms(d, set(d.atoms))]
    assert not any('OH2+' in f or '[O+]' in f for f in forms)    # no oxygen short of its octet
    assert sum('[CH+]' in f for f in forms) == 3                  # ortho, ortho, para
    benzene = _doc('c1ccccc1')
    assert len(reactions.resonance_forms(benzene, set(benzene.atoms))) == 1   # the other Kekulé form only
    c = win.canvas
    c.set_doc(_doc('[O-][N+](=O)c1ccccc1'))
    c.set_selection([])
    c.add_forms('resonance')
    assert not overlaps.find(c.doc)


def test_scheme_with_conditions_is_clear(win):
    c = win.canvas
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CC(=O)O'), (3, 4)); d.insert_rdkit(chem.mol_from_text('CCO'), (8, 6))
    ar = d.add_arrow(11, 4, 14, 4); d.insert_rdkit(chem.mol_from_text('CC(=O)OCC'), (18, 4))
    c.set_doc(d)
    c.set_conditions(ar, 'H2SO4 (cat.), Dean-Stark trap', 'toluene, reflux, 12 h')
    c.set_selection([]); c.arrange_scheme(); c.map_reaction()
    assert not overlaps.find(c.doc)


# ------------------------------------------------------------------------------------------- pictures
def _picture(w=40, h=20, color='#C03030'):
    img = QtGui.QImage(w, h, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor(color))
    return img


def test_paste_a_copied_picture(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    QtWidgets.QApplication.clipboard().setImage(_picture())
    win.paste()
    shapes = [s for s in c.doc.shapes.values() if s.kind == 'image']
    assert len(shapes) == 1 and ('shape', next(iter(c.doc.shapes))) in c.sel
    sh = shapes[0]
    assert (sh.x2 - sh.x1) / (sh.y2 - sh.y1) == pytest.approx(2.0)
    img = render.to_image(c.doc, 150, False)               # it is in the exported picture
    reds = sum(1 for x in range(0, img.width(), 3) for y in range(0, img.height(), 3)
               if QtGui.QColor(img.pixel(x, y)).red() > 150 and QtGui.QColor(img.pixel(x, y)).green() < 100)
    assert reds > 20


def test_smiles_text_still_pastes_as_a_structure(win):
    c = win.canvas
    c.set_doc(Doc())
    md = QtCore.QMimeData()
    md.setText('c1ccccc1')
    md.setImageData(_picture())
    QtWidgets.QApplication.clipboard().setMimeData(md)
    win.paste()
    assert len(c.doc.atoms) == 6 and not c.doc.shapes


def test_paste_a_copied_image_file(win, tmp_path):
    path = str(tmp_path / 'logo.png')
    _picture(30, 30, '#2050C0').save(path)
    md = QtCore.QMimeData()
    md.setUrls([QtCore.QUrl.fromLocalFile(path)])
    QtWidgets.QApplication.clipboard().setMimeData(md)
    win.canvas.set_doc(Doc())
    win.paste()
    assert [s.kind for s in win.canvas.doc.shapes.values()] == ['image']


def test_pictures_are_kept_moved_flipped_and_copied(win, tmp_path):
    import pictures
    c = win.canvas
    c.set_doc(Doc())
    i = c.insert_image(_picture())
    back = Doc.from_dict(json.loads(json.dumps(c.doc.to_dict())))
    assert pictures.decode(back.shapes[i].params['png']).size() == QtCore.QSize(40, 20)
    c.set_selection([('shape', i)])
    c.transform_selection('mirror_h')
    assert c.doc.shapes[i].x2 < c.doc.shapes[i].x1          # drawn mirrored
    assert not render.to_image(c.doc, 96, False).isNull()
    win.copy(); win.paste()
    assert sum(1 for s in c.doc.shapes.values() if s.kind == 'image') == 2
    assert c.doc.find_shape((c.doc.shapes[i].x1 + c.doc.shapes[i].x2) / 2,
                            (c.doc.shapes[i].y1 + c.doc.shapes[i].y2) / 2, 0.1) is not None   # click anywhere on it
