"""t over an atom, live analysis captions, autosave / recovery, lower-case labels, groups in .cdxml, polymer
brackets."""
import json
import os

import pytest
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtTest import QTest

import autosave
import cdxml
import chem
import model
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


def _doc(smiles, at=(3, 3)):
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles), at); return d


def _end_atom(doc):
    return next(a for a in doc.atoms if len(doc.neighbors(a)) == 1)


# ------------------------------------------------------------------------------------------- t over an atom
def test_t_over_an_atom_types_its_label(win):
    c = win.canvas
    c.set_doc(_doc('CCC'))
    a = _end_atom(c.doc)
    c.set_tool('bond1')
    c.hover = ('atom', a)
    QTest.keyClick(c, 't')
    assert c.editor and c.editor[1] == 'label' and c.editor[2] == a
    c.editor[0].setText('SO3K')
    c.commit_editor()
    assert c.doc.atoms[a].label == 'SO3K'


def test_text_tool_over_an_atom_takes_what_is_typed(win):
    c = win.canvas
    c.set_doc(_doc('CCC'))
    a = _end_atom(c.doc)
    c.set_tool('text')
    c.hover = ('atom', a)
    QTest.keyClick(c, 'o')                                  # the first letter starts the label
    assert c.editor and c.editor[0].text() == 'o'
    QTest.keyClicks(c.editor[0], 'me')
    c.commit_editor()
    assert c.doc.atoms[a].label == 'OMe'


def test_t_away_from_atoms_is_still_the_text_tool(win):
    c = win.canvas
    c.set_doc(Doc())
    picked = []
    c.tool_key.connect(picked.append)
    c.hover = None
    QTest.keyClick(c, 't')
    assert picked == ['text']


# ------------------------------------------------------------------------------------------- live captions
def test_pasted_analysis_follows_the_structure(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    win.info.refresh()
    win.info.paste_all()
    t = next(iter(c.doc.texts))
    assert 'C2H6O' in c.doc.texts[t].text
    o = next(a for a, at in c.doc.atoms.items() if at.el == 'O')
    c.snapshot(); c.doc.atoms[o].el = 'N'; c.edited()       # ethanol -> ethylamine
    assert 'C2H7N' in c.doc.texts[t].text and 'C2H6O' not in c.doc.texts[t].text
    c.snapshot(); c.doc.sprout(o); c.edited()               # a carbon more
    assert 'C3H9N' in c.doc.texts[t].text
    c.undo()
    assert 'C2H7N' in c.doc.texts[t].text


def test_pasted_properties_and_name_follow_the_structure(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    win.info.tabs.setCurrentIndex(1)
    win.info.refresh(keep_choice=True)
    win.info.props.item(3, 0).setCheckState(QtCore.Qt.CheckState.Checked)       # H-bond donors
    win.info.paste_properties()
    t = next(i for i, x in c.doc.texts.items() if x.text.startswith('H-bond donors'))
    assert c.doc.texts[t].text == 'H-bond donors: 1'
    o = next(a for a, at in c.doc.atoms.items() if at.el == 'O')
    c.snapshot(); c.doc.atoms[o].el = 'C'; c.edited()       # propane: no donors
    assert c.doc.texts[t].text == 'H-bond donors: 0'


def test_show_analysis_caption_is_live_and_typing_over_it_stops_that(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    c.add_analysis()
    t = next(iter(c.doc.texts))
    assert c.doc.texts[t].live
    c.edit_text(t)
    c.editor[0].setPlainText('my note')
    c.commit_editor()
    assert c.doc.texts[t].live == []
    o = next(a for a, at in c.doc.atoms.items() if at.el == 'O')
    c.snapshot(); c.doc.atoms[o].el = 'N'; c.edited()
    assert c.doc.texts[t].text == 'my note'


def test_live_captions_survive_saving():
    d = _doc('CCO')
    t = d.add_text(3, 6, 'x', live=['formula'])
    back = Doc.from_dict(json.loads(json.dumps(d.to_dict())))
    assert back.texts[t].live == ['formula']


# ------------------------------------------------------------------------------------------- autosave
def test_autosave_writes_a_copy_and_clears_it_on_save(win, tmp_path):
    win.autosave = autosave.Slot(str(tmp_path))
    c = win.canvas
    c.set_doc(_doc('CCO')); c.snapshot(); c.edited()
    win.autosave_now()
    assert os.path.exists(win.autosave.sed)
    with open(win.autosave.sed, encoding='utf-8') as f:
        assert len(Doc.from_dict(json.load(f)).atoms) == 3
    win.path = str(tmp_path / 'drawing.sed')
    assert win.save()
    assert not os.path.exists(win.autosave.sed)


def test_recovery_offers_a_copy_left_by_a_dead_program(win, tmp_path, monkeypatch):
    slot = autosave.Slot(str(tmp_path))
    slot.write(_doc('c1ccccc1'), r'C:\somewhere\benzene.sed')
    meta = json.load(open(slot.meta, encoding='utf-8'))
    meta['pid'] = 999999                                     # a program that is no longer running
    json.dump(meta, open(slot.meta, 'w', encoding='utf-8'))
    assert [o[0] for o in autosave.orphans(str(tmp_path))] == [slot.sed]
    monkeypatch.setattr(autosave, 'DIR', str(tmp_path))
    monkeypatch.setattr(QtWidgets.QMessageBox, 'exec', lambda self: 0)
    monkeypatch.setattr(QtWidgets.QMessageBox, 'clickedButton', lambda self: self.buttons()[0])  # Recover
    win.autosave = autosave.Slot(str(tmp_path / 'mine'))
    assert win.recover()
    assert len(win.canvas.doc.atoms) == 6 and win.dirty
    assert win.path == r'C:\somewhere\benzene.sed'
    assert not os.path.exists(slot.sed)                     # offered once


def test_a_running_programs_copy_is_not_an_orphan(tmp_path):
    slot = autosave.Slot(str(tmp_path))
    slot.write(_doc('CC'))
    assert autosave.orphans(str(tmp_path)) == []


# ------------------------------------------------------------------------------------------- lower-case labels
@pytest.mark.parametrize('typed, meant', [
    ('so3k', 'SO3K'), ('ch2cl', 'CH2Cl'), ('co2et', 'CO2Et'), ('n(ch3)2', 'N(CH3)2'), ('cf3', 'CF3'),
    ('no2', 'NO2'), ('boc', 'Boc'), ('ome', 'OMe'), ('so3na', 'SO3Na'), ('cooh', 'COOH'), ('nme2', 'NMe2')])
def test_lower_case_labels(typed, meant):
    assert model.parse_label(typed) == model.parse_label(meant)


def test_lower_case_single_elements_and_unknowns():
    assert model.parse_label('nh2') == ('element', 'N', 0)  # not nihonium
    assert model.parse_label('oh') == ('element', 'O', 0)
    assert model.parse_label('cl') == ('element', 'Cl', 0)
    assert model.parse_label('xyz') is None
    assert model.parse_label('Co') == ('element', 'Co', 0)  # capitals are still read as written


# ------------------------------------------------------------------------------------------- groups in .cdxml
def test_groups_round_trip_through_cdxml(tmp_path):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('C[NH3+]'), (2, 2))
    d.insert_rdkit(chem.mol_from_text('[Cl-]'), (5, 2))
    d.insert_rdkit(chem.mol_from_text('c1ccccc1'), (12, 2))
    salt = [f for f in d.fragments() if len(f) <= 2]
    d.make_group([('atom', a) for f in salt for a in f])
    t = d.add_text(3, 4, 'salt')
    d.join_group(('text', t), salt[0])
    path = str(tmp_path / 'g.cdxml')
    cdxml.write(d, path)
    assert open(path, encoding='utf-8').read().count('<group ') == 1
    back, _notes = cdxml.read(path)
    assert len(back.groups) == 1
    g = next(iter(back.groups))
    items = back.group_items(g)
    assert len({a for k, a in items if k == 'atom'}) == 3 and any(k == 'text' for k, _i in items)
    benzene = [f for f in back.fragments() if len(f) == 6][0]
    assert back.group_of(('atom', benzene[0])) is None


# ------------------------------------------------------------------------------------------- polymer brackets
def _polyethylene(rep):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CCCC'), (2, 2))
    xs = sorted(d.atoms, key=lambda a: d.atoms[a].x)
    mid = xs[1:3]
    x1 = (d.atoms[xs[0]].x + d.atoms[xs[1]].x) / 2
    x2 = (d.atoms[xs[2]].x + d.atoms[xs[3]].x) / 2
    ys = [d.atoms[a].y for a in d.atoms]
    s = d.add_shape(x1, min(ys) - 0.8, x2, max(ys) + 0.8, 'brackets')
    if rep:
        d.shapes[s].params['repeat'] = rep
    return d, s, mid


def test_polymer_brackets_with_a_number_count_the_unit():
    d, _s, _mid = _polyethylene('3')
    info = chem.analyse(d)
    assert info.formula == 'C8H18'                          # CH3-(CH2CH2)3-CH3
    plain, _s, _m = _polyethylene('')
    assert chem.analyse(plain).formula == 'C4H10'


def test_polymer_brackets_with_n_show_the_repeat_unit():
    d, _s, _mid = _polyethylene('n')
    info = chem.analyse(d)
    assert info.formula == 'C2H6(C2H4)n'
    assert info.repeat == 'n'
    assert '(n = 1)' in chem.live_caption(d, None, ['exact'])


def test_polymer_repeat_label_is_drawn(win):
    import render
    d, s, _mid = _polyethylene('n')
    r = render.Renderer(d)
    path = r.repeat_path(s)
    assert path is not None and not path.isEmpty()
    assert r.item_rect('shape', s).right() >= path.boundingRect().right() - 0.01
    d.shapes[s].params.pop('repeat')
    assert render.Renderer(d).repeat_path(s) is None
