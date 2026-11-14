"""Pressing an element's key again over its atom takes one H away each time, then gives them all back - as ChemDraw
23.1 does (checked there: N on a chain end NH2 -> NH -> N -> NH2, O: OH -> O -> OH, C: CH3 -> CH2 -> CH -> C -> CH3,
N in a chain NH -> N -> NH; a C pressed on a plain vertex first shows CH3; fewer H than usual get a red box)."""
import pytest
from PyQt6.QtTest import QTest

import chem
import render
from model import Doc


@pytest.fixture
def win(qapp):
    import kekulex as se
    w = se.MainWindow()
    w.resize(1100, 800)
    w.show()
    qapp.processEvents()
    yield w
    w.dirty = False
    w.close()


def _shown(c, a):
    """The label drawn on atom a, as one string ('' for a plain vertex)."""
    lab = render.Renderer(c.doc).labels.get(a)
    return ''.join(p[0] for p in lab['pieces']) if lab else ''


def _press(c, a, key, times=1):
    c.hover = ('atom', a)
    out = []
    for _ in range(times):
        QTest.keyClick(c, key)
        c.hover = ('atom', a)
        out.append((_shown(c, a), c.doc.valence_error(a)))
    return out


def _chain_end(c):
    return next(a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 1)


def test_n_on_a_chain_end_cycles_its_h(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CC'), (4, 4)); c.set_doc(d)
    a = _chain_end(c)
    assert _press(c, a, 'n', 4) == [('NH2', False), ('NH', True), ('N', True), ('NH2', False)]


def test_o_cycles_oh_o_oh(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CC'), (4, 4)); c.set_doc(d)
    a = _chain_end(c)
    assert _press(c, a, 'o', 3) == [('OH', False), ('O', True), ('OH', False)]


def test_c_shows_ch3_then_cycles_down(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CC'), (4, 4)); c.set_doc(d)
    a = _chain_end(c)
    assert _press(c, a, 'c', 5) == [('CH3', False), ('CH2', True), ('CH', True), ('C', True), ('CH3', False)]


def test_n_in_a_chain(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CCC'), (4, 4)); c.set_doc(d)
    mid = next(a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 2)
    shown = [s for s, _e in _press(c, mid, 'n', 3)]
    assert [s.replace('H', 'H') for s in shown] == ['NH', 'N', 'NH'] or shown == ['HN', 'N', 'HN']


def test_l_dopa_amine_and_phenol(win):
    """L-DOPA: N repeatedly over its NH2 gives NH then N (a nitrene-like N, boxed red), then NH2 again; the same for
    O over a phenol OH. The formula follows each step."""
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O'), (6, 5)); c.set_doc(d)
    n = next(a for a, at in c.doc.atoms.items() if at.el == 'N')
    formulas = []
    c.hover = ('atom', n)
    for _ in range(3):
        QTest.keyClick(c, 'n'); c.hover = ('atom', n)
        formulas.append(chem.analyse(c.doc).formula)
    assert formulas == ['C9H10NO4', 'C9H9NO4', 'C9H11NO4']
    phenol = next(a for a, at in c.doc.atoms.items() if at.el == 'O' and len(c.doc.neighbors(a)) == 1
                  and c.doc.atoms[c.doc.neighbors(a)[0][0]].el == 'C' and c.doc.implicit_h(a) == 1
                  and any(c.doc.bonds[b].order == 2 for _n, b in c.doc.neighbors(c.doc.neighbors(a)[0][0])))
    assert [s for s, _e in _press(c, phenol, 'o', 2)] in (['O', 'OH'], ['O', 'HO'])


def test_another_element_resets_the_h_and_undo_steps_back(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CC'), (4, 4)); c.set_doc(d)
    a = _chain_end(c)
    _press(c, a, 'n', 2)                                     # NH
    assert _press(c, a, 'o')[0] == ('OH', False)             # O with all its H again
    c.undo()
    assert c.doc.atoms[a].el == 'N' and c.doc.atoms[a].hcount == 1


def test_check_structure_names_the_open_valence(win):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text('CC'), (4, 4)); c.set_doc(d)
    a = _chain_end(c)
    _press(c, a, 'n', 2)
    assert any('open valence' in p for p in c.check_structure())
