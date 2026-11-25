"""Stoichiometry grid: grouped structures count as one compound; the grid's cells are typed into in place."""
import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtTest import QTest

import chem
import stoich
from model import Doc


@pytest.fixture
def win(qapp):
    from kekulex import MainWindow
    w = MainWindow()
    yield w
    w.dirty = False
    w.close()


def _esterification(c, group_reagent=True):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CC(=O)O'), (3, 4))
    d.insert_rdkit(chem.mol_from_text('CCO'), (7, 4))
    d.insert_rdkit(chem.mol_from_text('O=S(=O)(O)O'), (7, 6.5))
    d.add_arrow(10, 4, 14, 4)
    d.insert_rdkit(chem.mol_from_text('CC(=O)OCC'), (17, 4))
    c.set_doc(d)
    if group_reagent:
        frags = c.doc.fragments()
        pick = [f for f in frags if chem.analyse(c.doc, f).formula in ('C2H6O', 'H2O4S')]
        c.doc.make_group([('atom', a) for f in pick for a in f])


def test_grouped_structures_are_one_compound(win):
    _esterification(win.canvas)
    reac, prod = stoich.reaction_columns(win.canvas.doc)
    doc = win.canvas.doc
    assert [chem.analyse(doc, f).formula for f in reac] == ['C2H4O2', 'C2H8O5S']
    assert len(prod) == 1


def test_whole_reaction_grouped_still_splits(win):
    c = win.canvas
    _esterification(c, group_reagent=False)
    c.doc.make_group(c.doc.all_items())                     # everything, arrow included: not one compound
    reac, prod = stoich.reaction_columns(c.doc)
    assert len(reac) == 3 and len(prod) == 1


def test_grouping_after_the_grid_merges_its_columns(win):
    c = win.canvas
    _esterification(c, group_reagent=False)
    gid = c.analyze_stoichiometry()
    g = c.doc.grids[gid]
    assert len(g.reactants) == 3
    g.reactants[1]['equiv'] = 2.0
    pick = [f for f in c.doc.fragments() if chem.analyse(c.doc, f).formula in ('C2H6O', 'H2O4S')]
    c.doc.make_group([('atom', a) for f in pick for a in f])
    t = stoich.compute(c.doc, g)
    assert [col.formula for col in t['reactants']] == ['C2H4O2', 'C2H8O5S']
    assert g.reactants[1]['equiv'] == 2.0                   # the typed value is kept


def test_bare_number_is_the_base_unit():
    assert stoich.parse('0.6', 'mass') == pytest.approx(0.6)
    assert stoich.parse('2', 'volume') == pytest.approx(2.0)


def test_typing_in_a_cell_in_place_live_tab_and_escape(win, qapp):
    c = win.canvas
    _esterification(c)
    win.show()
    gid = c.analyze_stoichiometry()
    g = c.doc.grids[gid]
    undo = len(c.undo_stack)
    c.edit_grid_cell(gid, 'reactants', 0, 'mass')
    QTest.keyClicks(c.editor[0], '600 mg')
    assert g.reactants[0]['mass'] == pytest.approx(0.6)     # live, before Enter
    QTest.keyClick(c.editor[0], Qt.Key.Key_Tab)
    assert c.editor[3][:3] == ('reactants', 0, 'pct')       # Tab: the next cell that takes a value
    QTest.keyClick(c.editor[0], Qt.Key.Key_Escape)
    assert 'pct' not in g.reactants[0]
    assert len(c.undo_stack) == undo + 1                    # one step for the mass, none for the cancelled cell
    t = stoich.compute(c.doc, g)
    assert t['reactants'][1].cells['moles'].value == pytest.approx(0.6 / 60.052, rel=1e-3)
    c.edit_grid_cell(gid, 'reactants', 1, 'equiv')
    QTest.keyClicks(c.editor[0], '3')
    QTest.keyClick(c.editor[0], Qt.Key.Key_Return)
    assert c.editor is None and g.reactants[1]['equiv'] == 3.0
    c.edit_grid_cell(gid, 'reactants', 1, 'equiv')
    c.editor[0].setText('')
    c.grid_live('')
    QTest.keyClick(c.editor[0], Qt.Key.Key_Return)
    assert 'equiv' not in g.reactants[1]                    # emptied: calculated again


def test_double_click_header_sets_limiting(win, qapp):
    c = win.canvas
    _esterification(c)
    win.show()
    gid = c.analyze_stoichiometry()
    c.fit()
    rend = c.renderer()
    cells = rend.grid_layout(gid)[2]
    r = cells[('reactants', 1, 'header')]
    pos = c.rect_px(r).center()
    c.grid_click(gid, pos)
    g = c.doc.grids[gid]
    assert g.reactants[1]['limiting'] and not g.reactants[0]['limiting']
    r = c.renderer().grid_layout(gid)[2][('reactants', 0, 'mass')]
    c.grid_click(gid, c.rect_px(r).center())
    assert c.editor is not None and c.editor[1] == 'grid'
    c.commit_editor(cancel=True)
