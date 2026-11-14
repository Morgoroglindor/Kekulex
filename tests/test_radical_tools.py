"""The charge palette also holds the radical tools: radical •, diradical ••, radical cation •+ and radical anion •−,
each clicked on an atom (clicked again: taken away)."""
import pytest
from PyQt6 import QtCore
from PyQt6.QtTest import QTest

import chem
import icons
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


def click(c, a):
    at = c.doc.atoms[a]
    QTest.mouseClick(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(at.x, at.y).toPoint())


def _setup(win, smiles='CC'):
    c = win.canvas
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles), (4, 4)); c.set_doc(d)
    a = next(a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 1)
    return c, a


def test_the_palette_has_the_radicals_beside_the_charges(win):
    assert 'charges' in win.group_buttons
    assert {'charge+', 'charge-', 'radical', 'diradical', 'radical+', 'radical-'} <= set(win.tool_actions)
    for n in ('radical', 'diradical', 'radical+', 'radical-'):
        assert not icons.icon(n).isNull() and n in icons.DRAW
        win.set_tool(n)
        assert win.canvas.hint()


def test_radical_tool_toggles(win):
    c, a = _setup(win)
    win.set_tool('radical')
    click(c, a)
    assert c.doc.atoms[a].radical == 1 and c.doc.atoms[a].charge == 0
    assert chem.analyse(c.doc).formula.startswith('C2H5')
    click(c, a)
    assert c.doc.atoms[a].radical == 0


def test_diradical_tool(win):
    c, a = _setup(win)
    win.set_tool('diradical')
    click(c, a)
    assert c.doc.atoms[a].radical == 2
    assert chem.analyse(c.doc).formula.startswith('C2H4')
    win.set_tool('radical')                                 # another radical tool: replaces the dots
    click(c, a)
    assert c.doc.atoms[a].radical == 1


def test_radical_cation_and_anion(win):
    c, a = _setup(win, 'C=CC')
    win.set_tool('radical+')
    click(c, a)
    at = c.doc.atoms[a]
    assert (at.radical, at.charge) == (1, 1)
    assert chem.analyse(c.doc).radicals == 1
    win.set_tool('radical-')                                # the anion tool on the cation: becomes the anion
    click(c, a)
    assert (at.radical, at.charge) == (1, -1)
    click(c, a)                                             # again: neutral
    assert (at.radical, at.charge) == (0, 0)


def test_radical_tools_undo(win):
    c, a = _setup(win)
    win.set_tool('radical+')
    click(c, a)
    c.undo()
    assert (c.doc.atoms[a].radical, c.doc.atoms[a].charge) == (0, 0)


def test_radical_tools_ignore_empty_space(win):
    c, _a = _setup(win)
    n = len(c.doc.atoms)
    win.set_tool('radical')
    QTest.mouseClick(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(20, 20).toPoint())
    assert len(c.doc.atoms) == n and all(at.radical == 0 for at in c.doc.atoms.values())
