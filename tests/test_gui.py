"""The window with real mouse and key events, exports and the clipboard."""
import json

import pytest
from PyQt6 import QtCore
from PyQt6.QtTest import QTest

import chem
import render
from model import Doc


@pytest.fixture
def win(qapp):
    import structure_editor as se
    w = se.MainWindow()
    w.resize(1200, 800)
    w.show()
    qapp.processEvents()
    yield w
    w.dirty = False
    w.close()


def click(c, x, y):
    QTest.mouseClick(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(x, y).toPoint())


def drag(c, x0, y0, x1, y1):
    QTest.mousePress(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(x0, y0).toPoint())
    QTest.mouseMove(c, c.to_px((x0 + x1) / 2, (y0 + y1) / 2).toPoint())
    QTest.mouseMove(c, c.to_px(x1, y1).toPoint())
    QTest.mouseRelease(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(x1, y1).toPoint())


def test_draw_phenol_with_mouse_and_keys(win):
    c = win.canvas
    win.set_tool('benzene'); click(c, 0, 0)
    top = min(c.doc.atoms, key=lambda a: c.doc.atoms[a].y)
    win.set_tool('bond1'); click(c, c.doc.atoms[top].x, c.doc.atoms[top].y)
    new = max(c.doc.atoms)
    c.hover = ('atom', new)
    QTest.keyClick(c, QtCore.Qt.Key.Key_O)
    assert chem.analyse(c.doc).smiles == 'Oc1ccccc1'
    c.undo()
    assert chem.analyse(c.doc).smiles == 'Cc1ccccc1'
    c.redo()
    assert chem.analyse(c.doc).formula == 'C6H6O'


def test_bond_click_cycles_order(win):
    c = win.canvas
    win.set_tool('bond1'); click(c, 0, 0)
    b = next(iter(c.doc.bonds))
    bd = c.doc.bonds[b]
    mid = ((c.doc.atoms[bd.a].x + c.doc.atoms[bd.b].x) / 2, (c.doc.atoms[bd.a].y + c.doc.atoms[bd.b].y) / 2)
    click(c, *mid); assert bd.order == 2
    click(c, *mid); assert bd.order == 3
    click(c, *mid); assert bd.order == 1


def test_drag_chain_and_close_ring_by_drag(win):
    c = win.canvas
    win.set_tool('chain'); drag(c, 0, 0, 5 * 0.866, 0)
    assert len(c.doc.bonds) == 5
    win.set_tool('bond1')
    ends = [a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 1]
    a, b = ends
    drag(c, c.doc.atoms[a].x, c.doc.atoms[a].y, c.doc.atoms[b].x, c.doc.atoms[b].y)
    assert len(c.doc.bonds) == 6 and len(c.doc.atoms) == 6    # the drag ended on an atom: a ring, no new atom


def test_wedge_gives_stereo(win):
    c = win.canvas
    c.doc.insert_rdkit(chem.mol_from_text('CC(N)C(=O)O'))
    c.edited()
    centre = next(a for a in c.doc.atoms if len(c.doc.neighbors(a)) == 3 and c.doc.atoms[a].el == 'C'
                  and all(c.doc.bonds[b].order == 1 for _n, b in c.doc.neighbors(a)))
    methyl = next(n for n, _b in c.doc.neighbors(centre) if c.doc.atoms[n].el == 'C' and len(c.doc.neighbors(n)) == 1)
    bid = c.doc.bond_between(centre, methyl)
    if c.doc.bonds[bid].a != centre:
        bd = c.doc.bonds[bid]; bd.a, bd.b = bd.b, bd.a
    win.set_tool('wedge')
    bd = c.doc.bonds[bid]
    click(c, (c.doc.atoms[bd.a].x * 0.6 + c.doc.atoms[bd.b].x * 0.4), (c.doc.atoms[bd.a].y * 0.6 + c.doc.atoms[bd.b].y * 0.4))
    assert c.doc.bonds[bid].stereo == 'wedge'
    assert '@' in chem.analyse(c.doc).smiles


def test_select_move_delete(win):
    c = win.canvas
    win.set_tool('ring6'); click(c, 0, 0); click(c, 6, 0)
    win.set_tool('select')
    QTest.mousePress(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(-2, -2).toPoint())
    QTest.mouseMove(c, c.to_px(0, 0).toPoint())
    QTest.mouseMove(c, c.to_px(2, 2).toPoint())
    QTest.mouseRelease(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(2, 2).toPoint())
    assert len(c.selected_atoms()) == 6
    QTest.keyClick(c, QtCore.Qt.Key.Key_Delete)
    assert len(c.doc.atoms) == 6 and len(c.doc.fragments()) == 1


def test_smiles_box_and_info_panel(win):
    win.smiles_box.setText('Cn1cnc2c1c(=O)n(C)c(=O)n2C')
    win.insert_smiles_box()
    win.info.refresh()
    assert win.info.formula.text() == 'C<sub>8</sub>H<sub>10</sub>N<sub>4</sub>O<sub>2</sub>'
    assert win.info.exact.text() == '194.08038'
    assert win.info.table.item(0, 0).text() == '[M+H]+' and win.info.table.item(0, 1).text() == '195.0877'


def test_copy_paste(win, qapp):
    c = win.canvas
    c.doc.insert_rdkit(chem.mol_from_text('CCO'))
    c.edited()
    c.select_all()
    win.copy()
    md = qapp.clipboard().mimeData()
    assert md.text() == 'CCO'
    assert json.loads(bytes(md.data('application/x-structure-editor')))['atoms']
    win.paste()
    assert len(c.doc.fragments()) == 2


def test_flip_keeps_molecule(win):
    c = win.canvas
    c.doc.insert_rdkit(chem.mol_from_text('C[C@@H](N)C(=O)O'))
    c.edited()
    before = chem.analyse(c.doc).smiles
    c.transform_selection('flip_h')
    assert chem.analyse(c.doc).smiles == before
    c.clean_up()
    assert chem.analyse(c.doc).smiles == before


def test_exports(tmp_path, qapp):
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('OC(=O)C1CC(=CC=O)C=C(N1)C(=O)O'))
    doc.add_arrow(4, 0, 6, 0)
    doc.add_text(7, 0, 'betalamic acid')
    render.save_png(doc, str(tmp_path / 'a.png'), dpi=300)
    render.save_svg(doc, str(tmp_path / 'a.svg'))
    render.save_pdf(doc, str(tmp_path / 'a.pdf'))
    from PyQt6.QtGui import QImage
    img = QImage(str(tmp_path / 'a.png'))
    assert img.width() > 300 and img.height() > 100
    svg = (tmp_path / 'a.svg').read_text(encoding='utf-8')
    assert '<svg' in svg and '<path' in svg
    assert (tmp_path / 'a.pdf').read_bytes().startswith(b'%PDF')


def test_open_mol_and_save_sed(win, tmp_path):
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('c1ccncc1'))
    block, _ = doc.to_molblock()
    p = tmp_path / 'pyr.mol'
    p.write_text(block)
    win.open(str(p))
    assert chem.analyse(win.canvas.doc).formula == 'C5H5N'
    win.path = str(tmp_path / 'pyr.sed')
    assert win.save()
    win.open(str(tmp_path / 'pyr.sed'))
    assert chem.analyse(win.canvas.doc).smiles == 'c1ccncc1'


def key(c, text, item):
    c.hover = item
    QTest.keyClicks(c, text)


@pytest.mark.parametrize('keys, smiles', [('O', 'COC'), ('N', 'C[N+](=O)[O-]'), ('P', 'Cc1ccccc1'),
                                          ('2', 'CC=O'), ('4', 'CCC'), ('9', 'CC(C)C'), ('K', 'CCC(C)(C)C'),
                                          ('l', 'CCl'), ('C', 'CCl'), ('b', 'CBr'), ('y', 'CC(=O)OC(C)(C)C')])
def test_chemdraw_atom_hotkeys(win, keys, smiles):
    """Ethane, then the key over its last atom (letters replace the atom, digits add a group to it)."""
    c = win.canvas
    win.set_tool('bond1'); click(c, 0, 0)
    end = max(c.doc.atoms)
    key(c, keys, ('atom', end))
    assert chem.analyse(c.doc).smiles == chem.analyse(_doc(smiles)).smiles


def _doc(smiles):
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles)); return d


def test_chemdraw_bond_hotkeys(win):
    c = win.canvas
    win.set_tool('ring6'); click(c, 0, 0)
    b = next(iter(c.doc.bonds))
    key(c, 'a', ('bond', b))                               # fused benzene
    assert len(c.doc.atoms) == 10
    key(c, '2', ('bond', b))
    assert c.doc.bonds[b].order == 2
    key(c, 'l', ('bond', b)); assert c.doc.bonds[b].side == -1
    key(c, 'c', ('bond', b)); assert c.doc.bonds[b].side == 0
    key(c, 'w', ('bond', b)); assert c.doc.bonds[b].stereo == 'wedge' and c.doc.bonds[b].order == 1


def test_chemdraw_tool_hotkeys(win):
    c = win.canvas
    for k, tool in (('j', 'benzene'), ('x', 'bond1'), ('X', 'chain'), ('e', 'arrow'), ('t', 'text'), (' ', 'select')):
        key(c, k, None)
        assert c.tool == tool and win.tool_actions[tool].isChecked()
