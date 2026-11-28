"""ChemDraw-style dragging (Ctrl: a copy, Shift: straight across or down, both: a straight copy) and captions whose
width can be dragged so their lines wrap (the m/z line in two)."""
import pytest
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QMouseEvent
from PyQt6.QtWidgets import QApplication

import chem
import render
from model import Doc

CTRL, SHIFT, NONE = Qt.KeyboardModifier.ControlModifier, Qt.KeyboardModifier.ShiftModifier, Qt.KeyboardModifier.NoModifier
LEFT = Qt.MouseButton.LeftButton


@pytest.fixture
def win(qapp):
    from kekulex import MainWindow
    w = MainWindow()
    w.resize(1200, 800)
    w.show()
    yield w
    w.dirty = False
    w.close()


def send(c, kind, pos, mods=NONE):
    t = {'press': QMouseEvent.Type.MouseButtonPress, 'move': QMouseEvent.Type.MouseMove,
         'release': QMouseEvent.Type.MouseButtonRelease}[kind]
    buttons = Qt.MouseButton.NoButton if kind == 'release' else LEFT
    ev = QMouseEvent(t, QPointF(pos), QPointF(c.mapToGlobal(pos.toPoint())), LEFT, buttons, mods)
    QApplication.sendEvent(c, ev)


def drag(c, start, end, mods=NONE):
    send(c, 'press', start, mods)
    steps = 6
    for k in range(1, steps + 1):
        p = start + (end - start) * (k / steps)
        send(c, 'move', p, mods)
    send(c, 'release', end, mods)


def setup(win):
    c = win.canvas
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('c1ccccc1O'), (8, 6))
    c.set_doc(d)
    c.set_tool('select') if hasattr(c, 'set_tool') else win.set_tool('select')
    c.select_fragment(next(iter(d.atoms)))
    a = next(iter(d.atoms))
    return c, a


def atom_px(c, a):
    at = c.doc.atoms[a]
    return c.to_px(at.x, at.y)


def test_plain_drag_moves(win):
    c, a = setup(win)
    n = len(c.doc.atoms)
    p = atom_px(c, a)
    drag(c, p, p + QPointF(80, 40))
    assert len(c.doc.atoms) == n


def test_ctrl_drag_makes_a_selected_copy(win):
    c, a = setup(win)
    n = len(c.doc.atoms)
    before = {k: (at.x, at.y) for k, at in c.doc.atoms.items()}
    p = atom_px(c, a)
    drag(c, p, p + QPointF(120, 60), CTRL)
    assert len(c.doc.atoms) == 2 * n
    assert all((c.doc.atoms[k].x, c.doc.atoms[k].y) == xy for k, xy in before.items())   # the original stays
    sel = c.selected_atoms()
    assert len(sel) == n and not sel & set(before)                              # the copy is selected
    assert len(c.doc.fragments()) == 2
    c.undo() if hasattr(c, 'undo') else win.undo()
    assert len(c.doc.atoms) == n                                                # one undo step


def test_shift_drag_stays_straight(win):
    c, a = setup(win)
    y0 = c.doc.atoms[a].y
    x0 = c.doc.atoms[a].x
    p = atom_px(c, a)
    drag(c, p, p + QPointF(150, 30), SHIFT)                 # mostly across: only across
    assert c.doc.atoms[a].y == pytest.approx(y0, abs=1e-9)
    assert c.doc.atoms[a].x > x0 + 1
    p = atom_px(c, a)
    x1 = c.doc.atoms[a].x
    drag(c, p, p + QPointF(20, -140), SHIFT)                # mostly up: only up
    assert c.doc.atoms[a].x == pytest.approx(x1, abs=1e-9)


def test_ctrl_shift_drag_is_a_straight_copy(win):
    c, a = setup(win)
    n = len(c.doc.atoms)
    y0 = c.doc.atoms[a].y
    p = atom_px(c, a)
    drag(c, p, p + QPointF(140, 25), CTRL | SHIFT)
    assert len(c.doc.atoms) == 2 * n
    ys = sorted(round(c.doc.atoms[k].y, 6) for k in c.selected_atoms())
    orig = sorted(round(at.y, 6) for k, at in c.doc.atoms.items() if k not in c.selected_atoms())
    assert ys == orig                                       # the copy sits level with the original


def test_wrap_line_breaks_after_commas():
    line = 'm/z: 300.03 (100.0%), 301.03 (16.5%), 302.03 (2.1%), 303.04 (0.2%)'
    out = render.wrap_line(line, 40, len)
    assert len(out) == 2 and out[0].endswith('(16.5%),') and out[1].startswith('302.03')
    assert render.wrap_line('short', 100, len) == ['short']


def test_caption_width_handle_wraps_the_mz_line(win):
    c, _a = setup(win)
    atoms = set(c.doc.atoms)
    c.add_caption_below(atoms, chem.live_caption(c.doc, atoms, ['mz']), formula=True, live=['mz'])
    t = next(i for i, tx in c.doc.texts.items() if tx.text.startswith('m/z'))
    c.set_selection([('text', t)])
    rend = c.renderer()
    one = rend.text_rect(t)
    _i, handle = c._width_handle(rend)
    hp = c.rect_px(handle).center()
    drag(c, hp, hp - QPointF(one.width() * c.px_per_pt() * 0.45, 0))
    tx = c.doc.texts[t]
    assert tx.width > 0
    two = c.renderer().text_rect(t)
    assert two.height() > one.height() * 1.8 and two.width() < one.width()
    assert tx.text.startswith('m/z')                         # the text itself is unchanged: only drawn in two lines
    c.edited()                                              # a live caption keeps its width when it updates
    assert c.doc.texts[t].width == tx.width
    hp = c.rect_px(c._width_handle(c.renderer())[1]).center()
    drag(c, hp, hp + QPointF(2000, 0))                      # dragged out past its longest line: one line again
    assert c.doc.texts[t].width == 0.0
    d2 = Doc.from_dict(c.doc.to_dict()) if hasattr(Doc, 'from_dict') else None
    if d2 is not None:
        assert all(hasattr(x, 'width') for x in d2.texts.values())
