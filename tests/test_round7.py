"""Drawing elements (styles, arcs), orbitals, TLC plates, mass fragmentation and (retro)synthesis."""
import pytest
from PyQt6 import QtCore, QtGui
from rdkit import Chem

import canvas
import cdxml
import elements
import fragment
import render
import tlc
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


def _mol(smiles):
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(smiles))
    return d


def _cut_line(d, i, reach=0.8):
    """A line crossing bond i at its middle."""
    b = d.bonds[i]
    pa, pb = d.atoms[b.a], d.atoms[b.b]
    mx, my = (pa.x + pb.x) / 2, (pa.y + pb.y) / 2
    dx, dy = pb.x - pa.x, pb.y - pa.y
    return (mx - dy * reach, my + dx * reach, mx + dy * reach, my - dx * reach)


def test_element_tools_cover_chemdraws_palette():
    names = canvas.ELEMENT_TOOLS
    assert len(names) == 32                                                        # 4 rows of 8
    for kind in ('circle', 'oval', 'roundrect', 'rect'):
        assert {canvas.SHAPE_TOOLS[n][1] for n in names if canvas.SHAPE_TOOLS[n][0] == kind} == \
            {'', 'dashed', 'shaded', 'filled', 'shadow'}
    assert {canvas.SHAPE_TOOLS[n][1] for n in names if canvas.SHAPE_TOOLS[n][0] == 'line'} == \
        {'', 'dashed', 'wavy', 'bold'}
    arcs = {(canvas.SHAPE_TOOLS[n][2]['angle'], canvas.SHAPE_TOOLS[n][1]) for n in names if n.startswith('arc')}
    assert arcs == {(a, s) for a in (90, 120, 180, 270) for s in ('', 'dashed')}
    assert len(canvas.ORBITAL_TOOLS) == 21


def test_elements_draw_and_hit(qapp):
    d = Doc()
    for n, name in enumerate(canvas.ELEMENT_TOOLS + canvas.ORBITAL_TOOLS):
        kind, style, params = canvas.SHAPE_TOOLS[name]
        x = n * 2.0
        i = d.add_shape(x, 0, x + 0.8, 0.6 if kind not in ('circle', 'orbital') else 0, kind)
        d.shapes[i].style, d.shapes[i].params = style, dict(params)
    r = render.Renderer(d)
    assert render.to_image(d, dpi=40).width() > 200
    for i in d.shapes:
        assert not r.shape_paths(i).boundingRect().isEmpty()
    c = d.add_shape(50, 50, 51, 50, 'circle'); d.shapes[c].style = 'filled'
    assert d.find_shape(50.2, 50.1, 0.1) == c                                      # filled: the inside counts
    p = d.add_shape(60, 60, 60, 59, 'orbital'); d.shapes[p].params = {'type': 'p'}
    assert d.find_shape(60, 59.5, 0.1) == p


def test_arc_bulges_up_when_dragged_left_to_right():
    d = Doc()
    i = d.add_shape(0, 0, 2, 0, 'arc'); d.shapes[i].params = {'angle': 180}
    pts = elements.arc_points(d.shapes[i], 10)
    assert abs(pts[0][0]) < 1e-6 and abs(pts[-1][0] - 20) < 1e-6
    assert min(y for _x, y in pts) < -9.9                                          # a half circle, above (y down)
    j = d.add_shape(2, 0, 0, 0, 'arc'); d.shapes[j].params = {'angle': 90}
    assert max(y for _x, y in elements.arc_points(d.shapes[j], 10)) > 0           # right to left: below


def test_mass_fragmentation_homolytic():
    d = _mol('CCc1ccccc1')
    ethyl = next(i for i, b in d.bonds.items() if len(d.neighbors(b.a)) == 1 or len(d.neighbors(b.b)) == 1)
    frags = sorted(fragment.mass_fragments(d, _cut_line(d, ethyl)), key=lambda f: f[1])
    assert [f[0] for f in frags] == ['CH3', 'C7H7']
    assert abs(frags[0][1] - 15.0235) < 1e-3 and abs(frags[1][1] - 91.0548) < 1e-3
    ring = next(i for i, b in d.bonds.items() if d.ring_of_bond(i))
    one = fragment.mass_fragments(d, _cut_line(d, ring, 0.3))
    assert len(one) == 1 and one[0][0] == 'C8H10'                                 # a ring bond alone: one piece
    assert fragment.mass_fragments(d, (100, 100, 101, 101)) == []


def test_retrosynthesis_and_synthesis():
    d = _mol('CC(=O)OCc1ccccc1')
    before = len(d.fragments())
    ester = next(i for i, b in d.bonds.items() if {d.atoms[b.a].el, d.atoms[b.b].el} == {'C', 'O'} and b.order == 1
                 and any(d.bonds[j].order == 2 for _n, j in d.neighbors(b.a) + d.neighbors(b.b)))
    new = fragment.disconnect(d, _cut_line(d, ester), retro=True)
    assert any(k == 'arrow' and d.arrows[i].kind == 'retro' for k, i in new)
    assert len(d.fragments()) == before + 2
    import chem
    formulas = sorted(chem.analyse(d, f).formula for f in d.fragments())
    assert formulas == ['C2H4O', 'C7H8O', 'C9H10O2']                               # the pieces get their H back
    d2 = _mol('CC(=O)OCc1ccccc1')
    new = fragment.disconnect(d2, _cut_line(d2, ester), retro=False)
    arrow = next(d2.arrows[i] for k, i in new if k == 'arrow')
    assert arrow.kind == 'forward' and arrow.x2 < min(a.x for a in d2.atoms.values() if True) + 20


def test_tlc_plate(qapp):
    d = Doc()
    i = tlc.new_plate(d, 0, 0, 3.0, 4.0)
    plate = d.tlcs[i]
    assert len(plate.lanes) == 3
    r = render.Renderer(d)
    L = tlc.layout(r, i)
    assert L['front'] < L['origin']
    sr = L['spots'][(0, 0)]
    assert tlc.spot_at(r, i, sr.center().x(), sr.center().y()) == (0, 0)
    assert abs(tlc.rf_at(r, i, L['origin']) - 0) < 1e-9 and abs(tlc.rf_at(r, i, L['front']) - 1) < 1e-9
    assert abs(tlc.rf_at(r, i, sr.center().y()) - 0.5) < 1e-6
    plate.lanes[1][0]['show_rf'] = True
    assert render.to_image(d, dpi=40).width() > 20
    import json
    again = Doc.from_dict(json.loads(json.dumps(d.to_dict())))
    assert again.tlcs[i].lanes == plate.lanes


def test_cdxml_keeps_styles(tmp_path, qapp):
    d = Doc()
    for k, (kind, style) in enumerate((('rect', 'shaded'), ('roundrect', 'dashed'), ('oval', 'filled'),
                                       ('circle', 'shadow'), ('arc', ''), ('orbital', 'shaded'))):
        i = d.add_shape(k * 3, 0, k * 3 + 1, 1 if kind != 'circle' else 0, kind)
        d.shapes[i].style = style
    out = tmp_path / 's.cdxml'
    cdxml.write(d, str(out))
    back, _notes = cdxml.read(str(out))
    got = sorted((s.kind, s.style) for s in back.shapes.values())
    assert got == [('oval', 'filled'), ('oval', 'shadow'), ('rect', 'shaded'), ('roundrect', 'dashed')]


def test_window_tools(win, qapp):
    c = win.canvas
    c.set_doc(_mol('CCc1ccccc1'))
    atom = next(iter(c.doc.atoms))
    win.set_tool('orb:p:shaded')
    c.click(('atom', atom), c.doc.atoms[atom].x, c.doc.atoms[atom].y)
    orb = [s for s in c.doc.shapes.values() if s.kind == 'orbital']
    assert orb and orb[0].params['type'] == 'p' and (orb[0].x1, orb[0].y1) == (c.doc.atoms[atom].x,
                                                                                c.doc.atoms[atom].y)
    win.set_tool('tlc')
    c.click(None, 20, 20)
    assert len(c.doc.tlcs) == 1
    ethyl = next(i for i, b in c.doc.bonds.items() if len(c.doc.neighbors(b.a)) == 1 or len(c.doc.neighbors(b.b)) == 1)
    win.set_tool('massfrag')
    c.snapshot()
    assert c._cut('massfrag', [_cut_line(c.doc, ethyl)[:2], _cut_line(c.doc, ethyl)[2:]])
    assert any('91.0548' in t.text for t in c.doc.texts.values())
    for name in ('circle:shaded', 'orb:dz2:filled', 'massfrag', 'tlc', 'arc:120:dashed'):
        assert name in win.tool_actions


def test_nmr_drag_zooms_ppm_only(qapp):
    from types import SimpleNamespace as NS
    from nmr_dialog import SpectrumView
    v = SpectrumView()
    v.ax.set_xlim(10, 0); v.ax.set_ylim(0, 1)
    px = v.ax.transData.transform
    (x0, y0), (x1, y1) = px((8, 0.8)), px((6, 0.4))
    ev = lambda x, y, xd, yd: NS(inaxes=v.ax, x=x, y=y, xdata=xd, ydata=yd, button=1, dblclick=False, key=None)
    v._down(ev(x0, y0, 8, 0.8))
    v._move(ev(x1, y1, 6, 0.4))
    v._up(ev(x1, y1, 6, 0.4))
    assert v.ax.get_xlim() == pytest.approx((8, 6))
    assert v.ax.get_ylim() == pytest.approx((0, 1))


def test_pages_up_to_a0(win, qapp):
    assert render.page_size({'page': 'A0'}) == pytest.approx((2383.94, 3370.39))
    assert list(render.PAGES)[:5] == ['A0', 'A1', 'A2', 'A3', 'A4']
    c = win.canvas
    d = Doc(); d.style = {**d.style, 'page': 'A0'}
    c.set_doc(d)
    c.show_page()
    W, H = c.page_size()
    assert H * c.px_per_pt() <= c.height() - 40 + 1e-6                             # the whole A0 page fits
    assert win._printer().pageLayout().pageSize().id() == QtGui.QPageSize.PageSizeId.A0


def test_clean_up_has_a_second_key(win):
    assert [k.toString() for k in win.a_clean.shortcuts()] == ['Ctrl+Shift+K', 'Ctrl+Alt+K']


def _press(c, kind, pos, shift=False):
    from PyQt6.QtCore import QEvent, QPointF, Qt
    from PyQt6.QtGui import QMouseEvent
    mods = Qt.KeyboardModifier.ShiftModifier if shift else Qt.KeyboardModifier.NoModifier
    btn = Qt.MouseButton.LeftButton
    held = Qt.MouseButton.NoButton if kind == QEvent.Type.MouseButtonRelease else btn
    e = QMouseEvent(kind, QPointF(pos), QPointF(pos), QPointF(pos), btn, held, mods)
    {QEvent.Type.MouseButtonPress: c.mousePressEvent, QEvent.Type.MouseMove: c.mouseMoveEvent,
     QEvent.Type.MouseButtonRelease: c.mouseReleaseEvent}[kind](e)


def _drag(c, a, b, shift=False):
    from PyQt6.QtCore import QEvent, QPointF
    _press(c, QEvent.Type.MouseButtonPress, a, shift)
    _press(c, QEvent.Type.MouseMove, (a + b) / 2, shift)
    _press(c, QEvent.Type.MouseMove, b, shift)
    _press(c, QEvent.Type.MouseButtonRelease, b, shift)


def test_drag_inside_the_selection_box_moves_it(win, qapp):
    from PyQt6.QtCore import QPointF
    c = win.canvas
    c.set_doc(_mol('c1ccccc1'))
    win.set_tool('select')
    c.select_all()
    xs = [a.x for a in c.doc.atoms.values()]; ys = [a.y for a in c.doc.atoms.values()]
    middle = c.to_px((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    assert c.item_at(middle) is None                         # the empty middle of the ring
    x0 = xs[0]
    n = len(c.sel)
    _drag(c, middle, middle + QPointF(60, 0))
    assert c.doc.atoms[next(iter(c.doc.atoms))].x > x0 + 0.3
    assert len(c.sel) == n                                   # still selected, and one undo step
    c.undo()
    assert c.doc.atoms[next(iter(c.doc.atoms))].x == pytest.approx(x0)


def test_click_outside_the_box_still_starts_a_new_selection(win, qapp):
    from PyQt6.QtCore import QPointF
    c = win.canvas
    c.set_doc(_mol('c1ccccc1'))
    win.set_tool('select')
    c.select_all()
    far = c.to_px(max(a.x for a in c.doc.atoms.values()) + 3, min(a.y for a in c.doc.atoms.values()) - 3)
    x0 = c.doc.atoms[next(iter(c.doc.atoms))].x
    _drag(c, far, far + QPointF(5, 5))
    assert not c.sel
    assert c.doc.atoms[next(iter(c.doc.atoms))].x == x0


def test_shift_drag_inside_the_box_adds_with_a_band(win, qapp):
    from PyQt6.QtCore import QPointF
    c = win.canvas
    c.set_doc(_mol('c1ccccc1'))
    win.set_tool('select')
    c.select_all()
    xs = [a.x for a in c.doc.atoms.values()]; ys = [a.y for a in c.doc.atoms.values()]
    middle = c.to_px((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    x0 = c.doc.atoms[next(iter(c.doc.atoms))].x
    _drag(c, middle, middle + QPointF(30, 0), shift=True)
    assert c.doc.atoms[next(iter(c.doc.atoms))].x == x0      # Shift keeps its meaning: no move
