"""Charges and deleting as in ChemDraw, style sheets, page size and orientation, printing, NMR prediction (CSP5)."""
import os

import pytest
from PyQt6 import QtCore, QtWidgets
from PyQt6.QtTest import QTest

import chem
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
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles)); return d


def ethylamine():
    d = Doc()
    a, b, n = d.add_atom(0, 0), d.add_atom(0.87, -0.5), d.add_atom(1.73, 0, 'N')
    d.add_bond(a, b); d.add_bond(b, n)
    return d, b, n


# ------------------------------------------------------------------------------------------------ charges (octet)
def test_charges_follow_chemdraw():
    d, b, n = ethylamine()
    counts = []
    for q in (0, 1, 2, -1):
        d.atoms[n].charge = q
        counts.append(d.implicit_h(n))
    assert counts == [2, 3, 2, 1]                           # NH2, NH3+, NH2 2+, NH-
    d.atoms[n].charge = 0
    d.atoms[b].charge = 1
    assert d.implicit_h(b) == 1                             # CH+ (carbocation), CH- likewise
    d.atoms[b].charge = -1
    assert d.implicit_h(b) == 1
    assert chem.analyse(d).formula == 'C2H6N'


def test_charged_carbon_shows_its_label():
    d, b, _n = ethylamine()
    d.atoms[b].charge = 1
    lab = render.Renderer(d).labels[b]
    assert lab and [p[0] for p in lab['pieces']][:2] == ['C', 'H'] and lab['pieces'][-1][0] == '+'


def test_charge_stays_top_right():
    d = _doc('C[N+](=O)[O-]')
    r = render.Renderer(d)
    o = next(a for a, at in d.atoms.items() if at.charge == -1)
    pieces = r.labels[o]['pieces']
    assert pieces[-1][0] == '−' and pieces[-1][1] > pieces[0][1]   # the minus after the O, not before it


# ------------------------------------------------------------------------------------------------ deleting
def test_delete_over_heteroatom_turns_it_to_carbon(win):
    c = win.canvas
    d, b, n = ethylamine()
    d.atoms[n].charge = 1
    c.set_doc(d)
    c.hover = ('atom', n)
    QTest.keyClick(c, QtCore.Qt.Key.Key_Delete)
    at = c.doc.atoms[n]
    assert (at.el, at.charge) == ('C', 0) and len(c.doc.bonds) == 2   # bonds kept
    c.hover = ('atom', n)
    QTest.keyClick(c, QtCore.Qt.Key.Key_Delete)
    assert n not in c.doc.atoms                             # a plain carbon goes, with its bond


def test_delete_nickname_and_eraser(win):
    c = win.canvas
    d, b, n = ethylamine()
    d.set_nickname(n, 'OMe')
    c.set_doc(d)
    win.set_tool('erase')
    QTest.mouseClick(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(d.atoms[n].x, d.atoms[n].y).toPoint())
    assert c.doc.atoms[n].label == '' and c.doc.atoms[n].el == 'C'
    c.undo()
    assert c.doc.atoms[n].label == 'OMe'


# ------------------------------------------------------------------------------------------------ styles, page, print
def test_style_sheets_apply(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    win.apply_style_sheet('RSC Document')
    s = render.style_of(c.doc)
    assert s['bond_length'] == pytest.approx(12.2) and s['font_size'] == 7 and s['line_width'] == pytest.approx(0.45)
    c.doc.style['show_stereo'] = True
    win.apply_style_sheet('New Slide')
    s = render.style_of(c.doc)
    assert s['bond_length'] == 30 and s['caption_size'] == 16 and s['show_stereo']   # display options kept
    assert len(render.STYLE_SHEETS) >= 12
    for name, sheet in render.STYLE_SHEETS.items():
        assert set(sheet) <= set(render.DEFAULT_STYLE), name


def test_new_from_style_sheet(win):
    win.new_from_style_sheet('Wiley Document')
    assert render.style_of(win.canvas.doc)['bond_length'] == 17.0


def test_page_size_and_orientation():
    assert render.page_size({'page': 'A4', 'orientation': 'portrait'}) == pytest.approx((595.28, 841.89))
    assert render.page_size({'page': 'A4', 'orientation': 'landscape'}) == pytest.approx((841.89, 595.28))
    assert render.page_size({'page': 'Letter'}) == (612.0, 792.0)


def test_document_settings_keep_display_options(win, monkeypatch):
    import kekulex as se
    c = win.canvas
    c.doc.style['show_atom_numbers'] = True
    monkeypatch.setattr(se.StyleDialog, 'exec', lambda self: (self.orientation.setCurrentText('Landscape'),
                                                              self.page.setCurrentText('A3'), True)[-1])
    win.style_dialog()
    s = render.style_of(c.doc)
    assert s['page'] == 'A3' and s['orientation'] == 'landscape' and s['show_atom_numbers']


def test_print_to_pdf(win, tmp_path):
    from PyQt6.QtPrintSupport import QPrinter
    win.canvas.set_doc(_doc('c1ccccc1O'))
    printer = win._printer()
    printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
    out = str(tmp_path / 'print.pdf')
    printer.setOutputFileName(out)
    win._print_on(printer)
    data = open(out, 'rb').read()
    assert data.startswith(b'%PDF') and len(data) > 1000


# ------------------------------------------------------------------------------------------------ NMR (CSP5)
nmr = pytest.importorskip('nmr')
needs_csp5 = pytest.mark.skipif(not nmr.available(), reason='CSP5 is not installed')


@needs_csp5
def test_nmr_1h_ethanol():
    r = nmr.predict(_doc('CCO'), None, '1H', 'cdcl3')
    sig = {s.n_h: s for s in r.signals}
    assert sig[3].multiplicity == 't' and 1.0 < sig[3].shift < 1.5             # CH3
    assert sig[2].multiplicity == 'q' and 3.4 < sig[2].shift < 4.0             # CH2
    assert sig[1].multiplicity == 'br s'                                       # OH
    assert nmr.journal_text(r).startswith('¹H NMR (400 MHz, CDCl₃, predicted) δ ')


@needs_csp5
def test_nmr_13c_symmetry_and_atoms():
    d = _doc('O=[N+]([O-])c1ccccc1')
    r = nmr.predict(d, None, '13C', 'cdcl3')
    assert len(r.signals) == 4 and sorted(len(s.atoms) for s in r.signals) == [1, 1, 2, 2]
    assert 145 < r.signals[0].shift < 152                                       # ipso C first (highest)
    carbons = {a for a, at in d.atoms.items() if at.el == 'C'}
    assert set().union(*[set(s.atoms) for s in r.signals]) == carbons


@needs_csp5
def test_nmr_counts_nicknames():
    d = Doc()
    a, b = d.add_atom(0, 0), d.add_atom(1, 0)
    d.add_bond(a, b)
    d.set_nickname(b, 'OMe')
    r = nmr.predict(d, None, '1H', 'cdcl3')
    assert [(s.n_h, s.multiplicity) for s in r.signals] == [(6, 's')]           # dimethyl ether


@needs_csp5
def test_nmr_dialog(win, qapp):
    from nmr_dialog import NMRDialog
    c = win.canvas
    c.set_doc(_doc('CC(=O)O'))
    d = NMRDialog(c, None, '13C', win)
    d.worker.wait(120000)
    qapp.processEvents()
    assert d.table.rowCount() == 2 and '¹³C NMR' in d.text.text()
    d.table.setCurrentCell(0, 0)
    assert c.sel
    d.label_atoms()
    assert len(c.doc.texts) == 2
    d.close()


def test_couplings_known_patterns():
    """Expected J without CSP5: chains, aromatic ortho / meta, trans alkene, aldehyde, a pyranose in its chair."""
    from rdkit import Chem

    def pattern(smiles, heavy):
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
        h = next(n.GetIdx() for n in m.GetAtomWithIdx(heavy).GetNeighbors() if n.GetSymbol() == 'H')
        c = nmr.couplings(m, h, ranks, nmr._conformer(m))
        return nmr.multiplet_name(c), c
    assert pattern('CCCO', 1) == ('sext', [(7.0, 5)])
    assert pattern('COc1ccc(Cl)cc1Cl', 4) == ('dd', [(8.0, 1), (2.0, 1)])         # H5 of 2,4-dichloroanisole
    name, c = pattern('O=C/C=C/c1ccccc1', 2)
    assert name == 'dd' and c[0] == (16.0, 1) and c[1] == (2.5, 1)
    name, c = pattern('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O', 4)         # beta-glucose H1: diaxial
    assert name == 'd' and c[0][0] > 7.5
    s = nmr.Signal(7.0, 0, [1], 1, [(8.0, 1), (2.0, 1)], 'dd')
    assert s.j_text() == '8.0, 2.0' and len(nmr.lines(s)) == 4


@needs_csp5
def test_nmr_dialog_betanin_reopen_and_zoom(win, qapp):
    """The structure from the user's report (betanin, zwitterion): the window opens twice in a row (the second used
    to raise 'wrapped C/C++ object of type NMRDialog has been deleted'), J goes into the table and the text, and the
    spectrum zooms with the wheel."""
    c = win.canvas
    c.set_doc(_doc('O=C([O-])[C@@H]1Cc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc2[N+]1=CC=C1C[C@@H]'
                   '(C(=O)O)NC(C(=O)[O-])=C1'))
    for _ in range(2):
        win.open_nmr('1H')
        d = [w for w in win.findChildren(QtWidgets.QDialog) if type(w).__name__ == 'NMRDialog'][-1]
        d.worker.wait(180000)
        qapp.processEvents()
        assert d.table.rowCount() > 10
        js = [d.table.item(r, 3).text() for r in range(d.table.rowCount())]
        assert any(js) and 'J = ' in d.text.text()
        lo, hi = d.view.ax.get_xlim()
        from types import SimpleNamespace
        d.view._wheel(SimpleNamespace(inaxes=d.view.ax, xdata=5.0, ydata=0.0, x=0, y=0, button='up', key=None))
        lo2, hi2 = d.view.ax.get_xlim()
        assert abs(hi2 - lo2) < abs(hi - lo)
        d.close()
        qapp.processEvents()
        QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete.value)


# ------------------------------------------------------------------------------------------------ moving charges, bond counts
def test_drag_charge_sign(win):
    c = win.canvas
    d, b, n = ethylamine()
    d.atoms[n].charge = 1
    c.set_doc(d)
    win.set_tool('select')
    r = c.rect_px(c.renderer().charge_rect(n))
    start = r.center()
    end = c.to_px(d.atoms[n].x, d.atoms[n].y - 0.6)          # drop it straight above the N
    QTest.mousePress(c, QtCore.Qt.MouseButton.LeftButton, pos=start.toPoint())
    QTest.mouseMove(c, ((start + end) / 2).toPoint())
    QTest.mouseMove(c, end.toPoint())
    QTest.mouseRelease(c, QtCore.Qt.MouseButton.LeftButton, pos=end.toPoint())
    at = c.doc.atoms[n]
    assert at.charge_pos and abs(at.charge_pos[0]) < 0.1 and at.charge_pos[1] == pytest.approx(-0.6, abs=0.1)
    assert (at.x, at.y) == (d.atoms[n].x, d.atoms[n].y)      # the atom itself did not move
    moved = c.renderer().charge_rect(n).center()
    assert abs(moved.x() - at.x * c.renderer().BL) < 1.0
    c.select_all(); c.transform_selection('flip_v')
    assert c.doc.atoms[n].charge_pos                         # it follows the atom
    c.set_selection([('atom', n)]); c.reset_charge_positions()
    assert c.doc.atoms[n].charge_pos is None
    c.undo()
    assert c.doc.atoms[n].charge_pos


def test_bond_count_problems():
    d = Doc()
    c = d.add_atom(0, 0)
    for k in range(5):
        d.add_bond(c, d.add_atom(k, 1))
    assert d.valence_problem(c) == 'C takes at most 4 bonds, but has 5.'
    n = d.add_atom(9, 9, 'N')
    for k in range(4):
        d.add_bond(n, d.add_atom(9 + k, 10))
    assert 'N takes at most 3 bonds, but has 4' in d.valence_problem(n)
    d.atoms[n].charge = 1
    assert d.valence_problem(n) == ''                        # NR4+ is fine
    d.atoms[n].charge = 2
    assert 'N +2 takes at most 3' in d.valence_problem(n)
    o = d.add_atom(20, 0)
    d.add_bond(o, d.add_atom(21, 0)); d.add_bond(o, d.add_atom(19, 0))
    d.set_nickname(o, 'OMe')
    assert d.valence_problem(o) == 'OMe attaches by one bond, but has 2.'
    assert any('OMe attaches' in e for e in chem.analyse(d).errors)


def test_problem_tooltip(win):
    c = win.canvas
    d = Doc()
    a = d.add_atom(5, 5)
    for k in range(5):
        d.add_bond(a, d.add_atom(5 + k * 0.3, 6))
    c.set_doc(d)
    msgs = []
    c.message.connect(msgs.append)
    QTest.mouseMove(c, c.to_px(5, 5).toPoint())
    assert msgs and 'at most 4 bonds' in msgs[-1]
