"""The ChemDraw-like features: nicknames, CDXML, isotope patterns, bond / arrow / shape tools, typing in place,
arranging, templates."""
import math
import os

import pytest
from PyQt6 import QtCore
from PyQt6.QtTest import QTest

import cdxml
import chem
import model
import templates
from model import Doc

SAMPLE = os.path.join(os.path.expanduser('~'), 'Desktop', 'forclaude.cdxml')


@pytest.fixture
def win(qapp):
    import structure_editor as se
    w = se.MainWindow()
    w.resize(1300, 850)
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


def ethyl_with(label):
    d = Doc()
    a, b = d.add_atom(0, 0), d.add_atom(1, 0)
    d.add_bond(a, b)
    d.set_nickname(b, label)
    return d, b


# ------------------------------------------------------------------------------------------------ nicknames
@pytest.mark.parametrize('label, smiles', [('OMe', 'COC'), ('CO2H', 'CC(=O)O'), ('Boc', 'CC(=O)OC(C)(C)C'),
                                           ('NO2', 'C[N+](=O)[O-]'), ('Ph', 'Cc1ccccc1'), ('OTs', 'COS(=O)(=O)c1ccc(C)cc1'),
                                           ('TMS', 'C[Si](C)(C)C')])
def test_nickname_counts_in_full(label, smiles):
    d, _b = ethyl_with(label)
    assert chem.analyse(d).smiles == chem.analyse(_doc(smiles)).smiles
    assert len(d.atoms) == 2                               # still drawn as one label


def _doc(smiles):
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles)); return d


def test_parse_label():
    assert model.parse_label('OMe') == ('nickname', 'OMe')
    assert model.parse_label('ome') == ('nickname', 'OMe')
    assert model.parse_label('NH2') == ('element', 'N', 0)
    assert model.parse_label('O-') == ('element', 'O', -1)
    assert model.parse_label('NH3+') == ('element', 'N', 1)
    assert model.parse_label('Na+') == ('element', 'Na', 1)
    assert model.parse_label('Xyz') is None


def test_label_turns_when_bonded_on_the_right():
    assert model.label_order('CO2Me', True)[0] == ['Me', 'O2', 'C']
    assert model.label_order('NO2', True)[0] == ['O2', 'N']
    assert ''.join(model.label_order('TMS', True)[0]) == 'TMS'   # not chemical tokens: stays as typed


def test_expand_label(win):
    c = win.canvas
    d, b = ethyl_with('CO2Me')
    c.set_doc(d)
    c.expand_labels()
    assert len(c.doc.atoms) == 5 and not any(a.label for a in c.doc.atoms.values())
    assert chem.analyse(c.doc).smiles == 'COC(C)=O'                # methyl acetate, now drawn out


def test_type_label_in_place(win, qapp):
    c = win.canvas
    win.set_tool('bond1'); click(c, 3, 3)
    end = max(c.doc.atoms)
    c.edit_label(end)
    w = c.editor[0]
    w.setText('CO2H')
    QTest.keyClick(w, QtCore.Qt.Key.Key_Return)
    assert c.doc.atoms[end].label == 'CO2H'
    assert chem.analyse(c.doc).formula == 'C2H4O2'


def test_type_caption_in_place(win):
    c = win.canvas
    win.set_tool('text'); click(c, 4, 4)
    w = c.editor[0]
    w.setPlainText('H2SO4\ncat.')
    c.set_text_props(formula=True)
    QTest.keyClick(w, QtCore.Qt.Key.Key_Escape)
    t = next(iter(c.doc.texts.values()))
    assert t.text == 'H2SO4\ncat.' and t.formula


# ------------------------------------------------------------------------------------------------ CDXML
@pytest.mark.skipif(not os.path.exists(SAMPLE), reason='the sample ChemDraw file is not on this computer')
def test_read_users_chemdraw_file(qapp):
    doc, notes = cdxml.read(SAMPLE)
    formulas = sorted(chem.analyse(doc, f).formula for f in doc.fragments())
    assert formulas == ['C6H5NO2', 'C6H6', 'H2O4S', 'HNO3']
    assert {t.text for t in doc.texts.values()} >= {'nitric acid', 'sulfuric acid', 'nitrobenzene', '+'}
    assert len(doc.arrows) == 1 and 'the stoichiometry table' in notes
    assert doc.style['bond_length'] == pytest.approx(14.4)


def test_cdxml_round_trip(qapp, tmp_path):
    from rdkit import Chem
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('C[C@@H](N)C(=O)O'), (3, 3))
    a = d.add_atom(8, 3); b = d.add_atom(8.87, 2.5); d.add_bond(a, b); d.set_nickname(b, 'OMe')
    i = d.add_bond(d.add_atom(3, 6), d.add_atom(4, 6)); d.bonds[i].style = 'wavy'
    d.add_arrow(1, 8, 4, 8, 'equilibrium'); d.add_arrow(5, 8, 7, 8, 'curved', 0.3)
    d.add_text(2, 9, 'H2O, 25 °C', align='center', formula=True, bold=True)
    d.add_shape(0, 0, 6, 5, 'brackets')
    d.texts[max(d.texts)].color = '#E01E1E'
    p = str(tmp_path / 'x.cdxml')
    cdxml.write(d, p)
    rd = sorted(Chem.MolToSmiles(m) for m in Chem.MolsFromCDXMLFile(p))   # RDKit's own CDXML reader agrees
    assert rd == sorted([chem.analyse(d, f).smiles for f in d.fragments()])
    back, notes = cdxml.read(p)
    assert sorted(chem.analyse(back, f).smiles for f in back.fragments()) == rd
    assert any(at.label == 'OMe' for at in back.atoms.values())
    assert sorted(a.kind for a in back.arrows.values()) == ['curved', 'equilibrium']
    t = next(iter(back.texts.values()))
    assert (t.text, t.align, t.formula, t.bold, t.color) == ('H2O, 25 °C', 'center', True, True, '#E01E1E')
    assert next(iter(back.shapes.values())).kind == 'brackets'
    assert any(b.style == 'wavy' for b in back.bonds.values())
    assert not notes


def test_open_and_save_cdxml_in_window(win, tmp_path):
    d = _doc('c1ccccc1O')
    p = str(tmp_path / 'phenol.cdxml')
    cdxml.write(d, p)
    win.open(p)
    assert chem.analyse(win.canvas.doc).smiles == 'Oc1ccccc1'
    assert win.save() and open(p, encoding='utf-8').read().startswith('<?xml')


# ------------------------------------------------------------------------------------------------ masses
def test_isotope_pattern():
    peaks = chem.isotope_pattern({'C': 8, 'H': 11, 'N': 4, 'O': 2}, 1)       # caffeine [M+H]+
    assert peaks[0][0] == pytest.approx(195.087652, abs=2e-6) and peaks[0][1] == 100
    assert peaks[1][1] == pytest.approx(10.3, abs=0.4)
    cl = chem.isotope_pattern({'C': 1, 'H': 3, 'Cl': 1})                       # chlorine: M+2 about a third of M
    assert cl[1][1] == pytest.approx(32.0, abs=1.0) or cl[2][1] == pytest.approx(32.0, abs=1.0)


def test_elemental_analysis():
    info = chem.analyse(_doc('CC(=O)O'))
    assert info.analysis == 'C, 40.00; H, 6.71; O, 53.28'


def test_cip_labels():
    d = _doc('C[C@H](N)C(=O)O')                                                # L-alanine
    atoms, _bonds = chem.cip_labels(d)
    assert list(atoms.values()) == ['S']
    d = _doc('C/C=C/C')
    assert list(chem.cip_labels(d)[1].values()) == ['E']


def test_analysis_caption(win):
    c = win.canvas
    c.set_doc(_doc('CC(=O)O'))
    c.add_analysis()
    t = next(iter(c.doc.texts.values()))
    assert 'Chemical Formula: C2H4O2' in t.text and 'Exact Mass: 60.02' in t.text and t.formula


# ------------------------------------------------------------------------------------------------ tools
def test_arrow_and_shape_tools(win):
    c = win.canvas
    for k, tool in enumerate(('retro', 'resonance', 'curved', 'curved_half', 'arrow_dashed')):
        win.set_tool(tool); drag(c, 2, 2 + k, 5, 2 + k)
    assert sorted(a.kind for a in c.doc.arrows.values()) == sorted(['retro', 'resonance', 'curved', 'curved_half', 'dashed'])
    curved = next(i for i, a in c.doc.arrows.items() if a.kind == 'curved')
    bend = c.doc.arrows[curved].bend
    pts = model.arrow_points(c.doc.arrows[curved])
    mid = pts[len(pts) // 2]
    win.set_tool('curved'); click(c, *mid)
    assert c.doc.arrows[curved].bend == -bend                                  # clicked again: bends the other way
    win.set_tool('brackets'); drag(c, 8, 2, 11, 5)
    assert next(iter(c.doc.shapes.values())).kind == 'brackets'


def test_bond_style_tools_and_keys(win):
    c = win.canvas
    win.set_tool('bond1'); click(c, 3, 3)
    b = next(iter(c.doc.bonds))
    bd = c.doc.bonds[b]
    mid = ((c.doc.atoms[bd.a].x + c.doc.atoms[bd.b].x) / 2, (c.doc.atoms[bd.a].y + c.doc.atoms[bd.b].y) / 2)
    win.set_tool('bold'); click(c, *mid)
    assert bd.style == 'bold'
    c.hover = ('bond', b)
    QTest.keyClicks(c, 'y'); assert bd.style == 'wavy'
    QTest.keyClicks(c, '2'); assert bd.order == 2 and bd.style == ''
    QTest.keyClicks(c, 'd'); assert bd.order == 2 and bd.style == 'dashed'


def test_nickname_hotkey_keeps_label(win):
    c = win.canvas
    win.set_tool('bond1'); click(c, 3, 3)
    end = max(c.doc.atoms)
    c.hover = ('atom', end)
    QTest.keyClicks(c, 'O')                                                     # Shift+O = OMe, as in ChemDraw
    assert c.doc.atoms[end].label == 'OMe' and len(c.doc.atoms) == 2


def test_align_and_distribute(win):
    c = win.canvas
    for k, smi in enumerate(('C', 'CC', 'CCC')):
        c.doc.insert_rdkit(chem.mol_from_text(smi), (3 + 4 * k, 3 + k * 1.3))
    c.edited()
    c.select_all()
    c.align('top')
    rend = c.renderer()
    tops = [g[1].top() for g in c.object_groups()]
    assert max(tops) - min(tops) < 0.01
    c.distribute('h')
    rects = sorted((g[1] for g in c.object_groups()), key=lambda r: r.left())
    gaps = [rects[k + 1].left() - rects[k].right() for k in range(2)]
    assert abs(gaps[0] - gaps[1]) < 0.01
    del rend


def test_color_and_rotate(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    c.select_all()
    c.set_color('#E01E1E')
    assert all(a.color == '#E01E1E' for a in c.doc.atoms.values())
    before = {i: (a.x, a.y) for i, a in c.doc.atoms.items()}
    c.rotate(90)
    after = {i: (a.x, a.y) for i, a in c.doc.atoms.items()}
    d0 = [math.dist(before[i], before[j]) for i in before for j in before]
    d1 = [math.dist(after[i], after[j]) for i in after for j in after]
    assert d0 == pytest.approx(d1)
    assert chem.analyse(c.doc).smiles == 'CCO'


def test_templates_are_valid():
    for cat, items in templates.TEMPLATES.items():
        for name, smiles in items:
            assert chem.mol_from_text(smiles) is not None, (cat, name)


def test_template_menu_inserts(win):
    win.insert_text(templates.TEMPLATES['Amino acids (L)'][1][1])
    assert chem.analyse(win.canvas.doc).smiles == 'C[C@H](N)C(=O)O'
