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
    import kekulex as se
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
    assert len(doc.arrows) == 1 and not notes                  # the stoichiometry grid comes over too
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
    names = dict(templates.TEMPLATES['Reagents'])
    win.insert_text(names['DMAP'])
    assert chem.analyse(win.canvas.doc).formula == 'C7H10N2'
    assert not any('Betalain' in c for c in templates.TEMPLATES)     # bio templates are back (tests/test_bio.py)


# ------------------------------------------------------------------------------------------------ stoichiometry
def _nitration(c):
    """Benzene + HNO3 -> nitrobenzene, H2SO4 under the arrow, drawn through the canvas."""
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('c1ccccc1'), (2, 3))
    doc.insert_rdkit(chem.mol_from_text('O=[N+]([O-])O'), (6, 3))
    doc.insert_rdkit(chem.mol_from_text('O=S(=O)(O)O'), (11, 5))
    doc.add_arrow(9, 3, 13, 3)
    doc.insert_rdkit(chem.mol_from_text('O=[N+]([O-])c1ccccc1'), (17, 3))
    c.set_doc(doc)


def test_reaction_columns(win):
    import stoich
    _nitration(win.canvas)
    reac, prod = stoich.reaction_columns(win.canvas.doc)
    doc = win.canvas.doc
    assert [chem.analyse(doc, f).formula for f in reac] == ['C6H6', 'HNO3', 'H2O4S']
    assert [chem.analyse(doc, f).formula for f in prod] == ['C6H5NO2']


def test_parse_units():
    import stoich
    assert stoich.parse('250 mg', 'mass') == pytest.approx(0.25)
    assert stoich.parse('1,5mmol', 'moles') == pytest.approx(0.0015)
    assert stoich.parse('2 mL', 'volume') == pytest.approx(2.0)
    assert stoich.parse('500 µL', 'volume') == pytest.approx(0.5)
    assert stoich.parse('0.5 M', 'molarity') == pytest.approx(0.5)
    assert stoich.parse('20 mM', 'molarity') == pytest.approx(0.02)
    assert stoich.parse('1.2 g/mL', 'density') == pytest.approx(1.2)
    assert stoich.parse('95 %', 'percent') == pytest.approx(95)
    assert stoich.parse('', 'mass') is None
    with pytest.raises(ValueError):
        stoich.parse('ten grams', 'mass')
    with pytest.raises(ValueError):
        stoich.parse('5 furlongs', 'mass')


def test_stoichiometry_calculations(win):
    import stoich
    c = win.canvas
    _nitration(c)
    gid = c.analyze_stoichiometry()
    g = c.doc.grids[gid]
    g.reactants[0]['mass'] = 7.8114                         # 0.1 mol benzene (limiting)
    g.reactants[1]['equiv'] = 1.2
    g.reactants[2].update(volume=1.0, density=1.83, pct=96.0)  # H2SO4 typed as a volume of 96 % acid
    g.products[0]['measured'] = 10.0
    t = stoich.compute(c.doc, g)
    r, p = t['reactants'], t['products'][0]
    assert r[0].cells['moles'].value == pytest.approx(0.1, rel=1e-4)
    assert r[1].cells['moles'].value == pytest.approx(0.12, rel=1e-4)
    assert r[1].cells['mass'].value == pytest.approx(0.12 * 63.012, rel=1e-3)
    h2so4 = 1.0 * 1.83 * 0.96 / 98.079
    assert r[2].cells['moles'].value == pytest.approx(h2so4, rel=1e-4)
    assert r[2].cells['equiv'].value == pytest.approx(h2so4 / 0.1, rel=1e-4)   # its equivalents follow
    assert p.cells['exp_mass'].value == pytest.approx(12.3111, rel=1e-3)
    assert p.cells['yield'].value == pytest.approx(10.0 / 12.3111 * 100, rel=1e-3)
    assert r[0].cells['mass'].typed and not r[1].cells['mass'].typed
    g.reactants[0]['limiting'], g.reactants[1]['limiting'] = False, True      # HNO3 limiting at 1 eq
    g.reactants[1].pop('equiv'); g.reactants[1]['moles'] = 0.05
    t = stoich.compute(c.doc, g)
    assert t['reactants'][0].cells['equiv'].value == pytest.approx(2.0, rel=1e-4)
    assert t['products'][0].cells['exp_moles'].value == pytest.approx(0.05)


def test_grid_draws_and_moves(win):
    c = win.canvas
    _nitration(c)
    gid = c.analyze_stoichiometry()
    rect = c.renderer().grid_rect(gid)
    assert rect.width() > 200 and rect.height() > 100
    assert c.item_at(c.rect_px(rect).center()) == ('grid', gid)
    c.set_selection([('grid', gid)])
    c.transform_selection('flip_h')                         # grids move with the selection
    assert c.doc.grids[gid].x is not None
    assert c.doc.grids[gid] is not None and render_ok(c)


def render_ok(c):
    import render
    img = render.to_image(c.doc, 150)
    return img.width() > 100


def test_grid_dialog_typing(win, qapp):
    from stoich_dialog import GridDialog
    c = win.canvas
    _nitration(c)
    gid = c.analyze_stoichiometry()
    d = GridDialog(c, gid, win)
    t = d.tables['reactants']
    row = t.rows.index('mass')
    t.item(row, 0).setText('781 mg')
    assert c.doc.grids[gid].reactants[0]['mass'] == pytest.approx(0.781)
    assert d.tables['reactants'].item(t.rows.index('moles'), 0).text().startswith('10,00') or \
        d.tables['reactants'].item(t.rows.index('moles'), 0).text().startswith('10.00')
    t.item(row, 0).setText('lots')
    assert 'not a number' in d.error.text()
    d.accept()
    assert c.doc.grids[gid].reactants[0]['mass'] == pytest.approx(0.781)
    c.undo()
    assert 'mass' not in c.doc.grids[gid].reactants[0]


def test_grid_survives_save_and_copy(win, tmp_path):
    import json
    c = win.canvas
    _nitration(c)
    gid = c.analyze_stoichiometry()
    c.doc.grids[gid].reactants[0]['mass'] = 1.0
    back = Doc.from_dict(json.loads(json.dumps(c.doc.to_dict())))
    assert back.grids[gid].reactants[0]['mass'] == 1.0
    c.select_all()
    win.copy(); win.paste()
    assert len(c.doc.grids) == 2
    new = [g for i, g in c.doc.grids.items() if i != gid][0]
    assert set(new.reactants[0]['atoms']).isdisjoint(c.doc.grids[gid].reactants[0]['atoms'])


@pytest.mark.skipif(not os.path.exists(SAMPLE) or 'stoichiometrygrid' not in open(SAMPLE, encoding='utf-8').read(),
                    reason='the sample ChemDraw file is missing or has no stoichiometry grid')
def test_grid_from_chemdraw_file(qapp):
    import stoich
    doc, notes = cdxml.read(SAMPLE)
    assert len(doc.grids) == 1 and not notes
    g = next(iter(doc.grids.values()))
    t = stoich.compute(doc, g)
    assert [c.formula for c in t['reactants']] == ['C6H6', 'HNO3', 'H2O4S']
    assert [c.formula for c in t['products']] == ['C6H5NO2']
    assert t['reactants'][0].cells['limiting'].value


def test_check_structure(win):
    c = win.canvas
    d = Doc()
    a = d.add_atom(0, 0)
    for k in range(5):
        d.add_bond(a, d.add_atom(math.cos(k), math.sin(k)))
    c.set_doc(d)
    assert c.check_structure()
    c.set_doc(_doc('CCO'))
    assert c.check_structure() == []


def test_properties():
    rows = dict((n, v) for n, v, _t in chem.properties(_doc('CC(=O)Oc1ccccc1C(=O)O'), None))   # aspirin
    assert rows['LogP'] == '1.31' and rows['H-bond donors'] == '1' and rows['Rule of five'] == 'passes'


def test_paste_buttons(win):
    c = win.canvas
    c.set_doc(_doc('CC(=O)O'))
    win.info.refresh()
    win.info.paste('exact')
    win.info.paste('formula')
    texts = sorted(t.text for t in c.doc.texts.values())
    assert texts == ['Chemical Formula: C2H4O2', 'Exact Mass: 60.02']
    ys = sorted(t.y for t in c.doc.texts.values())
    assert ys[1] > ys[0]                                    # the second goes under the first


def test_pubchem_lookups(win, monkeypatch):
    import pubchem
    monkeypatch.setattr(pubchem, 'smiles_for_name', lambda name: 'c1ccccc1O' if name == 'phenol' else None)
    monkeypatch.setattr(pubchem, 'name_for_inchikey', lambda key: ('phenol', 'phenol'))
    from PyQt6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QInputDialog, 'getText', staticmethod(lambda *a, **k: ('phenol', True)))
    win.name_to_structure()
    assert chem.analyse(win.canvas.doc).smiles == 'Oc1ccccc1'
    win.structure_to_name()
    assert any(t.text == 'phenol' for t in win.canvas.doc.texts.values())


# ------------------------------------------------------------------------------------------------ lasso, colours, menu
def test_lasso_selects_inside_outline(win):
    c = win.canvas
    c.doc.insert_rdkit(chem.mol_from_text('c1ccccc1'), (4, 4))
    c.doc.insert_rdkit(chem.mol_from_text('CCO'), (12, 4))
    c.edited()
    win.set_tool('lasso')
    pts = [(1.5, 1.5), (7, 1.0), (7.5, 7), (1.0, 7.5), (1.5, 1.5)]   # around the benzene only
    QTest.mousePress(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(*pts[0]).toPoint())
    for x, y in pts[1:]:
        QTest.mouseMove(c, c.to_px(x, y).toPoint())
    QTest.mouseRelease(c, QtCore.Qt.MouseButton.LeftButton, pos=c.to_px(*pts[-1]).toPoint())
    assert len(c.selected_atoms()) == 6
    assert chem.analyse(c.doc, c.selected_atoms()).formula == 'C6H6'


def test_highlight_and_ring_fill(win, tmp_path):
    import json
    c = win.canvas
    c.set_doc(_doc('c1ccc2ccccc2c1'))                       # naphthalene: two rings
    c.select_all()
    c.set_highlight('#2EE62E')
    assert all(a.highlight == '#2EE62E' for a in c.doc.atoms.values())
    assert all(b.highlight == '#2EE62E' for b in c.doc.bonds.values())
    c.set_ring_fill('#FFF7B0')
    assert len(c.doc.fills) == 2 and all(len(f.atoms) == 6 for f in c.doc.fills.values())
    c.set_ring_fill('#C9DDF7')                              # recolouring does not add more fills
    assert len(c.doc.fills) == 2 and {f.color for f in c.doc.fills.values()} == {'#C9DDF7'}
    back = Doc.from_dict(json.loads(json.dumps(c.doc.to_dict())))
    assert len(back.fills) == 2 and all(a.highlight for a in back.atoms.values())
    assert render_ok(c)
    c.set_ring_fill('')
    c.set_highlight('')
    assert not c.doc.fills and not any(a.highlight for a in c.doc.atoms.values())
    c.undo()
    assert not c.doc.fills                                   # undo brings back the highlight, not the fill
    c.set_ring_fill('#FFF7B0')
    victim = c.doc.fills[next(iter(c.doc.fills))].atoms[0]
    c.doc.remove_atom(victim)
    assert len(c.doc.fills) < 2                             # a fill goes with its ring


def test_flip_mirrors_rotate_180_keeps(win):
    c = win.canvas
    c.set_doc(_doc('C[C@@H](N)C(=O)O'))
    before = chem.analyse(c.doc).smiles
    c.transform_selection('flip_h')                         # Rotate 180° Horizontal
    assert chem.analyse(c.doc).smiles == before
    c.transform_selection('mirror_h')                       # Flip Horizontal: the enantiomer
    assert chem.analyse(c.doc).smiles == before.replace('@@', '@')


def test_add_frame(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    c.select_all()
    c.add_frame('brackets')
    sh = next(iter(c.doc.shapes.values()))
    xs = [a.x for a in c.doc.atoms.values()]
    assert sh.kind == 'brackets' and min(sh.x1, sh.x2) < min(xs) and max(sh.x1, sh.x2) > max(xs)


def test_right_click_menu(win, qapp, monkeypatch):
    from PyQt6 import QtWidgets
    c = win.canvas
    c.set_doc(_doc('CC(=O)O'))
    shown = []
    monkeypatch.setattr(QtWidgets.QMenu, 'exec', lambda self, *a: shown.append(
        [act.text() for act in self.actions()]))
    atom = next(iter(c.doc.atoms))
    pos = c.to_px(c.doc.atoms[atom].x, c.doc.atoms[atom].y).toPoint()
    QTest.mouseClick(c, QtCore.Qt.MouseButton.RightButton, pos=pos)
    assert shown and {'&Align', 'Add &Frame', 'A&tom', '&Molecule', 'A&nalysis'} <= set(shown[-1])
    assert ('atom', atom) in c.sel                          # the clicked object is what the menu acts on
    c.set_selection([])
    QTest.mouseClick(c, QtCore.Qt.MouseButton.RightButton, pos=c.to_px(30, 30).toPoint())
    assert '&Paste' in shown[-1] and 'A&tom' not in shown[-1]


def test_right_drag_still_pans(win):
    c = win.canvas
    o = QtCore.QPointF(c.origin)
    QTest.mousePress(c, QtCore.Qt.MouseButton.RightButton, pos=QtCore.QPoint(300, 300))
    QTest.mouseMove(c, QtCore.QPoint(340, 330))
    QTest.mouseMove(c, QtCore.QPoint(380, 360))
    QTest.mouseRelease(c, QtCore.Qt.MouseButton.RightButton, pos=QtCore.QPoint(380, 360))
    assert c.origin != o


def test_edit_selected_from_menu(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    c.select_all()
    c.edit_selected('bond', order=2, stereo='', style='', side=None)
    assert all(b.order == 2 for b in c.doc.bonds.values())
    c.edit_selected('atom', charge=1)
    assert all(a.charge == 1 for a in c.doc.atoms.values())


# ------------------------------------------------------------------------------------------------ names (Openclatura)
@pytest.mark.parametrize('smiles, name', [
    ('CCO', 'ethanol'), ('OC(=O)c1ccccc1', 'benzoic acid'), ('O=[N+]([O-])c1ccccc1', 'nitrobenzene'),
    ('C[C@H](N)C(=O)O', '(2S)-2-aminopropanoic acid'), ('CN(C)c1ccncc1', 'N,N-dimethylpyridin-4-amine'),
    ('CC(C)Cc1ccc(cc1)C(C)C(=O)O', '2-[4-(2-methylpropyl)phenyl]propanoic acid'),      # ( ) inside [ ]
])
def test_iupac_name(smiles, name):
    assert chem.iupac_name(_doc(smiles)) == name


def test_iupac_name_counts_nicknames_and_stereo():
    d, _b = ethyl_with('CO2H')                              # CH3-CO2H drawn with a label
    assert chem.iupac_name(d) == 'acetic acid'
    assert chem.iupac_name(_doc('C/C=C/C')) == '(2E)-but-2-ene'


def test_nest_enclosures():
    assert chem.nest_enclosures('1-(2-(3-(4-methylphenyl)propyl)phenyl)ethanol') == \
        '1-{2-[3-(4-methylphenyl)propyl]phenyl}ethanol'
    assert chem.nest_enclosures('(2S)-2-aminopropanoic acid') == '(2S)-2-aminopropanoic acid'
    assert chem.nest_enclosures('2,6(1H)-dicarboxylic acid') == '2,6(1H)-dicarboxylic acid'


def test_structure_to_name_caption(win):
    c = win.canvas
    c.set_doc(_doc('Oc1ccc(cc1)[N+](=O)[O-]'))
    win.info.refresh()
    assert win.info.name.text() == '4-nitrophenol'
    win.structure_to_iupac()
    assert [t.text for t in c.doc.texts.values()] == ['4-nitrophenol']
    win.info.paste('name')
    assert len(c.doc.texts) == 2
