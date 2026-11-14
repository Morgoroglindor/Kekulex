"""MS fragments, m/z matching, attached arrows and lone pairs, find, reaction maps / balance, tautomers and
resonance, stereo checks, snap guides, batch export, reaction conditions and schemes, polymer brackets in .cdxml."""
import os

import pytest
from PyQt6 import QtCore
from PyQt6.QtTest import QTest

import batch
import cdxml
import chem
import ms
import reactions
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


def _doc(smiles, at=(3, 3)):
    d = Doc(); d.insert_rdkit(chem.mol_from_text(smiles), at); return d


def _acetone_mplus():
    d = _doc('CC(C)=O')
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    d.atoms[o].charge, d.atoms[o].radical = 1, 1
    co = next(a for a, at in d.atoms.items() if at.el == 'C' and len(d.neighbors(a)) == 3)
    b = next(bb for n, bb in d.neighbors(co) if d.atoms[n].el == 'C')
    return d, b


# ------------------------------------------------------------------------------------------- 1 cleave
def test_alpha_cleavage_of_acetone_radical_cation():
    d, b = _acetone_mplus()
    found = {}
    for spec, _how in ms.options(d, b):
        found[spec[:2]] = [ms.describe(p) for p, _cut in ms.pieces(d, b, spec)]
    texts = ' | '.join(' + '.join(v) for v in found.values())
    assert 'C2H3O+  m/z 43.0178' in texts                   # acylium
    assert 'CH3•  15.0235 (neutral)' in texts


def test_ring_bond_cannot_be_cleaved():
    d = _doc('C1CCCCC1')
    assert ms.options(d, next(iter(d.bonds))) == []


def test_inductive_cleavage_of_a_protonated_alcohol():
    d = _doc('CC[OH2+]')
    co = next(b for b, bd in d.bonds.items() if 'O' in (d.atoms[bd.a].el, d.atoms[bd.b].el))
    descs = [[ms.describe(p) for p, _c in ms.pieces(d, co, spec)] for spec, how in ms.options(d, co)
             if how.startswith('inductive')]
    assert descs and 'C2H5+  m/z 29.0386' in descs[0] and 'H2O  18.0106 (neutral)' in descs[0]


def test_cleave_draws_both_fragments_with_live_captions(win):
    c = win.canvas
    d, b = _acetone_mplus()
    c.set_doc(d)
    k = next(k for k, (spec, _h) in enumerate(ms.options(d, b))
             if 'C2H3O+' in ' '.join(ms.describe(p) for p, _c in ms.pieces(d, b, spec)))
    n_atoms = len(c.doc.atoms)
    new = c.cleave_bond(b, k)
    assert len(c.doc.atoms) == n_atoms * 2                  # the two pieces hold the parent's atoms again
    caps = [c.doc.texts[i].text for kk, i in new if kk == 'text']
    assert any('m/z 43.0178' in t for t in caps) and '+' in caps
    assert any(kk == 'arrow' for kk, _i in new)


# ------------------------------------------------------------------------------------------- 2 m/z
def test_compositions_find_glucose_mh():
    found = ms.compositions(181.0707, 1, 3)
    assert found[0][0] == 'C6H13O6' and found[0][4] == 'even'


def test_drawn_matches_and_limits():
    d = _doc('CC(=O)O')
    units = [('acetic acid', set(d.atoms))]
    hits = ms.drawn_matches(d, units, 61.0284, 5)              # [M+H]+
    assert hits and hits[0][1] == '[M+H]+'
    assert ms.limits_from_text('C0-40 H N4 Cl') == {'C': 40, 'H': 160, 'N': 4, 'Cl': 6}


def test_match_dialog(win):
    from ms_dialog import MatchDialog
    win.canvas.set_doc(_doc('CC(=O)O'))
    dlg = MatchDialog(win)
    dlg.mz.setText('61.0284')
    dlg.search()
    assert dlg.drawn.rowCount() >= 1 and dlg.comps.rowCount() >= 1
    assert dlg.comps.item(0, 0).text() == 'C2H5O2+'


# ------------------------------------------------------------------------------------------- 3 arrows, lone pairs
def test_curved_arrow_ends_stay_on_their_atoms(win):
    c = win.canvas
    d = _doc('C=CC')
    c.set_doc(d)
    a, b = list(d.atoms)[:2]
    j = d.add_arrow(d.atoms[a].x + 0.2, d.atoms[a].y - 0.3, d.atoms[b].x, d.atoms[b].y - 0.3, 'curved', 0.3)
    ar = d.arrows[j]
    d.attach_arrow_end(ar, 'tail', ('atom', a))
    d.attach_arrow_end(ar, 'head', ('bond', next(iter(d.bonds))))
    c.snapshot()
    d.transform([('atom', x) for x in d.atoms], lambda x, y: (x + 2, y + 1))
    c.edited()
    assert (ar.x1, ar.y1) == pytest.approx((c.doc.atoms[a].x + 0.2, c.doc.atoms[a].y - 0.3))
    back = Doc.from_dict(c.doc.to_dict())
    assert back.arrows[j].tail_on[:2] == ['atom', a]


def test_lone_pairs_counted_and_drawn():
    d = _doc('CO')
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    assert d.lone_pair_count(o) == 2
    d.atoms[o].lone_pairs = -1
    dots = render.Renderer(d).lone_pair_dots(o)
    assert len(dots) == 4
    n = _doc('C[NH3+]')
    assert n.lone_pair_count(next(a for a, at in n.atoms.items() if at.el == 'N')) == 0


# ------------------------------------------------------------------------------------------- 4 find
def test_find_substructure_label_and_text(win):
    c = win.canvas
    d = _doc('CC(=O)O')
    d.insert_rdkit(chem.mol_from_text('c1ccccc1'), (10, 3))
    t = d.add_text(3, 6, 'Step 2: hydrolysis')
    c.set_doc(d)
    assert c.find('c1ccccc1') == 6
    assert c.find('C(=O)[OH]') == 3
    assert c.find('hydrolysis') == 1 and ('text', t) in c.sel
    assert c.find('xyz123') == 0


# ------------------------------------------------------------------------------------------- 5 map, balance
def _esterification():
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CC(=O)O'), (0, 0))
    d.insert_rdkit(chem.mol_from_text('CCO'), (4, 0))
    d.add_arrow(6, 0, 9, 0)
    d.insert_rdkit(chem.mol_from_text('CC(=O)OCC'), (12, 0))
    return d


def test_atom_map_and_balance():
    d = _esterification()
    assert reactions.map_atoms(d) >= 5
    maps = [at.map for at in d.atoms.values() if at.map]
    assert all(maps.count(m) == 2 for m in set(maps))       # each number once on each side
    text, diff = reactions.balance(d)
    assert diff == {'H': -2, 'O': -1} and 'H2O' in text
    d.insert_rdkit(chem.mol_from_text('O'), (16, 0))
    assert reactions.balance(d)[0].startswith('Balanced')


def test_map_numbers_are_drawn(win):
    c = win.canvas
    c.set_doc(_esterification())
    assert c.map_reaction() >= 5
    img = render.to_image(c.doc, 150, False)
    assert not img.isNull()
    c.clear_map()
    assert not any(at.map for at in c.doc.atoms.values())


# ------------------------------------------------------------------------------------------- 6 forms
def test_tautomers_and_resonance(win):
    from rdkit import Chem
    d = _doc('CC(=O)CC(C)=O')
    smis = {Chem.MolToSmiles(m) for m in reactions.tautomers(d, set(d.atoms))}
    assert 'CC(=O)C=C(C)O' in smis
    a = _doc('CC(=O)[O-]')
    assert len(reactions.resonance_forms(a, set(a.atoms))) == 1
    c = win.canvas
    c.set_doc(_doc('CC(=O)[O-]'))
    new = c.add_forms('resonance')
    assert any(k == 'arrow' and c.doc.arrows[i].kind == 'resonance' for k, i in new)
    assert len(c.doc.fragments()) == 2


# ------------------------------------------------------------------------------------------- 7 stereo
def test_stereo_problems():
    d = _doc('CC(O)CC')
    assert any('stereocentre' in p for p in reactions.stereo_problems(d, set(d.atoms)))
    e = _doc('C[C@H](O)CC')
    assert not reactions.stereo_problems(e, set(e.atoms))
    f = _doc('CCCC')
    b = next(iter(f.bonds))
    f.bonds[b].stereo = 'wedge'
    assert any('not a stereocentre' in p for p in reactions.stereo_problems(f, set(f.atoms)))


# ------------------------------------------------------------------------------------------- 8 snap
def test_snap_lines_up_a_dragged_structure(win):
    c = win.canvas
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CCCC'), (3, 3))
    d.insert_rdkit(chem.mol_from_text('CCCC'), (10, 5.07))
    c.set_doc(d)
    second = [('atom', a) for f in d.fragments() for a in f if d.atoms[a].x > 7]
    targets = c._snap_targets(second)
    (dx, dy), guides = c._snap(targets, (0.0, -2.0))      # dragged up to almost level: snaps level
    assert guides and abs(dy + 2.07) < 1e-6
    (_dx, dy2), none = c._snap(targets, (0.0, -1.0))      # far off: free
    assert dy2 == -1.0 and not [g for g in none if g[0] == 'h']


# ------------------------------------------------------------------------------------------- 9 batch
def test_batch_export(tmp_path):
    src = tmp_path / 'list.smi'
    src.write_text('CCO ethanol\nc1ccccc1 benzene\nCC(=O)O acetic acid\nnot_a_smiles x\n', encoding='utf-8')
    items = batch.read(str(src))
    assert [n for n, _m in items] == ['ethanol', 'benzene', 'acetic acid']
    files = batch.export_each(items, str(tmp_path / 'out'), '.svg')
    assert sorted(os.path.basename(f) for f in files) == ['acetic acid.svg', 'benzene.svg', 'ethanol.svg']
    grid = batch.grid_doc(items, columns=2)
    assert len(grid.fragments()) == 3 and {t.text for t in grid.texts.values()} == {'ethanol', 'benzene', 'acetic acid'}
    out = tmp_path / 'grid.png'
    batch.export_grid(items, str(out), 2, dpi=96)
    assert out.stat().st_size > 1000


# ------------------------------------------------------------------------------------------- 10 schemes
def test_conditions_sit_over_and_under_the_arrow_and_move_with_it(win):
    c = win.canvas
    d = _esterification()
    c.set_doc(d)
    arrow = next(iter(d.arrows))
    c.set_conditions(arrow, 'H2SO4 (cat.)', 'reflux, 3 h')
    have = c.conditions_of(arrow)
    ar = c.doc.arrows[arrow]
    above, below = c.doc.texts[have['above']], c.doc.texts[have['below']]
    assert above.y < ar.y1 < below.y and above.x == pytest.approx((ar.x1 + ar.x2) / 2)
    assert c.doc.group_of(('text', have['above'])) == c.doc.group_of(('arrow', arrow)) is not None


def test_arrange_scheme_lines_structures_up(win):
    c = win.canvas
    d = _esterification()
    frags = d.fragments()
    d.transform([('atom', a) for a in frags[1]], lambda x, y: (x, y + 3))     # one out of line
    c.set_doc(d)
    c.arrange_scheme()
    rend = c.renderer()
    centres = []
    for f in c.doc.fragments():
        ys = [c.doc.atoms[a].y for a in f]
        centres.append((min(ys) + max(ys)) / 2)
    assert max(centres) - min(centres) < 0.6
    ar = next(iter(c.doc.arrows.values()))
    assert ar.y1 == ar.y2


# ------------------------------------------------------------------------------------------- 11 polymer .cdxml
def test_polymer_brackets_round_trip_through_cdxml(tmp_path):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text('CCCC'), (2, 2))
    xs = sorted(d.atoms, key=lambda a: d.atoms[a].x)
    x1 = (d.atoms[xs[0]].x + d.atoms[xs[1]].x) / 2
    x2 = (d.atoms[xs[2]].x + d.atoms[xs[3]].x) / 2
    ys = [d.atoms[a].y for a in d.atoms]
    s = d.add_shape(x1, min(ys) - 0.8, x2, max(ys) + 0.8, 'brackets')
    d.shapes[s].params['repeat'] = 'n'
    path = str(tmp_path / 'p.cdxml')
    cdxml.write(d, path)
    text = open(path, encoding='utf-8').read()
    assert '<bracketedgroup' in text and 'SRULabel="n"' in text
    back, _notes = cdxml.read(path)
    assert [sh.params.get('repeat') for sh in back.shapes.values()] == ['n']
    assert chem.analyse(back).formula == 'C2H6(C2H4)n'
