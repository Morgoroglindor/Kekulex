"""Isotope labels (13C, 15N, D), Add missing charges, partial charges, running curved arrows, and annotating a measured
mass spectrum."""
import base64
import struct

import pytest
from PyQt6 import QtCore
from PyQt6.QtTest import QTest

import cdxml
import chem
import model
import ms
import msfile
import reactions
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


def _doc(smiles):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text(smiles), (4, 4))
    return d


def _type_label(c, a, text):
    c.edit_label(a)
    c.editor[0].setText(text)
    c.commit_editor()


def _shown(doc, a):
    lab = render.Renderer(doc).labels.get(a)
    return ''.join(p[0] for p in sorted(lab['pieces'], key=lambda p: p[1])) if lab else ''


# ------------------------------------------------------------------------------------------------------- isotopes
@pytest.mark.parametrize('typed, iso, el', [('13C', 13, 'C'), ('[15N]H2', 15, 'N'), ('¹⁸O', 18, 'O'), ('D', 2, 'H'),
                                            ('13ch3', 13, 'C'), ('37Cl', 37, 'Cl')])
def test_label_isotope(typed, iso, el):
    got, rest = model.label_isotope(typed)
    assert got == iso and model.parse_label(rest)[1] == el


def test_no_such_isotope():
    assert model.label_isotope('40C') == (0, '40C')
    assert model.label_isotope('CH3') == (0, 'CH3')


def test_typed_13c_counts_in_formula_mass_smiles_and_shows(win):
    c = win.canvas
    c.set_doc(_doc('CCO'))
    end = next(a for a, at in c.doc.atoms.items() if at.el == 'C' and len(c.doc.neighbors(a)) == 1)
    plain = chem.analyse(c.doc)
    _type_label(c, end, '13CH3')
    at = c.doc.atoms[end]
    assert (at.el, at.iso) == ('C', 13)
    info = chem.analyse(c.doc)
    assert info.formula == 'C¹³CH6O' and '<sup>13</sup>C' in info.formula_html
    assert info.exact == pytest.approx(plain.exact + 1.003355, abs=1e-5)
    assert '[13CH3]' in info.smiles
    assert _shown(c.doc, end) in ('13CH3', 'H313C')
    assert chem.isotope_pattern(info.counts)[0][0] == pytest.approx(info.exact, abs=1e-6)
    c.edit_label(end)
    assert c.editor[0].text().startswith('13C')
    c.commit_editor(cancel=True)


def test_deuterium_label_and_mass(win):
    c = win.canvas
    d = Doc()
    h = d.add_atom(5, 4, 'H')
    d.add_bond(d.add_atom(4, 4), h)
    c.set_doc(d)
    _type_label(c, h, 'D')
    assert c.doc.atoms[h].iso == 2 and _shown(c.doc, h) == 'D'
    info = chem.analyse(c.doc)
    assert info.formula == 'CH3D'
    assert info.exact == pytest.approx(12.0 + 3 * 1.007825 + 2.014102, abs=1e-5)


def test_isotope_survives_molfile_and_cdxml(tmp_path):
    d = _doc('[13CH3]C(=O)[15NH2]')
    block = d.to_molblock()[0]
    assert 'M  ISO' in block
    back = Doc(); back.insert_rdkit(chem.mol_from_text(block), (0, 0))
    assert sorted(at.iso for at in back.atoms.values() if at.iso) == [13, 15]
    out = tmp_path / 'iso.cdxml'
    cdxml.write(d, str(out))
    assert 'Isotope="13"' in out.read_text(encoding='utf-8')
    again, _notes = cdxml.read(str(out))
    assert sorted(at.iso for at in again.atoms.values() if at.iso) == [13, 15]
    assert chem.analyse(again).formula == chem.analyse(d).formula


def test_isotope_menu_choices():
    import kekulex as se
    assert 13 in se.isotope_choices('C') and 14 in se.isotope_choices('C')
    assert 2 in se.isotope_choices('H') and 3 in se.isotope_choices('H')
    assert 15 in se.isotope_choices('N') and 18 in se.isotope_choices('O')


def test_another_element_clears_the_isotope(win):
    c = win.canvas
    c.set_doc(_doc('CC'))
    a = next(iter(c.doc.atoms))
    c.doc.atoms[a].iso = 13
    c.hover = ('atom', a)
    QTest.keyClick(c, 'o')
    assert c.doc.atoms[a].el == 'O' and c.doc.atoms[a].iso == 0


# ------------------------------------------------------------------------------------------------- missing charges
def test_add_missing_charges(win):
    c = win.canvas
    d = _doc('C[N+](C)(C)C.CC')
    n = next(a for a, at in d.atoms.items() if at.el == 'N')
    d.atoms[n].charge = 0                                   # drawn without its charge
    end = next(a for a, at in d.atoms.items() if at.el == 'C' and len(d.neighbors(a)) == 1
               and all(d.atoms[x].el == 'C' for x, _b in d.neighbors(a)) and d.bond_sum(a) == 1
               and not any(d.atoms[x].el == 'N' for x, _b in d.neighbors(a)))
    c.set_doc(d)
    _type_label(c, end, 'O')                                # ethanol-like O ...
    o = end
    c.doc.atoms[o].hcount = 0                               # ... typed with no H: an alkoxide
    assert c.doc.needed_charge(n) == 1 and c.doc.needed_charge(o) == -1
    assert any('needs a charge' in p for p in c.check_structure())
    done = c.fix_charges()
    assert done == {n: 1, o: -1}
    assert c.doc.atoms[n].charge == 1 and c.doc.atoms[o].charge == -1
    assert not c.doc.valence_problem(n) and not c.doc.valence_problem(o)
    assert chem.analyse(c.doc, {o} | {x for x, _b in c.doc.neighbors(o)}).charge == -1
    c.undo()
    assert c.doc.atoms[n].charge == 0


def test_typed_nh3_becomes_ammonium():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.atoms[a].el, d.atoms[a].hcount = 'N', 3
    assert d.fix_charges() == {a: 1}
    assert d.atoms[a].hcount is None and d.implicit_h(a) == 3


def test_carbon_is_left_alone():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.atoms[a].hcount = 2                                   # CH2 on a chain end: a cation or an anion?
    assert d.needed_charge(a) is None and d.fix_charges() == {}


# ------------------------------------------------------------------------------------------------- partial charges
def test_partial_charges():
    d = _doc('CC(=O)OC')
    q = chem.partial_charges(d)
    assert len(q) == len(d.atoms) and abs(sum(q.values())) < 1e-6
    carbonyl_o = next(a for a, at in d.atoms.items() if at.el == 'O' and len(d.neighbors(a)) == 1)
    carbonyl_c = d.neighbors(carbonyl_o)[0][0]
    assert q[carbonyl_o] < -0.2 and q[carbonyl_c] > 0.2


def test_partial_charges_with_a_nickname():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.set_nickname(a, 'OMe')
    q = chem.partial_charges(d)
    assert set(q) == set(d.atoms) and abs(sum(q.values())) < 1e-6


def test_partial_charge_view(win):
    c = win.canvas
    c.set_doc(_doc('CC(=O)Cl'))
    win.a_partial.setChecked(True)
    assert c.show_partial
    c.repaint()
    first = c.partial_charges()
    assert c.partial_charges() is first                     # cached while the drawing is the same
    win.a_partial.setChecked(False)


# ------------------------------------------------------------------------------------------------- curved arrows
def _arrow(d, kind, tail, head):
    i = d.add_arrow(0, 0, 1, 1, kind)
    d.attach_arrow_end(d.arrows[i], 'tail', tail)
    d.attach_arrow_end(d.arrows[i], 'head', head)
    d.update_attached_arrows()
    return i


def _sn2():
    d = _doc('[OH-].CBr')
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    cc = next(a for a, at in d.atoms.items() if at.el == 'C')
    br = next(a for a, at in d.atoms.items() if at.el == 'Br')
    _arrow(d, 'curved', ('atom', o), ('atom', cc))
    _arrow(d, 'curved', ('bond', d.bond_between(cc, br)), ('atom', br))
    return d


def test_sn2_step():
    d = _sn2()
    new, changed, problems = reactions.run_arrows(d, list(d.arrows))
    assert not problems
    assert chem.analyse(new).smiles == 'CO.[Br-]'
    assert len(changed) == 3


def test_homolysis_and_radical_recombination():
    d = _doc('BrBr')
    a, b = list(d.atoms)
    bb = d.bond_between(a, b)
    for h in (a, b):
        _arrow(d, 'curved_half', ('bond', bb), ('atom', h))
    new, _c, problems = reactions.run_arrows(d, list(d.arrows))
    assert not problems and all(at.radical == 1 for at in new.atoms.values()) and not new.bonds
    d2 = new
    d2.arrows.clear()
    _arrow(d2, 'curved_half', ('atom', a), ('atom', b))
    _arrow(d2, 'curved_half', ('atom', b), ('atom', a))
    back, _c, problems = reactions.run_arrows(d2, list(d2.arrows))
    assert not problems and chem.analyse(back).smiles == 'BrBr'


def test_carbonyl_addition_needs_both_arrows():
    d = _doc('CC(C)=O.[C-]#N')
    nc = next(a for a, at in d.atoms.items() if at.el == 'C' and at.charge == -1)
    co = next(a for a, at in d.atoms.items() if at.el == 'C' and len(d.neighbors(a)) == 3)
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    _arrow(d, 'curved', ('atom', nc), ('atom', co))
    _new, _c, problems = reactions.run_arrows(d, list(d.arrows))
    assert problems and 'bonds' in problems[0]
    _arrow(d, 'curved', ('bond', d.bond_between(co, o)), ('atom', o))
    new, _c, problems = reactions.run_arrows(d, list(d.arrows))
    assert not problems and chem.analyse(new).smiles == 'CC(C)([O-])C#N'


def test_wrong_arrows_are_named():
    d = _doc('C.O')
    cc, o = list(d.atoms)
    _arrow(d, 'curved', ('atom', cc), ('atom', o))           # methane has no lone pair
    assert 'no lone pair' in reactions.run_arrows(d, list(d.arrows))[2][0]
    d = _doc('BrBr')
    a, b = list(d.atoms)
    _arrow(d, 'curved_half', ('bond', d.bond_between(a, b)), ('atom', a))    # one fish-hook alone
    assert any('partner' in p for p in reactions.run_arrows(d, list(d.arrows))[2])
    d = _doc('CC')
    i = d.add_arrow(30, 30, 31, 31, 'curved')               # far from everything
    assert 'does not start' in reactions.run_arrows(d, [i])[2][0]


def test_unattached_arrow_ends_are_found():
    d = _doc('[OH-].CBr')
    o = next(a for a, at in d.atoms.items() if at.el == 'O')
    cc = next(a for a, at in d.atoms.items() if at.el == 'C')
    po, pc = d.atoms[o], d.atoms[cc]
    i = d.add_arrow(po.x + 0.3, po.y - 0.3, pc.x - 0.2, pc.y - 0.1, 'curved')   # drawn from O's lone pair, loose
    assert reactions.arrow_ends(d, d.arrows[i]) == (('atom', o), ('atom', cc))


def test_run_arrows_draws_the_next_step(win):
    c = win.canvas
    c.set_doc(_sn2())
    before = len(c.doc.atoms)
    items, problems = c.run_arrows()
    assert not problems
    new_atoms = {i for k, i in items if k == 'atom'}
    assert len(c.doc.atoms) == before * 2 and len(new_atoms) == before
    assert chem.analyse(c.doc, new_atoms).smiles == 'CO.[Br-]'
    arrows = [c.doc.arrows[i] for k, i in items if k == 'arrow']
    assert arrows and arrows[0].kind == 'forward'
    old_right = max(c.doc.atoms[a].x for a in c.doc.atoms if a not in new_atoms)
    assert min(c.doc.atoms[a].x for a in new_atoms) > old_right
    c.undo()
    assert len(c.doc.atoms) == before


def test_check_structure_names_bad_arrows(win):
    c = win.canvas
    d = _doc('C.O')
    cc, o = list(d.atoms)
    _arrow(d, 'curved', ('atom', cc), ('atom', o))
    c.set_doc(d)
    assert any(p.startswith('Curved arrows:') for p in c.check_structure())


# ------------------------------------------------------------------------------------------------- spectra
def _mzml(tmp_path, scans):
    def arr(vals):
        return base64.b64encode(struct.pack(f'<{len(vals)}d', *vals)).decode()
    body = []
    for k, (rt, mz, inten) in enumerate(scans):
        body.append(f'''<spectrum index="{k}" id="scan={k + 1}" defaultArrayLength="{len(mz)}">
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
<scanList count="1"><scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt}"
 unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/></scan></scanList>
<binaryDataArrayList count="2">
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/><binary>{arr(mz)}</binary></binaryDataArray>
<binaryDataArray encodedLength="0"><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/><binary>{arr(inten)}</binary>
</binaryDataArray></binaryDataArrayList></spectrum>''')
    text = ('<?xml version="1.0" encoding="utf-8"?><mzML xmlns="http://psi.hupo.org/ms/mzml"><run id="r">'
            f'<spectrumList count="{len(scans)}">' + ''.join(body) + '</spectrumList></run></mzML>')
    p = tmp_path / 'run.mzML'
    p.write_text(text, encoding='utf-8')
    return str(p)


def test_spectrum_files(tmp_path):
    p = tmp_path / 'a.csv'
    p.write_text('m/z;Intensity\n43,018;100\n58,041;30\n', encoding='utf-8')
    s = msfile.read(str(p))
    assert list(s[0].mz) == [43.018, 58.041]
    jdx = msfile.parse_text('##TITLE=x\n##PEAK TABLE=(XY..XY)\n15,90 43,999\n58,300\n##END=')
    assert list(jdx[0].mz) == [15, 43, 58] and msfile.unit_resolution(msfile.centroid(jdx[0]))
    msp = msfile.parse_text('Name: acetone\nNum Peaks: 3\n15 90; 43 999; 58 300\n')
    assert msp[0].title == 'acetone' and len(msp[0].mz) == 3
    mgf = msfile.parse_text('BEGIN IONS\nTITLE=t\nPEPMASS=73.06\n43.01 10\n55.05 5\nEND IONS\n')
    assert mgf[0].precursor == pytest.approx(73.06) and mgf[0].level == 2
    path = _mzml(tmp_path, [(1.5, [73.0648, 95.0467], [100.0, 20.0]), (2.0, [73.0648], [1000.0])])
    run = msfile.read(path)
    assert len(run) == 2 and run[1].rt == pytest.approx(2.0) and run[1].tic == 1000.0


def test_profile_spectrum_is_centroided():
    import numpy as np
    x = np.arange(70.0, 76.0, 0.002)
    y = 1000 * np.exp(-((x - 73.0648) / 0.004) ** 2) + 300 * np.exp(-((x - 74.068) / 0.004) ** 2)
    peaks = msfile.centroid(msfile.Spectrum('p', x, y), 0.05)
    assert len(peaks) == 2 and peaks[0][0] == pytest.approx(73.0648, abs=5e-4)


def test_annotate_butanone_ei():
    d = _doc('CCC(C)=O')
    cands = ms.candidate_ions(d, [('1', set(d.atoms))], positive=True)
    rows = ms.annotate([(43.0178, 100), (57.0335, 40), (72.0570, 30), (99.9, 5)], cands, 10)
    labels = [r[2][2] if r[2] else None for r in rows]
    assert labels[:3] == ['C2H3O', 'C3H5O', 'C4H8O'] and labels[3] is None
    assert rows[2][2][1] == '[M]•+'


def test_spectrum_dialog(win, tmp_path):
    c = win.canvas
    c.set_doc(_doc('CCC(C)=O'))
    path = _mzml(tmp_path, [(1.0, [43.0178, 57.0335, 73.0648], [100.0, 40.0, 30.0]),
                            (2.0, [73.0648, 95.0467], [1000.0, 300.0])])
    dlg = win.annotate_spectrum(path)
    assert dlg.which.currentIndex() == 1                    # the scan with the most ions
    matched = {r[2][1] for r in dlg.rows if r[2]}
    assert {'[M+H]+', '[M+Na]+'} <= matched
    assert 'C4H9O' in dlg.table_text()
    dlg.which.setCurrentIndex(0)
    assert sum(1 for r in dlg.rows if r[2]) == 3
    n = len(c.doc.shapes)
    dlg.insert()
    assert len(c.doc.shapes) == n + 1 and next(iter(c.doc.shapes.values())).kind == 'image'
    dlg.close()


def _menu_actions(menu):
    out = {}
    for a in menu.actions():
        if a.menu():
            out.update({a.text().replace('&', '') + '/' + k: v for k, v in _menu_actions(a.menu()).items()})
        elif a.text():
            out[a.text().replace('&', '')] = a
    return out


def test_right_click_isotope_and_arrow_entries(win, monkeypatch):
    from PyQt6 import QtWidgets
    seen = {}
    monkeypatch.setattr(QtWidgets.QMenu, 'exec', lambda self, *a: seen.update(_menu_actions(self)))
    c = win.canvas
    c.set_doc(_sn2())
    cc = next(a for a, at in c.doc.atoms.items() if at.el == 'C')
    c.sel = {('atom', cc)}
    win.object_menu(QtCore.QPoint(0, 0), ('atom', cc))
    assert 'Atom/Isotope/¹³C' in seen and 'Atom/Isotope/¹⁴C' in seen
    seen['Atom/Isotope/¹³C'].trigger()
    assert c.doc.atoms[cc].iso == 13
    assert 'Molecule/Add missing charges' in seen
    seen.clear()
    arrow = next(i for i in c.doc.arrows)
    c.sel = {('arrow', arrow)}
    win.object_menu(QtCore.QPoint(0, 0), ('arrow', arrow))
    assert 'Run the curved arrows (next step)' in seen
