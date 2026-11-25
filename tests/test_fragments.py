"""The Fragments tab: rule-based MS fragments (fragments.py), drawing them, and the tab in the formula panel."""
import pytest

import chem
import fragments
import ms
from model import Doc


def doc_of(smiles):
    d = Doc()
    d.insert_rdkit(chem.mol_from_text(smiles), (6, 6))
    return d


def mzs(smiles, mode, steps=2):
    pre, frags = fragments.predict(doc_of(smiles), None, mode, steps)
    return pre, frags, {round(f.mz, 3): f for f in frags}


def has(table, mz, tol=0.002):
    return next((f for k, f in table.items() if abs(k - mz) <= tol), None)


# ------------------------------------------------------------------------------------------- EI
def test_ei_2_hexanone_alpha_cleavage_and_mclafferty():
    pre, _f, t = mzs('CCCCC(C)=O', 'EI')
    assert pre.ion == 'C6H12O•+' and pre.mz == pytest.approx(100.0883, abs=1e-4)
    assert has(t, 43.0178).how.startswith('α-cleavage at C=O')           # CH3C≡O+
    assert has(t, 85.0648).ion == 'C5H9O+'                               # loss of CH3•
    mc = has(t, 58.0413)
    assert mc.how == 'McLafferty rearrangement' and mc.odd and mc.lost == 'C3H6'
    assert mc.likelihood == 'likely'


def test_ei_propylbenzene_tropylium_and_mclafferty_92():
    _p, _f, t = mzs('CCCc1ccccc1', 'EI')
    trop = has(t, 91.0542)
    assert 'tropylium' in trop.how and trop.lost == 'C2H5•'
    assert has(t, 92.0621).odd
    assert has(t, 65.0386) is not None                                   # C7H7+ − C2H2


def test_ei_amine_iminium_then_alkene():
    _p, _f, t = mzs('CCN(CC)CC', 'EI')
    assert has(t, 86.0964).how.startswith('α-cleavage at N')
    assert 'alkene' in has(t, 58.0651).how                               # 86 − C2H4


def test_ei_retro_diels_alder_cyclohexene():
    _p, _f, t = mzs('C1CC=CCC1', 'EI')
    f = has(t, 54.0464)
    assert f.how.startswith('retro-Diels–Alder') and f.lost == 'C2H4' and f.odd


def test_ei_halide_and_nitro():
    _p, _f, t = mzs('CCCCBr', 'EI')
    assert has(t, 57.0699).how.startswith('inductive cleavage at C–Br')
    _p, _f, t = mzs('O=[N+]([O-])c1ccccc1', 'EI')
    assert has(t, 77.0386).lost == 'NO2•'


# ------------------------------------------------------------------------------------------- ESI
def test_esi_flavone_rda_both_modes():
    apigenin = 'O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12'
    pre, _f, t = mzs(apigenin, 'ESI+')
    assert pre.mz == pytest.approx(271.0601, abs=1e-4)
    assert 'retro-Diels–Alder' in has(t, 153.0182).how                   # 1,3A+
    assert has(t, 119.0491) is not None                                  # 1,3B+
    _p, _f, t = mzs(apigenin, 'ESI-')
    assert 'retro-Diels–Alder' in has(t, 151.0037).how                   # 1,3A−
    assert has(t, 117.0346) is not None and has(t, 225.0557) is not None  # 1,3B−, −CO2


def test_esi_minus_chlorogenic_acid():
    _p, _f, t = mzs('O=C(/C=C/c1ccc(O)c(O)c1)OC1CC(O)(C(=O)O)CC(O)C1O', 'ESI-')
    assert has(t, 191.0561).how.startswith('ester cleavage (alkoxide')    # quinate
    assert has(t, 179.0350).how.startswith('ester cleavage (carboxylate')  # caffeate
    assert has(t, 135.0452) is not None                                  # caffeate − CO2


def test_esi_plus_glycoside_aglycone():
    _p, _f, t = mzs('CC1OC(Oc2c(-c3ccc(O)c(O)c3)oc3cc(O)cc(O)c3c2=O)C(O)C(O)C1O', 'ESI+')   # quercitrin
    y0 = has(t, 303.0499)
    assert 'glycosidic' in y0.how and y0.lost == 'C6H10O4' and y0.score >= 0.85


def test_esi_plus_peptide_b_and_y_ions():
    _p, _f, t = mzs('CC(N)C(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O', 'ESI+')      # Ala-Gly-Phe
    assert 'b-type' in has(t, 129.0659).how                              # b2
    assert 'y-type' in has(t, 166.0863).how                              # y1
    assert 'y-type' in has(t, 223.1077).how                              # y2
    assert has(t, 101.0709) is not None                                  # a2 = b2 − CO


def test_esi_neutral_losses():
    _p, _f, t = mzs('CCCCO', 'ESI+')
    assert has(t, 57.0699).lost == 'H2O'
    _p, _f, t = mzs('CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'ESI-')                # ibuprofen [M−H]− − CO2
    assert has(t, 161.1336).lost == 'CO2'
    _p, _f, t = mzs('CC(C)(C)OC(=O)NCC(=O)O', 'ESI+')                    # Boc: − C4H8
    assert has(t, 120.0291).lost == 'C4H8'


def test_steps_add_consecutive_losses():
    one = mzs('CC(N)C(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O', 'ESI+', 1)[1]
    two = mzs('CC(N)C(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O', 'ESI+', 2)[1]
    assert all(f.step == 1 for f in one)
    assert len(two) > len(one) and any(f.step == 2 and ' + ' in f.lost for f in two)


def test_drawn_ion_used_as_it_is():
    d = doc_of('C[N+](C)(C)CCO')                                         # choline
    pre, frags = fragments.predict(d, None, 'EI', 2)
    assert pre.ion == 'C5H14NO+' and pre.drawn_ion
    assert frags and all(f.z == 1 and not f.odd for f in frags)


SMILES = ['CCCCC(C)=O', 'CCCc1ccccc1', 'CCN(CC)CC', 'C1CC=CCC1', 'O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12',
          'O=C(/C=C/c1ccc(O)c(O)c1)OC1CC(O)(C(=O)O)CC(O)C1O', 'CC(N)C(=O)NCC(=O)NC(Cc1ccccc1)C(=O)O',
          'CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'OP(=O)(O)OCC(N)C(=O)O', 'COc1cc(C=CC(=O)O)ccc1O',
          'CC1OC(Oc2c(-c3ccc(O)c(O)c3)oc3cc(O)cc(O)c3c2=O)C(O)C(O)C1O', 'c1ccc2c(c1)CCCC2', 'C', '[CH3]']


@pytest.mark.parametrize('smiles', SMILES)
def test_drawn_fragments_have_the_listed_formula(smiles):
    """Every fragment that can be drawn has exactly the listed formula, charge and m/z, and no valence problem."""
    d = doc_of(smiles)
    for mode in ('EI', 'ESI+', 'ESI-'):
        pre, frags = fragments.predict(d, None, mode, 3)
        if pre is None:
            continue
        for f in frags:
            assert f.counts and f.mz > 0 and f.lost
            piece, why = fragments.draw_doc(d, f)
            if piece is None:
                assert why
                continue
            info = chem.analyse(piece)
            assert info.charge == f.z, (smiles, mode, f.how)
            assert chem.formula_text(info.counts) == f.formula, (smiles, mode, f.how)
            assert not any(piece.valence_problem(a) for a in piece.atoms), (smiles, mode, f.how)


def test_nickname_split_cannot_be_drawn_but_whole_one_can():
    d = Doc()
    c1 = d.add_atom(0, 0); c2 = d.add_atom(1, 0); o = d.add_atom(2, 0, 'O')
    me = d.add_atom(1.5, 0.9)
    d.add_bond(c1, c2); d.add_bond(c2, o, 2); d.add_bond(c2, me)
    ph = d.add_atom(-1, 0); d.atoms[ph].el, d.atoms[ph].label = 'C', 'Ph'
    d.add_bond(ph, c1)
    pre, frags = fragments.predict(d, None, 'EI', 1)
    assert pre.ion == 'C9H10O•+'                                        # phenylacetone
    trop = next(f for f in frags if f.ion == 'C7H7+')
    piece, _why = fragments.draw_doc(d, trop)
    assert piece is not None and ph in piece.atoms                       # the whole Ph group stays
    acyl = next(f for f in frags if f.ion == 'C2H3O+')
    assert ph in acyl.lost_atoms


def test_spectrum_annotation_uses_rule_fragments():
    d = doc_of('CCCCC(C)=O')
    units = [('Structure', set(d.atoms))]
    cands = ms.candidate_ions(d, units, True, cleave=False)
    hit = ms.annotate([(58.0413, 100.0)], cands, 5.0)[0]
    assert hit[2] is not None and ', fragment: McLafferty' in hit[2][0]


# ------------------------------------------------------------------------------------------- the tab
def test_fragments_tab_lists_highlights_and_draws(qapp):
    from kekulex import MainWindow
    w = MainWindow()
    d = doc_of('CCCCC(C)=O')
    w.canvas.set_doc(d)
    w.show_fragments()
    tab = w.info.frag_tab
    assert w.info.tabs.currentWidget() is tab
    k = tab.mode.findData('EI')
    tab.mode.setCurrentIndex(k)
    w.info.refresh(keep_choice=True)
    assert tab.table.rowCount() >= 3 and 'C6H12O•+' in tab.head.text()
    row = next(r for r in range(tab.table.rowCount()) if tab.table.item(r, 1).text() == 'C3H6O•+')
    tab.table.selectRow(row)
    keep, lost, cut = w.canvas.fragment_view
    assert keep and lost and len(cut) == 1
    n_atoms = len(w.canvas.doc.atoms)
    tab.draw()
    assert len(w.canvas.doc.atoms) == n_atoms + 4                       # C3H6O drawn: 3 C + 1 O
    assert any('−C3H6' in t.text for t in w.canvas.doc.texts.values())
    assert any('m/z 58.04' in t.text for t in w.canvas.doc.texts.values())
    tab.level.setCurrentIndex(0)                                         # likely only
    assert all(f.likelihood == 'likely' for f in tab.frags)
    w.info.tabs.setCurrentIndex(0)
    assert w.canvas.fragment_view is None
    w.dirty = False
    w.close()
