"""MS round 3: enol McLafferty ions and the six-membered transition state, sodium / ammonium adducts, isotope peaks,
ranking drawn candidates against a spectrum, searching an LC-MS run, the fragmentation tree, MSP / MGF export."""
import math

import numpy as np
import pytest

import chem
import fragments
import ms
import msfile
from model import Doc

CQA = 'O=C(/C=C/c1ccc(O)c(O)c1)OC1CC(O)(C(=O)O)CC(O)C1O'
CQA_PEAKS = [(353.0878, 18), (354.0911, 3.2), (191.0561, 100), (192.0594, 8.1), (179.0350, 45), (180.0375, 4.6),
             (173.0455, 12), (161.0244, 8), (135.0452, 30)]


def doc_of(*smiles):
    d = Doc()
    for k, s in enumerate(smiles):
        d.insert_rdkit(chem.mol_from_text(s), (6 + 14 * k, 6))
    return d


def test_version_is_0_11():
    import version
    assert version.__version__ == '0.11.0'


def test_mclafferty_ion_is_drawn_as_the_enol():
    from rdkit import Chem
    d = doc_of('CCCCC(C)=O')
    _p, frags = fragments.predict(d, None, 'EI', 1)
    f = next(f for f in frags if f.how == 'McLafferty rearrangement')
    piece, _why = fragments.draw_doc(d, f)
    mol, _o, _e = chem.to_mol(piece)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles('C=C(C)O'))
    assert any(t.text == '•+' for t in piece.texts.values())


def test_mclafferty_precursor_drawn_as_six_membered_ring():
    d = doc_of('CCCCC(C)=O')
    _p, frags = fragments.predict(d, None, 'EI', 1)
    f = next(f for f in frags if f.how == 'McLafferty rearrangement')
    states, _why = fragments.mechanism(d, f, set(d.atoms))
    first = states[0][0]
    h = next(a for a, at in first.atoms.items() if at.el == 'H')
    o = next(a for a, at in first.atoms.items() if at.el == 'O')
    ho = math.hypot(first.atoms[h].x - first.atoms[o].x, first.atoms[h].y - first.atoms[o].y)
    assert ho == pytest.approx(1.0, abs=0.05)               # the H sits next to the O, one ring edge away
    assert chem.analyse(first).formula == 'C6H12O'


def test_ammonium_adduct_goes_through_m_plus_h():
    d = doc_of('CC1OC(Oc2c(-c3ccc(O)c(O)c3)oc3cc(O)cc(O)c3c2=O)C(O)C(O)C1O')
    pre, frags = fragments.predict(d, None, 'ESI+NH4', 2)
    assert pre.label == '[M+NH4]+' and pre.mz == pytest.approx(466.1344, abs=1e-3)
    mh = next(f for f in frags if f.lost == 'NH3')
    assert mh.mz == pytest.approx(449.1078, abs=1e-3)
    y0 = next(f for f in frags if abs(f.mz - 303.0499) < 1e-3)
    assert y0.lost.startswith('NH3 + ') and fragments.chain_of(y0)[0] is mh
    states, _why = fragments.mechanism(d, y0, set(d.atoms))
    assert len(states) == 3                                 # [M+NH4]+ → [M+H]+ → Y0+


def test_sodium_adduct_keeps_sodium_and_cuts_across_sugars():
    d = doc_of('CC1OC(Oc2c(-c3ccc(O)c(O)c3)oc3cc(O)cc(O)c3c2=O)C(O)C(O)C1O')
    pre, frags = fragments.predict(d, None, 'ESI+Na', 2)
    assert pre.counts.get('Na') == 1
    t = {round(f.mz, 3): f for f in frags}
    assert any(abs(m - 325.032) < 0.002 for m in t)         # [Y0+Na]+
    assert any(f.lost == 'C2H4O2' for f in frags)           # cross-ring −60
    pc = doc_of('CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(C)=O')
    _p, frags = fragments.predict(pc, None, 'ESI+Na', 2)
    assert next(f for f in frags if f.lost == 'C3H9N').score >= 0.85   # −59


def test_isotope_peaks_are_recognised_and_fitted():
    d = doc_of(CQA)
    cands = ms.candidate_ions(d, [('1', set(d.atoms))], False, False)
    rows = ms.annotate(CQA_PEAKS, cands, 5.0)
    rows, fits = ms.isotopes(rows, 5.0)
    iso = [r for r in rows if r[2] and r[2][0].startswith('isotope of')]
    assert {round(r[0], 2) for r in iso} >= {354.09, 192.06, 180.04}
    assert 'M+1 8.1 %' in fits[191.0561]
    assert ms.ion_charge('[M+2H]2+') == 2 and ms.ion_charge('C7H11O6−') == -1


def test_ranking_prefers_the_right_candidate():
    d = doc_of(CQA, 'O=C(O)/C=C/c1ccc(OC2OC(CO)C(O)C(O)C2O)c(O)c1')
    units = [(f'{k}', set(u)) for k, u in enumerate(d.fragments(), 1)]
    ranked = ms.rank_structures(d, units, CQA_PEAKS, False, 5.0)
    best = ranked[0]
    assert chem.analyse(d, best['ids']).formula == 'C16H18O9'
    assert best['explained'] > 0.85 and best['score'] > ranked[1]['score']


def test_search_run_finds_msms_and_ms1():
    d = doc_of(CQA)
    ms2 = msfile.Spectrum('scan 2', np.array([p[0] for p in CQA_PEAKS]), np.array([float(p[1]) for p in CQA_PEAKS]),
                          rt=4.2, level=2, precursor=353.0878)
    ms1 = msfile.Spectrum('scan 1', np.array([200.0, 353.0878]), np.array([10.0, 5000.0]), rt=4.1, level=1)
    other = msfile.Spectrum('scan 3', np.array([100.0]), np.array([5.0]), rt=6.0, level=2, precursor=401.2)
    hits = ms.search_run(d, set(d.atoms), [ms1, ms2, other], False, 10.0)
    assert hits[0]['kind'] == 'MS2' and hits[0]['index'] == 1 and hits[0]['explained'] > 0.8
    assert any(h['kind'] == 'MS1' and h['index'] == 0 for h in hits)
    assert all(h['index'] != 2 for h in hits)


def test_tree_scheme_and_export(qapp, tmp_path):
    from kekulex import MainWindow
    import fragments_tab
    w = MainWindow()
    w.canvas.set_doc(doc_of(CQA))
    w.show_fragments()
    tab = w.info.frag_tab
    tab.mode.setCurrentIndex(tab.mode.findData('ESI-'))
    tab.level.setCurrentIndex(1)
    w.info.refresh(keep_choice=True)
    links = dict(tab.tree_links())
    k135 = next(k for k, f in enumerate(tab.frags) if abs(f.mz - 135.0452) < 1e-3)
    assert abs(tab.frags[links[k135]].mz - 179.0350) < 1e-3   # 135 hangs under 179
    tab.tree_b.setChecked(True)
    root = tab.tree.topLevelItem(0)
    assert root.childCount() >= 3
    texts = len(w.canvas.doc.texts)
    tab.draw_scheme()
    assert len(w.canvas.doc.texts) - texts >= 2 * len(tab.frags) + 1
    msp = fragments_tab.spectrum_text(tab.pre, tab.frags, 'C16H18O9', 'x', False, 'ESI-')
    assert 'Precursor_type: [M-H]-' in msp and f'Num Peaks: {len({round(f.mz, 4) for f in tab.frags})}' in msp
    mgf = fragments_tab.spectrum_text(tab.pre, tab.frags, 'C16H18O9', 'x', True, 'ESI-')
    assert mgf.startswith('BEGIN IONS') and 'PEPMASS=353.0878' in mgf and 'CHARGE=1-' in mgf
    assert msfile.parse_text(msp)[0].mz.size == len({round(f.mz, 4) for f in tab.frags})   # reads back
    w.dirty = False
    w.close()


def test_spectrum_dialog_ranks_and_offers_run_search(qapp):
    from kekulex import MainWindow
    from spectrum_dialog import SpectrumDialog
    w = MainWindow()
    w.canvas.set_doc(doc_of(CQA))
    dlg = SpectrumDialog(w)
    dlg.mode.setCurrentIndex(1)
    ms2 = msfile.Spectrum('s', np.array([p[0] for p in CQA_PEAKS]), np.array([float(p[1]) for p in CQA_PEAKS]),
                          rt=4.2, level=2, precursor=353.0878)
    ms1 = msfile.Spectrum('m', np.array([353.0878]), np.array([5000.0]), rt=4.1, level=1)
    dlg.set_spectra([ms1, ms2], 'test')
    dlg.which.setCurrentIndex(1)
    assert 'explains' in dlg.rank.text() and dlg.find_b.isEnabled()
    assert any(dlg.table.item(r, 7) and 'M+1' in dlg.table.item(r, 7).text() for r in range(dlg.table.rowCount()))
    dlg.find_in_run()
    t = dlg.run_dialog.findChild(__import__('PyQt6.QtWidgets', fromlist=['QTableWidget']).QTableWidget)
    assert t.rowCount() >= 2
    dlg.run_dialog.close()
    dlg.close()
    w.dirty = False
    w.close()
