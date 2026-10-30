"""Biopolymer sequences, BioDraw objects and the bio templates."""
import json

import pytest
from PyQt6 import QtWidgets
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem import rdMolDescriptors as D

import biodraw
import biopolymer
import model
import render
import templates
from model import Doc, Sequence


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


def _seq(text, polymer=None, **kw):
    pol, res, extras = biopolymer.parse(text, polymer)
    return Sequence(0, 0, pol, res, **{**extras, **kw})


# ------------------------------------------------------------------------------------------------ chemistry
def test_peptide_builder_matches_rdkit():
    for r in 'ACDEFGHIKLMNPQRSTVWY':
        mine, err = biopolymer.to_mol(Sequence(0, 0, 'PEPTIDE', [r]))
        assert Chem.MolToSmiles(mine) == Chem.MolToSmiles(Chem.MolFromSequence(r)), r
    seq = 'ACDEFGHIKLMNPQRSTVWY'
    mine, _ = biopolymer.to_mol(Sequence(0, 0, 'PEPTIDE', list(seq)))
    assert Chem.MolToSmiles(mine) == Chem.MolToSmiles(Chem.MolFromSequence(seq))


def test_rare_residues_and_links():
    def formula(**kw):
        mol, err = biopolymer.to_mol(Sequence(0, 0, **kw))
        assert mol is not None, err
        return D.CalcMolFormula(mol)
    assert formula(polymer='PEPTIDE', residues=['U']) == 'C3H7NO2Se'               # selenocysteine
    assert formula(polymer='PEPTIDE', residues=['O']) == 'C12H21N3O3'              # pyrrolysine
    assert formula(polymer='PEPTIDE', residues=list('GGGGG'), cyclic=True) == 'C10H15N5O5'
    linear = formula(polymer='PEPTIDE', residues=list('CAAAAC'))
    bridged = formula(polymer='PEPTIDE', residues=list('CAAAAC'), bridges=[[0, 5]])
    assert linear == 'C18H32N6O7S2' and bridged == 'C18H30N6O7S2'                 # S-S: two H fewer
    lactam = formula(polymer='PEPTIDE', residues=list('KAAD'), bridges=[[0, 3]])
    assert lactam == 'C16H27N5O6'                                                  # K-D amide: water lost
    assert formula(polymer='PEPTIDE', residues=list('AG'), n_term='Ac', c_term='NH2') == 'C7H13N3O3'
    mol, err = biopolymer.to_mol(Sequence(0, 0, 'PEPTIDE', list('AAG'), bridges=[[0, 2]]))
    assert mol is None and 'bridge' in err.lower()


def test_pyrrolysine_stereo():
    mol, _ = biopolymer.to_mol(Sequence(0, 0, 'PEPTIDE', ['O']))
    rdCIPLabeler.AssignCIPLabels(mol)
    assert sorted(a.GetProp('_CIPCode') for a in mol.GetAtoms() if a.HasProp('_CIPCode')) == ['R', 'R', 'S']


def test_nucleic_acids():
    dna = _seq('ACGT')
    assert dna.polymer == 'DNA'
    mol, _ = biopolymer.to_mol(dna)
    assert D.CalcMolFormula(mol) == 'C39H50N15O22P3'
    mol, _ = biopolymer.to_mol(_seq('ACGT', p5=True))
    assert D.CalcMolFormula(mol).endswith('P4')
    assert _seq('AUCG', 'DNA').residues == list('ATCG')                           # as ChemDraw: U -> T for DNA
    assert _seq('ATCG', 'RNA').residues == list('AUCG')


def test_parse_formats():
    q = _seq('H-Ala-Gly-Pro-NH2')
    assert (q.polymer, q.residues, q.n_term, q.c_term) == ('PEPTIDE', list('AGP'), 'H', 'NH2')
    assert _seq('Ac-Lys-Ala-OH').n_term == 'Ac'
    assert _seq('>sp|P01308| insulin A\nGIVEQCCTSI CSLYQLENYCN\n').residues == list('GIVEQCCTSICSLYQLENYCN')
    q = _seq('PEPTIDE1{C.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-4:R3$$$')
    assert q.bridges == [[0, 3]]
    q = _seq('PEPTIDE1{G.G.G}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$')
    assert q.cyclic
    q = _seq('RNA1{[dR](A)P.[dR](T)}$$$$')
    assert q.polymer == 'DNA' and q.residues == ['A', 'T']
    with pytest.raises(ValueError):
        biopolymer.parse('Ala-Foo-Gly')
    with pytest.raises(ValueError):
        biopolymer.parse('ACDQ1@', 'DNA')


def test_write_formats_round_trip():
    q = Sequence(0, 0, 'PEPTIDE', list('CAAC'), bridges=[[0, 3]])
    again = _seq(biopolymer.helm(q))
    assert again.residues == q.residues and again.bridges == q.bridges
    q = Sequence(0, 0, 'PEPTIDE', list('GGGG'), cyclic=True)
    assert _seq(biopolymer.helm(q)).cyclic
    dna = Sequence(0, 0, 'DNA', list('ACGT'))
    assert biopolymer.helm(dna) == 'RNA1{[dR](A)P.[dR](C)P.[dR](G)P.[dR](T)}$$$$'
    assert _seq(biopolymer.helm(dna)).residues == list('ACGT')
    assert biopolymer.three_letter(Sequence(0, 0, 'PEPTIDE', list('AG'), c_term='NH2')) == 'H-Ala-Gly-NH2'
    assert biopolymer.fasta(Sequence(0, 0, 'PEPTIDE', list('A' * 70))).count('\n') == 1


def test_analysis_charge_states():
    info = biopolymer.analyse(Sequence(0, 0, 'PEPTIDE', list('GRGDSPK')))
    assert info.formula == 'C28H49N11O11' and abs(info.exact - 715.3613) < 1e-3
    labels = [i.label for i in info.ions]
    assert labels.count('[M+2H]2+') == 1 and '[M+3H]3+' in labels
    m3 = next(i for i in info.ions if i.label == '[M+3H]3+')
    assert abs(m3.mz - (715.36126 + 3 * 1.007276) / 3) < 1e-3
    info = biopolymer.analyse(Sequence(0, 0, 'PEPTIDE', list('AXG')))
    assert info.errors


# ------------------------------------------------------------------------------------------------ drawing
def test_sequence_drawing_and_file_round_trip(tmp_path, qapp):
    d = Doc()
    i = d.add_sequence(1, 1, 'PEPTIDE', list('GIVEQCCTSICSLYQLENYCN'), bridges=[[5, 10]], show_ids=True)
    d.add_sequence(1, 12, 'DNA', list('ATGC'), p5=True)
    r = render.Renderer(d)
    lay = biodraw.seq_layout(r, i)
    assert len(lay['cells']) == 21
    rows = {round(c.top(), 3) for c in lay['cells']}
    assert len(rows) == 3                                                          # 10 per line
    k = biodraw.residue_at(r, i, lay['cells'][4].center().x(), lay['cells'][4].center().y())
    assert k == 4
    assert not r.bounds().isNull()
    img = render.to_image(d, dpi=72)
    assert img.width() > 50
    again = Doc.from_dict(json.loads(json.dumps(d.to_dict())))
    assert again.seqs[i].bridges == [[5, 10]] and again.seqs[i].residues == d.seqs[i].residues


def test_bio_objects_draw_and_hit(qapp):
    d = Doc()
    for k, kind in enumerate(model.BIO_KINDS):
        d.add_shape(*biodraw.default_shape(kind, k * 10, 0), kind)
    r = render.Renderer(d)
    img = render.to_image(d, dpi=40)
    assert img.width() > 100
    for i, sh in d.shapes.items():
        assert not r.shape_paths(i).boundingRect().isEmpty()
        cx, cy = (sh.x1 + sh.x2) / 2, (sh.y1 + sh.y2) / 2
        assert d.find_shape(cx, cy, 0.1) == i                                      # the inside counts
    pl = d.add_shape(*biodraw.default_shape('plasmid', 0, 20), 'plasmid')
    d.shapes[pl].params = {'bp': 4361, 'name': 'pBR322', 'regions': [{'start': 3293, 'end': 4153, 'label': 'AmpR'}],
                           'markers': [{'pos': 29, 'label': 'EcoRI'}]}
    render.to_image(d, dpi=40)


def test_cdxml_skips_bio(tmp_path, qapp):
    import cdxml
    d = Doc()
    d.add_shape(0, 0, 3, 3, 'cell')
    d.add_shape(0, 5, 3, 7, 'rect')
    d.add_sequence(0, 10, 'PEPTIDE', list('AG'))
    out = tmp_path / 'x.cdxml'
    cdxml.write(d, str(out))
    back, _notes = cdxml.read(str(out))
    assert [s.kind for s in back.shapes.values()] == ['rect']


# ------------------------------------------------------------------------------------------------ templates
def test_bio_templates():
    for cat in ('Amino acids', 'Sugars', 'Nucleobases and nucleotides'):
        assert templates.TEMPLATES[cat]
        for name, smi in templates.TEMPLATES[cat]:
            assert Chem.MolFromSmiles(smi) is not None, name
    sugars = dict(templates.TEMPLATES['Sugars'])

    def cip(smi):
        m = Chem.MolFromSmiles(smi)
        rdCIPLabeler.AssignCIPLabels(m)
        return sorted(a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode'))
    assert cip(sugars['β-D-Glucopyranose']) == sorted('RRSSR')
    assert cip(sugars['α-D-Glucopyranose']) == sorted('SRSSR')
    assert cip(sugars['β-D-Fructofuranose']) == sorted('RSSR')
    nuc = dict(templates.TEMPLATES['Nucleobases and nucleotides'])
    assert D.CalcMolFormula(Chem.MolFromSmiles(nuc['ATP'])) == 'C10H16N5O13P3'
    adenosine = Chem.MolFromSmiles(nuc['Adenosine'])
    for key in ('ATP', 'NAD⁺'):                                                    # the same ribose stereo
        assert Chem.MolFromSmiles(nuc[key]).HasSubstructMatch(adenosine, useChirality=True), key
    assert 'Betalains' not in templates.TEMPLATES                                  # left out on request


# ------------------------------------------------------------------------------------------------ window
def test_window_sequence_actions(win, qapp):
    c = win.canvas
    c.set_doc(Doc())
    i = c.add_sequence(dict(polymer='PEPTIDE', residues=list('CAAC'), bridges=[[0, 3]]))
    c.set_selection([('seq', i)])
    win.info.refresh()
    assert 'Sequence' in win.info.which.currentText()
    win.copy_sequence('helm')
    assert QtWidgets.QApplication.clipboard().text().startswith('PEPTIDE1{C.A.A.C}$')
    c.set_sequence_props(three_letter=True)
    assert c.doc.seqs[i].three_letter
    QtWidgets.QApplication.clipboard().setText('H-Ala-Gly-Pro-NH2')
    win.paste_sequence('PEPTIDE')
    assert len(c.doc.seqs) == 2
    QtWidgets.QApplication.clipboard().setText('PEPTIDE1{A.C}$$$$')
    win.paste()
    assert len(c.doc.seqs) == 3
    c.set_selection([('seq', i)])
    win.sequence_analysis_caption()
    assert any('C12H20N4O5S2' in t.text for t in c.doc.texts.values())
    c.set_selection([('seq', i)])
    c.expand_sequences()
    assert i not in c.doc.seqs and len(c.doc.seqs) == 2
    assert len(c.doc.atoms) == 23                                                   # C12 N4 O5 S2 heavy atoms
    c.undo()
    assert i in c.doc.seqs
    for kind in model.BIO_KINDS:                                                    # every bio tool by a click
        c.set_tool(kind)
        if kind != 'plasmid':
            c.click(None, 5.0, 5.0)
    assert sum(sh.kind in model.BIO_KINDS for sh in c.doc.shapes.values()) == len(model.BIO_KINDS) - 1
