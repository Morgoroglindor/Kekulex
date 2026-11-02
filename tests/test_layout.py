"""Structures from a name or SMILES are drawn the way chemists draw them, without changing the molecule."""
import numpy as np
import pytest
from rdkit import Chem

import chem
import layout2d
from model import Doc

MOLS = {'cholesterol': 'CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
        'testosterone': 'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C',
        'estradiol': 'C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc34)[C@@H]1CC[C@@H]2O',
        'morphine': 'CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5',
        'penicillin G': 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
        'cephalexin': 'CC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)[C@@H](c3ccccc3)N)SC1)C(=O)O',
        'glucose': 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
        'adenosine': 'Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O',
        'quinine': 'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O',
        'tryptophan': 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N',
        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'quercetin': 'O=C1C(O)=C(c2ccc(O)c(O)c2)Oc2cc(O)cc(O)c12',
        'E-alkene': 'C/C=C/C(=O)O', 'Z-alkene': 'C/C=C\\C(=O)O', 'L-alanine': 'C[C@H](N)C(=O)O', 'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
        'salt': 'CC(=O)[O-].[Na+]', 'atom': 'O'}


def _xy(mol, smarts):
    conf = mol.GetConformer()
    return [np.array([conf.GetAtomPosition(a).x, conf.GetAtomPosition(a).y]) for a in mol.GetSubstructMatch(Chem.MolFromSmarts(smarts))]


def _laid_out(smiles):
    m = Chem.MolFromSmiles(smiles)
    layout2d.compute(m)
    return m


@pytest.mark.parametrize('name', MOLS)
def test_the_drawing_is_still_the_same_molecule(name, qapp):
    """Mirroring and turning the layout must not change R/S or E/Z, through the wedges of the real drawing."""
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(MOLS[name]))
    mol, _order, error = chem.to_mol(d)
    assert mol is not None, error
    Chem.DetectBondStereoChemistry(mol, mol.GetConformer())
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(MOLS[name]))


@pytest.mark.parametrize('name', ['testosterone', 'cholesterol', 'estradiol'])
def test_steroid_ring_a_lower_left_and_ring_d_upper_right(name):
    m = _laid_out(MOLS[name])
    assert layout2d.template_for(m)[0] == 'steroid'
    query, _coords = next((q, c) for n, q, c in layout2d._compiled if n == 'steroid')
    conf = m.GetConformer()
    match = m.GetSubstructMatch(query)
    p = {k: conf.GetAtomPosition(match[k]) for k in (2, 14)}           # C3 and C16
    assert p[2].x < p[14].x and p[2].y < p[14].y


def test_steroid_side_chain_points_outwards():
    assert layout2d.clashes(_laid_out(MOLS['cholesterol'])) == 0


def test_purine_six_ring_left_and_imidazole_right():
    m = _laid_out(MOLS['caffeine'])
    n1, _c2, _n3, _c4, _c5, _c6, _n7, c8, _n9 = _xy(m, '[#7]1~[#6]~[#7]~[#6]2~[#6](~[#6]~1)~[#7]~[#6]~[#7]~2')
    assert n1[0] < c8[0]


def test_indole_nh_at_the_lower_right():
    m = _laid_out(MOLS['tryptophan'])
    assert layout2d.template_for(m)[0] == 'indole'
    ring = _xy(m, '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#7]~2')
    assert ring[8][1] < ring[6][1] and ring[8][0] > ring[0][0]


@pytest.mark.parametrize('smiles, ok', [('CCCCCCCCCC(=O)O', True), ('OC(=O)CCCCCCCCC', True)])
def test_chain_acid_on_the_right(smiles, ok):
    m = _laid_out(smiles)
    conf = m.GetConformer()
    o = [conf.GetAtomPosition(a.GetIdx()).x for a in m.GetAtoms() if a.GetSymbol() == 'O']
    c = [conf.GetAtomPosition(a.GetIdx()).x for a in m.GetAtoms() if a.GetSymbol() == 'C']
    assert (min(o) > np.mean(c)) == ok


def test_every_template_is_consistent():
    for name, query, coords in layout2d._compiled:
        assert query is not None and query.GetNumAtoms() == len(coords), name
        xy = np.array(coords)
        d = np.linalg.norm(xy[:, None] - xy[None], axis=2)
        np.fill_diagonal(d, 9)
        assert d.min() > 0.3, name
        lengths = [d[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in query.GetBonds()]
        if name not in ('morphinan', 'quinuclidine'):
            assert max(lengths) < 1.35 and min(lengths) > 0.7, name


def test_a_template_only_needs_its_skeleton():
    assert layout2d.template_for(Chem.MolFromSmiles('c1ccccc1')) is None
    assert layout2d.template_for(Chem.MolFromSmiles('CCO')) is None


BETALAINS = {'betalamic acid': 'O=CC=C1C[C@@H](C(=O)O)NC(C(=O)O)=C1',
             'betanidin': 'OC(=O)[C@@H]1Cc2cc(O)c(O)cc2[N+]1=CC=C1C[C@@H](C(=O)O)NC(C(=O)[O-])=C1',
             'betanin': 'O=C([O-])[C@@H]1Cc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc2[N+]1=CC=C1C[C@@H](C(=O)O)NC(C(=O)O)=C1',
             'indicaxanthin': 'OC(=O)[C@@H]1CCC[N+]1=CC=C1C[C@@H](C(=O)O)NC(C(=O)[O-])=C1',
             'vulgaxanthin I': 'NC(=O)CC[C@H](NC=CC1=CC(C(=O)O)=NC(C(=O)O)C1)C(=O)O',
             'miraxanthin V': 'OC(=O)C1CC(=CC=[NH+]CCc2ccc(O)c(O)c2)C=C(N1)C(=O)[O-]',
             'portulacaxanthin II': 'OC(=O)C(Cc1ccc(O)cc1)[NH+]=CC=C1CC(C(=O)O)NC(C(=O)[O-])=C1',
             'neobetanin': 'O=C([O-])[C@@H]1Cc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc2[N+]1=CC=C1C=C(C(=O)O)N=C(C(=O)O)C1'}


@pytest.mark.parametrize('name', BETALAINS)
def test_betalains_ring_n_at_the_bottom_and_the_chain_up(name, qapp):
    m = _laid_out(BETALAINS[name])
    assert layout2d.template_for(m)[0] == 'betalamic'
    assert layout2d.clashes(m) == 0
    query = next(q for n, q, _c in layout2d._compiled if n == 'betalamic')
    match = m.GetSubstructMatch(query)
    conf = m.GetConformer()
    y = {k: conf.GetAtomPosition(match[k]).y for k in (0, 7)}    # the ring N and the chain carbon
    assert y[0] < y[7]
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(BETALAINS[name]))
    mol, _order, error = chem.to_mol(d)
    assert mol is not None, error
    Chem.DetectBondStereoChemistry(mol, mol.GetConformer())
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for b in mol.GetBonds():                                 # these SMILES leave E/Z open; a drawing always shows one
        b.SetStereo(Chem.BondStereo.STEREONONE)
        b.SetBondDir(Chem.BondDir.NONE)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(BETALAINS[name]))


def test_betalains_stay_out_of_the_prepacked_molecules():
    import templates
    names = ' '.join(n.lower() for items in templates.TEMPLATES.values() for n, _s in items)
    assert not any(word in names for word in ('betal', 'betan', 'xanthin'))


# from PubChem by name
TRITERPENES = {
    "oleanolic acid": "C[C@]12CC[C@@H](C([C@@H]1CC[C@@]3([C@@H]2CC=C4[C@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)C(=O)O)C)C)(C)C)O",
    "ursolic acid": "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H](C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O",
    "beta-amyrin": "C[C@@]12CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H](C5(C)C)O)C)C)[C@@H]1CC(CC2)(C)C)C",
    "alpha-amyrin": "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CC[C@@H](C5(C)C)O)C)C)[C@@H]2[C@H]1C)C)C",
    "friedelin": "C[C@H]1C(=O)CC[C@@H]2[C@@]1(CC[C@H]3[C@]2(CC[C@@]4([C@@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)C)C)C)C)C",
    "glycyrrhetinic acid": "C[C@]12CC[C@](C[C@H]1C3=CC(=O)[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)(C)C(=O)O",
    "taraxasterol": "C[C@H]1[C@@H]2[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC[C@]2(CCC1=C)C)C)C)(C)C)O)C",
    "maslinic acid": "C[C@@]12CC[C@@H]3[C@@]([C@H]1CC=C4[C@]2(CC[C@@]5([C@H]4CC(CC5)(C)C)C(=O)O)C)(C[C@H]([C@@H](C3(C)C)O)O)C",
    "asiatic acid": "C[C@@H]1CC[C@@]2(CC[C@@]3(C(=CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(C[C@H]([C@@H]([C@@]5(C)CO)O)O)C)C)[C@@H]2[C@H]1C)C)C(=O)O",
    "betulin": "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)CO",
    "betulinic acid": "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C(=O)O",
    "lupeol": "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C",
    "hopane": "CC(C)[C@H]1CC[C@]2([C@H]1CC[C@@]3([C@@H]2CC[C@H]4[C@]3(CC[C@@H]5[C@@]4(CCCC5(C)C)C)C)C)C"
}
TETRACYCLIC = {
    "lanosterol": "C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C",
    "cycloartenol": "C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CC[C@]34[C@H]2CC[C@@H]5[C@]3(C4)CC[C@@H](C5(C)C)O)C)C",
    "dammarenediol II": "CC(=CCC[C@@](C)([C@H]1CC[C@@]2([C@@H]1CC[C@H]3[C@]2(CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C)O)C",
    "euphol": "C[C@H](CCC=C(C)C)[C@@H]1CC[C@]2([C@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C"
}

PUBCHEM_BETANIN = "C1[C@H](N=C(C=C1/C=C/N2[C@@H](CC3=CC(=C(C=C32)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)C(=O)O)C(=O)O)C(=O)O"
PUBCHEM_NEOBETANIN = "C1[C@H](N(C2=CC(=C(C=C21)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O)O)/C=C/C4=CC(=NC(=C4)C(=O)O)C(=O)O)C(=O)O"
PUBCHEM_INDICAXANTHIN = "C1C[C@H]([N+](=C/C=C/2\\C[C@H](NC(=C2)C(=O)O)C(=O)O)C1)C(=O)[O-]"


@pytest.mark.parametrize('name', TRITERPENES)
def test_pentacyclic_triterpene_rises_from_ring_a_to_ring_e(name, qapp):
    m = _laid_out(TRITERPENES[name])
    found = layout2d.template_for(m)
    assert found[0] in ('oleanane', 'lupane')
    assert layout2d.clashes(m) == 0
    query = next(q for n, q, _c in layout2d._compiled if n == found[0])
    match = next(k for k in (found[1],))
    conf = m.GetConformer()
    atoms = sorted(match, key=lambda a: match[a].x)
    left, right = conf.GetAtomPosition(atoms[0]), conf.GetAtomPosition(atoms[-1])
    assert left.x < right.x and left.y < right.y


@pytest.mark.parametrize('name', ['oleanolic acid', 'glycyrrhetinic acid', 'taraxasterol', 'maslinic acid', 'asiatic acid', 'friedelin'])
def test_ring_a_is_the_oxygenated_end(name):
    m = _laid_out(TRITERPENES[name])
    conf = m.GetConformer()
    oxygen = [a for a in m.GetAtoms() if a.GetAtomicNum() == 8 and any(n.IsInRing() and n.GetDegree() >= 3 for n in a.GetNeighbors())]
    fixed = layout2d.template_for(m)[1]
    c3 = next(a for a in fixed if abs(fixed[a].x + layout2d.R3) < 1e-6 and abs(fixed[a].y + 0.5) < 1e-6)
    assert any(n.GetAtomicNum() == 8 for n in m.GetAtomWithIdx(c3).GetNeighbors())
    assert oxygen


@pytest.mark.parametrize('name', TETRACYCLIC)
def test_tetracyclic_triterpenes_follow_the_steroid_skeleton(name):
    m = _laid_out(TETRACYCLIC[name])
    assert layout2d.template_for(m)[0] == 'steroid' and layout2d.clashes(m) == 0


@pytest.mark.parametrize('smiles', list(TRITERPENES.values()) + list(TETRACYCLIC.values()), ids=list(TRITERPENES) + list(TETRACYCLIC))
def test_triterpene_drawing_is_the_same_molecule(smiles, qapp):
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(smiles))
    mol, _order, error = chem.to_mol(d)
    assert mol is not None, error
    Chem.DetectBondStereoChemistry(mol, mol.GetConformer())
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_betanin_from_pubchem_gets_its_ring_nh():
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    new = layout2d.usual_form(PUBCHEM_BETANIN)
    before, after = Chem.MolFromSmiles(PUBCHEM_BETANIN), Chem.MolFromSmiles(new)
    ring_nh = Chem.MolFromSmarts('[#7;R;H1;X3]1-[#6;R](-C(=O)O)=[#6;R]-[#6;R](=[#6;!R]-[#6;!R]=[#7+])-[#6;R]-[#6;R]-1')
    assert not before.HasSubstructMatch(ring_nh) and after.HasSubstructMatch(ring_nh)
    assert CalcMolFormula(before) == CalcMolFormula(after)
    assert Chem.GetFormalCharge(after) == 0
    assert after.HasSubstructMatch(Chem.MolFromSmarts('[#7+]-[#6](-C(=O)[O-])'))    # the carboxylate next to the iminium N
    centres = lambda m: Chem.FindMolChiralCenters(m, useLegacyImplementation=False)
    assert len(centres(before)) == len(centres(after))
    assert [c for _a, c in centres(before)].count('S') + [c for _a, c in centres(before)].count('R') == len(centres(after))


def test_usual_form_keeps_the_extended_chain():
    # ChemSpider's iminium record of betanin (44750964), with the stereocentre it leaves open filled in
    expected = 'O=C(O)C1=C/C(=C/C=[N+]2/c3cc(O)c(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)cc3C[C@H]2C(=O)[O-])C[C@@H](C(=O)O)N1'
    assert Chem.CanonSmiles(layout2d.usual_form(PUBCHEM_BETANIN)) == Chem.CanonSmiles(expected)


@pytest.mark.parametrize('smiles', [layout2d.usual_form(PUBCHEM_BETANIN), PUBCHEM_NEOBETANIN,
                                    'O=C(O)C1=C/C(=C/C=[N+]2/CCC[C@H]2C(=O)[O-])C[C@@H](C(=O)O)N1',
                                    'C/C=C/c1ccccc1', 'C/C=C\\c1ccccc1', 'OC(=O)/C=C/C(=O)O', 'OC(=O)/C=C\\C(=O)O'])
def test_drawing_keeps_the_stated_double_bond_geometry(smiles, qapp):
    d = Doc()                                       # compared as SMILES: the InChIKey of the iminium forms ignores E/Z
    d.insert_rdkit(Chem.MolFromSmiles(smiles))
    drawn = Chem.RemoveHs(Chem.MolFromMolBlock(d.to_molblock()[0]))
    assert Chem.MolToSmiles(drawn) == Chem.CanonSmiles(smiles)


@pytest.mark.parametrize('smiles', [PUBCHEM_INDICAXANTHIN, PUBCHEM_NEOBETANIN, layout2d.usual_form(PUBCHEM_BETANIN),
                                    'O=C(O)C1=C/C(=C/C=[N+]2/CCC[C@H]2C(=O)[O-])C[C@@H](C(=O)O)N1'])
def test_betalains_are_turned_over_whatever_the_source(smiles, qapp):
    m = _laid_out(smiles)
    query = next(q for n, q, _c in layout2d._compiled if n == 'betalamic')
    match = m.GetSubstructMatch(query)
    conf = m.GetConformer()
    assert conf.GetAtomPosition(match[1]).x > conf.GetAtomPosition(match[9]).x    # the ring's C=C side on the right
    assert conf.GetAtomPosition(match[0]).y < conf.GetAtomPosition(match[7]).y    # NH still at the bottom


def test_usual_form_leaves_the_others_alone():
    for smiles in (PUBCHEM_NEOBETANIN, PUBCHEM_INDICAXANTHIN, 'CCO', 'c1ccncc1', 'CN=CC=CN(C)C'):
        assert layout2d.usual_form(smiles) == smiles


def _bond_lengths(mol):
    conf = mol.GetConformer()
    return [layout2d._length(conf, b) for b in mol.GetBonds()]


@pytest.mark.parametrize('smiles', [TRITERPENES['lupeol'], TRITERPENES['oleanolic acid'], BETALAINS['betanin'], BETALAINS['indicaxanthin'],
                                    TETRACYCLIC['lanosterol'], MOLS['testosterone'], MOLS['aspirin'], MOLS['caffeine'], MOLS['glucose']],
                         ids=['lupeol', 'oleanolic', 'betanin', 'indicaxanthin', 'lanosterol', 'testosterone', 'aspirin', 'caffeine', 'glucose'])
def test_bonds_have_one_length(smiles):
    m = _laid_out(smiles)
    lengths = _bond_lengths(m)
    assert max(lengths) / min(lengths) < 1.12


def test_lupane_rings_are_regular():
    m = _laid_out(TRITERPENES['lupeol'])
    fixed = layout2d.template_for(m)[1]
    conf = m.GetConformer()
    for a, p in fixed.items():
        q = conf.GetAtomPosition(a)
        assert (q.x, q.y) == pytest.approx((p.x, p.y), abs=1e-6)
    for b in m.GetBonds():
        if b.GetBeginAtomIdx() in fixed and b.GetEndAtomIdx() in fixed:
            assert layout2d._length(conf, b) == pytest.approx(1.0, abs=0.01)


def test_ring_fusion_stereocentres_show_their_h_and_no_ring_wedges(qapp):
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(TRITERPENES['lupeol']))
    h = [a for a in d.atoms.values() if a.el == 'H']
    assert len(h) == 4                                       # C5, C9, C13, C18, as in the usual drawing
    ring_bonds = [b for i, b in d.bonds.items() if d.ring_of_bond(i)]
    assert not [b for b in ring_bonds if b.stereo]
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    assert chem.analyse(d).formula == 'C30H50O'


def test_explicit_h_stays_out_of_nmr_trouble(qapp):
    import nmr
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(MOLS['testosterone']))
    r = nmr.predict(d, None, '13C', 'cdcl3')
    assert sum(len(s.atoms) for s in r.signals) == 19
    r = nmr.predict(d, None, '1H', 'cdcl3')
    assert sum(s.n_h for s in r.signals) == 28


@pytest.mark.parametrize('name', ['betanin', 'indicaxanthin', 'betanidin'])
def test_betalain_carboxylate_right_of_the_iminium_as_in_betanin(name, qapp):
    m = _laid_out(BETALAINS[name])
    match = m.GetSubstructMatch(Chem.MolFromSmarts('[#6;!R]=[#7+;R]-[#6;R]-[#6](~[#8])~[#8]'))
    conf = m.GetConformer()
    assert conf.GetAtomPosition(match[2]).x > conf.GetAtomPosition(match[1]).x
    assert layout2d.clashes(m) == 0
    lengths = _bond_lengths(m)
    assert max(lengths) / min(lengths) < 1.12
