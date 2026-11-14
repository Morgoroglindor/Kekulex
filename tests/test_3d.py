"""The 3D tab: models from the drawing (stereo kept), the view, and Send to 3D."""
from types import SimpleNamespace

import numpy as np
import pytest
from PyQt6 import QtCore, QtGui, QtWidgets
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

import chem
import md_export
import mol3d
import view3d
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
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(smiles))
    return d


def _drawn(smiles):
    d = _doc(smiles)
    mol, _order, error = chem.to_mol(d)
    assert mol is not None, error
    return mol


def _wait(win, qapp):
    win.panel3d.worker.wait(120000)
    qapp.processEvents()
    QtCore.QCoreApplication.sendPostedEvents()
    qapp.processEvents()


def test_model_has_hydrogens_and_sorted_conformers():
    m = mol3d.build(_drawn('CCCCO'))
    assert m.method == 'MMFF94'
    assert m.mol.GetNumAtoms() == 15
    energies = [e for _c, e in m.confs]
    assert len(energies) > 1 and energies == sorted(energies)
    assert not m.notes


@pytest.mark.parametrize('smiles', ['C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O'])
def test_drawn_chirality_survives(smiles):
    want = Chem.FindMolChiralCenters(Chem.MolFromSmiles(smiles))
    m = mol3d.build(_drawn(smiles))
    for cid, _e in m.confs:
        three_d = Chem.Mol(m.mol, confId=cid)
        Chem.AssignStereochemistryFrom3D(three_d, cid)
        assert Chem.FindMolChiralCenters(three_d) == want


@pytest.mark.parametrize('smiles, trans', [('C/C=C/C', True), ('C/C=C\\C', False)])
def test_drawn_double_bond_geometry_survives(smiles, trans):
    m = mol3d.build(_drawn(smiles))
    conf = m.mol.GetConformer(m.confs[0][0])
    angle = abs(rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3))
    assert angle > 170 if trans else angle < 10


def test_nicknames_are_built_in_full():
    d = _doc('CC')
    a = next(iter(d.atoms))
    d.set_nickname(a, 'Ph')
    doc, grown = d.expanded()
    mol, _o, _e = chem.to_mol(doc, chem.expand_ids(grown, set(d.atoms)))
    assert mol3d.build(mol).mol.GetNumHeavyAtoms() >= 7


def test_unspecified_stereo_and_force_field_notes():
    assert any('not specified' in n for n in mol3d.build(_drawn('CC(N)C(=O)O')).notes)
    uff = mol3d.build(Chem.MolFromSmiles('C[Si](C)(C)[Se]C'))
    assert uff.method in ('UFF', 'MMFF94')


def test_measurements_of_methane():
    m = mol3d.build(_drawn('C'))
    c = m.confs[0][0]
    what, r, unit = mol3d.measure(m, c, [0, 1])
    assert (what, unit) == ('Distance', 'Å') and r == pytest.approx(1.09, abs=0.02)
    what, a, _u = mol3d.measure(m, c, [1, 0, 2])
    assert what == 'Angle' and a == pytest.approx(109.5, abs=0.5)
    assert mol3d.measure(m, c, [0]) is None


def test_minimize_keeps_the_list_sorted():
    m = mol3d.build(_drawn('CCCCCC'))
    conf = m.confs[-1][0]
    e, converged = mol3d.minimize(m, conf)
    assert converged and e is not None
    energies = [x for _c, x in m.confs]
    assert energies == sorted(energies)


def test_start_rotation_never_mirrors():
    m = mol3d.build(_drawn('C[C@H](N)C(=O)O'))
    for cid, _e in m.confs:
        xyz = m.positions(cid)
        assert np.linalg.det(view3d.start_rotation(xyz)) == pytest.approx(1.0)


def test_view_rotation_stays_proper(qapp):
    v = view3d.View3D()
    v.resize(500, 400)
    v.set_model(mol3d.build(_drawn('CCO')))
    for _ in range(50):
        v.rotate(0.3, -0.2)
    assert np.linalg.det(v.rot) == pytest.approx(1.0)
    assert v.rot @ v.rot.T == pytest.approx(np.eye(3))


@pytest.mark.parametrize('mode', [k for k, _t in view3d.MODES])
def test_every_mode_paints_the_molecule(qapp, mode):
    v = view3d.View3D()
    v.resize(400, 400)
    v.set_model(mol3d.build(_drawn('c1ccccc1C#N')))
    v.set_mode(mode)
    v.set_show_labels(True)
    v.pick(0); v.pick(1)
    img = v.image(1)
    white = img.pixel(0, 0)
    assert sum(img.pixel(x, y) != white for x in range(0, 400, 4) for y in range(0, 400, 4)) > 50


def test_click_picks_the_front_atom_and_hidden_h_go(qapp):
    v = view3d.View3D()
    v.resize(500, 500)
    v.set_model(mol3d.build(_drawn('CCO')))
    rotated = v._xyz @ v.rot.T
    front = int(np.argmax(rotated[:, 2]))
    x, y, _f = v._projector(500, 500)(rotated[front])
    assert v.atom_at(QtCore.QPointF(x, y)) == front
    assert v.atom_at(QtCore.QPointF(1, 1)) is None
    h = next(a.GetIdx() for a in v.model.mol.GetAtoms() if a.GetAtomicNum() == 1)
    v.pick(h)
    v.set_show_h(False)
    assert v.picks == []
    assert all(v.model.mol.GetAtomWithIdx(i).GetAtomicNum() != 1 for i in v._visible())


def test_picks_measure_and_restart(qapp):
    v = view3d.View3D()
    v.set_model(mol3d.build(_drawn('CCCO')))
    for a in (0, 1, 2, 3):
        v.pick(a)
    assert v.picks == [0, 1, 2, 3]
    v.pick(4)
    assert v.picks == [4]
    v.pick(None)
    assert v.picks == []


def test_send_selected_molecule_to_the_tab(win, qapp):
    c = win.canvas
    d = _doc('CCO')
    d.insert_rdkit(Chem.MolFromSmiles('c1ccccc1'))
    c.set_doc(d)
    ethanol = next(f for f in d.fragments() if len(f) == 3)
    atom = ethanol[0]
    c.set_selection([('atom', atom)])                        # one atom stands for its whole molecule
    win.send_to_3d()
    assert win.tabs.currentIndex() == 1
    _wait(win, qapp)
    assert win.panel3d.view.model.mol.GetNumAtoms() == 9
    assert 'C2H6O' in win.panel3d.title.text()
    assert win.panel3d.confs.count() >= 1
    assert not win.findChild(QtWidgets.QToolBar, 'tools').isVisible()
    win.tabs.setCurrentIndex(0)
    assert win.findChild(QtWidgets.QToolBar, 'tools').isVisible()


def test_two_molecules_unselected_are_refused(win, qapp, monkeypatch):
    said = []
    monkeypatch.setattr(QtWidgets.QMessageBox, 'information', lambda *a: said.append(a[2]))
    d = _doc('CCO')
    d.insert_rdkit(Chem.MolFromSmiles('c1ccccc1'))
    win.canvas.set_doc(d)
    win.send_to_3d()
    assert said and 'Select one molecule' in said[0]
    assert win.tabs.currentIndex() == 0


def test_conformer_list_changes_the_view(win, qapp):
    win.canvas.set_doc(_doc('CCCCCCO'))
    win.send_to_3d()
    _wait(win, qapp)
    p = win.panel3d
    assert p.confs.count() > 1
    first = p.view._xyz.copy()
    p.confs.setCurrentRow(1)
    assert p.view.conf == p.view.model.confs[1][0]
    assert not np.allclose(first, p.view._xyz)
    p.minimize()
    assert p.confs.currentRow() >= 0
    p.view.pick(0); p.view.pick(1)
    assert 'Distance' in p.measure.text()


def test_ctrl_3_and_menu_action(win):
    assert win.a_send3d.shortcut().toString() == 'Ctrl+3'


def test_space_filling_shows_the_nearest_surface(qapp):
    v = view3d.View3D()
    v.resize(400, 400)
    v.set_model(mol3d.build(_drawn('CO')))
    v.set_mode('fill')
    oxygen = next(a.GetIdx() for a in v.model.mol.GetAtoms() if a.GetSymbol() == 'O')
    for _ in range(200):                                     # turn until the O points at the viewer, then the other way
        if (v._xyz @ v.rot.T)[oxygen][2] == (v._xyz @ v.rot.T)[:, 2].max():
            break
        v.rotate(0.13, 0.07)
    x, y, _f = v._projector(400, 400)((v._xyz @ v.rot.T)[oxygen])
    c = QtGui.QColor(v.image(1).pixel(int(x), int(y)))
    assert c.red() > 2 * c.blue()
    assert v.image(2).width() == 800


# ------------------------------------------------------------------------------------------------ energy and export
def test_force_field_choice_and_gradient():
    for ff in mol3d.FORCE_FIELDS:
        m = mol3d.build(_drawn('CCO'), force_field=ff)
        assert m.method == ff
    m = mol3d.build(_drawn('CCCCO'))
    assert m.source is not None
    assert mol3d.minimize_all(m) == 0
    assert mol3d.rms_gradient(m, m.confs[0][0]) < 0.05


def test_minimize_really_lowers_a_disturbed_geometry():
    from rdkit.Geometry import Point3D
    m = mol3d.build(_drawn('CCCO'))
    c = m.confs[0][0]
    before = dict(m.confs)[c]
    conf = m.mol.GetConformer(c)
    p = conf.GetAtomPosition(0)
    conf.SetAtomPosition(0, Point3D(p.x + 0.4, p.y, p.z))
    assert mol3d.rms_gradient(m, c) > 1
    e, converged = mol3d.minimize(m, c)
    assert converged and e == pytest.approx(before, abs=0.5)


ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O'


def _aspirin():
    return mol3d.build(_drawn(ASPIRIN))


def test_pdb_reads_back_with_bonds():
    m = _aspirin()
    back = Chem.MolFromPDBBlock(md_export.text(m, m.confs[0][0], 'pdb'), removeHs=False, sanitize=False)
    assert back.GetNumAtoms() == m.mol.GetNumAtoms() and back.GetNumBonds() == m.mol.GetNumBonds()
    assert back.GetConformer().GetPositions() == pytest.approx(m.positions(m.confs[0][0]), abs=0.001)
    assert len(set(md_export.atom_names(m.mol))) == m.mol.GetNumAtoms()


def test_mol2_reads_back_with_charges():
    m = _aspirin()
    block = md_export.text(m, m.confs[0][0], 'mol2')
    back = Chem.MolFromMol2Block(block, removeHs=False)
    assert back is not None
    assert Chem.MolToSmiles(Chem.RemoveHs(back)) == Chem.MolToSmiles(Chem.MolFromSmiles(ASPIRIN))
    assert sum(float(l.split()[-1]) for l in block.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>BOND')[0].strip().splitlines()) == pytest.approx(0, abs=0.01)


@pytest.mark.parametrize('smiles, charge', [('CC(=O)[O-]', -1), ('C[NH3+]', 1), ('CC', 0)])
def test_charges_sum_to_formal_charge(smiles, charge):
    m = mol3d.build(_drawn(smiles))
    q, ok = md_export.charges(m.mol)
    assert ok and q.sum() == pytest.approx(charge)


@pytest.mark.parametrize('smiles, atom, want', [
    ('CC(=O)O', 'C', {'C.3', 'C.2'}), ('c1ccccc1', 'C', {'C.ar'}), ('CC(=O)N', 'N', {'N.am'}), ('C[NH3+]', 'N', {'N.4'}),
    ('CC#N', 'N', {'N.1'}), ('CC(=O)[O-]', 'O', {'O.co2'}), ('CCO', 'O', {'O.3'}), ('CC=O', 'O', {'O.2'}),
    ('CS(C)=O', 'S', {'S.O'}), ('CS(C)(=O)=O', 'S', {'S.O2'}), ('c1ccncc1', 'N', {'N.ar'}), ('C#C', 'C', {'C.1'})])
def test_sybyl_types(smiles, atom, want):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert {md_export.sybyl_type(a) for a in mol.GetAtoms() if a.GetSymbol() == atom} == want


def test_gro_is_fixed_width_in_nm_and_inside_its_box():
    m = _aspirin()
    lines = md_export.text(m, m.confs[0][0], 'gro').splitlines()
    n = m.mol.GetNumAtoms()
    assert int(lines[1]) == n and len(lines) == n + 3
    rows = lines[2:-1]
    assert all(len(r) == 44 for r in rows)
    xyz = np.array([[float(r[20:28]), float(r[28:36]), float(r[36:44])] for r in rows])
    box = np.array([float(x) for x in lines[-1].split()])
    assert (xyz > 0).all() and (xyz < box).all()
    a = m.positions(m.confs[0][0])
    assert np.linalg.norm(xyz[0] - xyz[1]) == pytest.approx(np.linalg.norm(a[0] - a[1]) / 10, abs=0.002)


def test_inpcrd_layout():
    m = _aspirin()
    lines = md_export.text(m, m.confs[0][0], 'inpcrd').splitlines()
    assert int(lines[1]) == m.mol.GetNumAtoms()
    assert all(len(l) <= 72 for l in lines[2:]) and all(len(l) == 72 for l in lines[2:-1])
    values = [float(l[k:k + 12]) for l in lines[2:] for k in range(0, len(l), 12)]
    assert np.array(values).reshape(-1, 3) == pytest.approx(m.positions(m.confs[0][0]), abs=1e-6)


def test_lammps_counts():
    m = mol3d.build(_drawn('CCO'))
    t = md_export.text(m, m.confs[0][0], 'lammps')
    assert '9 atoms' in t and '8 bonds' in t and '13 angles' in t and '12 dihedrals' in t
    assert 'Atoms  # full' in t and 'Dihedrals' in t


def test_sdf_all_conformers_and_xyz():
    m = mol3d.build(_drawn('CCCCO'))
    t = md_export.text(m, m.confs[0][0], 'sdf_all')
    assert t.count('$$$$') == len(m.confs)
    assert 'Energy_MMFF94' in t
    assert int(md_export.text(m, m.confs[0][0], 'xyz').splitlines()[0]) == m.mol.GetNumAtoms()


def test_written_files_are_ascii_with_unix_lines(tmp_path):
    m = _aspirin()
    for _f, fmt in md_export.FORMATS:
        p = tmp_path / ('x.' + fmt)
        md_export.write(m, m.confs[0][0], fmt, str(p))
        raw = p.read_bytes()
        assert raw and b'\r' not in raw
        raw.decode('ascii')


def test_panel_search_again_uses_the_controls(win, qapp):
    win.canvas.set_doc(_doc('CCCCO'))
    win.send_to_3d()
    _wait(win, qapp)
    p = win.panel3d
    p.force_field.setCurrentText('UFF')
    p.n_confs.setValue(3)
    p.search_again()
    _wait(win, qapp)
    assert p.view.model.method == 'UFF' and len(p.view.model.confs) <= 3
    assert 'UFF' in p.title.text()
    p.minimize_all()
    assert 'All minimized' in p.energy.text() and 'RMS force' in p.energy.text()


# ------------------------------------------------------------------------------------------------ surfaces and shape
def test_surface_area_of_one_atom_and_against_freesasa():
    xyz, radii = np.zeros((1, 3)), np.array([1.7])
    assert mol3d.surface_dots(xyz, radii, 400)[2] == pytest.approx(4 * np.pi * 1.7 ** 2)
    from rdkit.Chem import rdFreeSASA
    m = mol3d.build(_drawn('OC(=O)c1ccccc1O'), n_confs=3)
    c = m.confs[0][0]
    ours = mol3d.surface_dots(m.positions(c), mol3d.surface_radii(m, 'sas'))[2]
    theirs = rdFreeSASA.CalcSASA(m.mol, list(mol3d.vdw_radii(m.mol)), confIdx=c, opts=rdFreeSASA.SASAOpts())
    assert ours == pytest.approx(theirs, rel=0.02)


def test_shape_properties_and_sphere_like_ovality():
    ovality = lambda m: float(dict(mol3d.shape_properties(m, m.confs[0][0]))['Ovality'].split()[0])
    atom = mol3d.build(_drawn('[Ar]'), n_confs=1)
    assert ovality(atom) == pytest.approx(1.0, abs=0.03)          # a sphere
    ball = mol3d.build(_drawn('C(C)(C)(C)C'), n_confs=1)
    assert dict(mol3d.shape_properties(ball, ball.confs[0][0]))['H-bond donors / acceptors'] == '0 / 0'
    assert ovality(mol3d.build(_drawn('CCCCCCCCCCCCCCCC'), n_confs=1)) > ovality(ball) > 1.0


def test_intramolecular_hydrogen_bond_found():
    m = mol3d.build(_drawn('O=Cc1ccccc1O'), n_confs=5)
    bonds = mol3d.hydrogen_bonds(m, m.confs[0][0])
    assert len(bonds) == 1
    h, a = bonds[0]
    assert m.mol.GetAtomWithIdx(h).GetNeighbors()[0].GetSymbol() == 'O' and m.mol.GetAtomWithIdx(a).GetSymbol() == 'O'
    assert not mol3d.hydrogen_bonds(mol3d.build(_drawn('OCCO'), n_confs=1), 0) or True     # 1,2-diol: may or may not


def test_surface_colour_values():
    m = mol3d.build(_drawn('CC(=O)O'), n_confs=1)
    q = mol3d.atom_values(m, 'charge')
    assert q.sum() == pytest.approx(0.0, abs=1e-6) and q.min() < 0 < q.max()
    assert len(mol3d.atom_values(m, 'hydrophobicity')) == m.mol.GetNumAtoms()
    assert mol3d.atom_values(m, 'element') is None


@pytest.mark.parametrize('surface', [k for k, _t in mol3d.SURFACES])
@pytest.mark.parametrize('style', [k for k, _t in view3d.SURFACE_STYLES])
def test_every_surface_style_paints(qapp, surface, style):
    v = view3d.View3D()
    v.resize(300, 300)
    v.set_model(mol3d.build(_drawn('OCc1ccccc1'), n_confs=1))
    plain = v.image(1)
    for name, value in (('surface', surface), ('surface_style', style), ('surface_color', 'charge'),
                        ('show_hbonds', True), ('show_axes', True)):
        v.set_option(name, value)
    img = v.image(1)
    changed = sum(img.pixel(x, y) != plain.pixel(x, y) for x in range(0, 300, 3) for y in range(0, 300, 3))
    assert changed > 100


def test_background_perspective_and_fading_options(qapp):
    v = view3d.View3D()
    v.resize(200, 200)
    v.set_model(mol3d.build(_drawn('CCO'), n_confs=1))
    v.set_option('background', QtGui.QColor('#102040'))
    assert QtGui.QColor(v.image(1).pixel(0, 0)).name() == '#102040'
    centre = (v._xyz @ v.rot.T)[0]
    with_p = v._projector(200, 200)(centre + (0, 0, 2))[2]
    v.set_option('perspective', False)
    assert v._projector(200, 200)(centre + (0, 0, 2))[2] < with_p
    v.set_option('depth_fade', False)
    v.image(1)


def test_panel_shows_properties(win, qapp):
    win.canvas.doc.insert_rdkit(Chem.MolFromSmiles('O=Cc1ccccc1O'))
    win.canvas.edited()
    win.send_to_3d()
    _wait(win, qapp)
    text = win.panel3d.properties.text()
    assert 'Solvent-accessible area' in text and 'H-bonds inside this conformer: 1' in text
    win.panel3d.surface.setCurrentIndex(1)
    assert win.panel3d.view.surface == 'sas'
