"""Drawing geometry: rings, fusion, chains, bond directions, double-bond sides, files."""
import json
import math

import chem
from model import Doc


def test_benzene_alternates():
    doc = Doc()
    doc.add_ring(doc.ring_points(6, (0, 0), -math.pi / 2), aromatic=True)
    assert sorted(b.order for b in doc.bonds.values()) == [1, 1, 1, 2, 2, 2]
    assert chem.analyse(doc).formula == 'C6H6'


def test_fused_benzene_is_naphthalene():
    doc = Doc()
    doc.add_ring(doc.ring_points(6, (0, 0), -math.pi / 2), aromatic=True)
    right = max(doc.bonds, key=lambda i: doc.atoms[doc.bonds[i].a].x + doc.atoms[doc.bonds[i].b].x)
    doc.add_ring(doc.ring_on_bond(right, 6), aromatic=True)
    assert len(doc.atoms) == 10 and len(doc.bonds) == 11
    info = chem.analyse(doc)
    assert info.formula == 'C10H8' and info.smiles == 'c1ccc2ccccc2c1'


def test_fused_ring_goes_outside():
    doc = Doc()
    doc.add_ring(doc.ring_points(6, (0, 0), -math.pi / 2))
    b = next(iter(doc.bonds))
    pts = doc.ring_on_bond(b, 5)
    cx = sum(x for x, _y in pts) / 5; cy = sum(y for _x, y in pts) / 5
    assert math.hypot(cx, cy) > 1.0                       # its centre is outside the first ring


def test_zigzag_growth():
    doc = Doc()
    a = doc.add_atom(0, 0)
    prev = a
    for _ in range(4):
        ang = doc.next_bond_angle(prev)
        at = doc.atoms[prev]
        nxt = doc.atom_at_or_new(at.x + math.cos(ang), at.y + math.sin(ang), exclude={prev})
        doc.add_bond(prev, nxt)
        prev = nxt
    xs = [doc.atoms[i].x for i in sorted(doc.atoms)]
    assert xs == sorted(xs)                                # a chain growing to the right, not curling up
    assert chem.analyse(doc).formula == 'C5H12'


def test_chain_points():
    doc = Doc()
    pts = doc.chain_points(0, 0, 4 * math.sqrt(3) / 2, 0)
    assert len(pts) == 4
    for (x1, y1), (x2, y2) in zip([(0, 0)] + pts, pts):
        assert math.hypot(x2 - x1, y2 - y1) == __import__('pytest').approx(1.0)


def test_double_bond_sides():
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('CC(C)=O'))
    co = next(i for i, b in doc.bonds.items() if b.order == 2)
    assert doc.double_side(co) == 0                        # carbonyl on a branch point: centred
    doc = Doc()
    doc.add_ring(doc.ring_points(6, (0, 0), -math.pi / 2), aromatic=True)
    for i, b in doc.bonds.items():
        if b.order == 2:
            pa, pb = doc.atoms[b.a], doc.atoms[b.b]
            nx, ny = -(pb.y - pa.y), pb.x - pa.x
            mx, my = (pa.x + pb.x) / 2, (pa.y + pb.y) / 2
            assert doc.double_side(i) == (1 if (0 - mx) * nx + (0 - my) * ny > 0 else -1)   # toward the centre


def test_merge_and_remove():
    doc = Doc()
    a, b, c = doc.add_atom(0, 0), doc.add_atom(1, 0), doc.add_atom(2, 0)
    doc.add_bond(a, b); bc = doc.add_bond(b, c)
    doc.remove_bond(bc)
    assert c not in doc.atoms                              # a lone carbon goes with its last bond
    d = doc.add_atom(5, 5)
    doc.add_bond(b, d)
    e = doc.add_atom(0, 1); doc.add_bond(d, e)
    doc.merge_atom(e, a)
    assert doc.bond_between(d, a) is not None and e not in doc.atoms


def test_file_round_trip(tmp_path):
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('C[C@@H](N)C(=O)O'))
    doc.add_arrow(0, 3, 2, 3, 'equilibrium')
    doc.add_text(3, 3, 'alanine\n(2)')
    doc.style = {'bond_length': 20.0}
    p = tmp_path / 'x.sed'
    p.write_text(json.dumps(doc.to_dict()))
    back = Doc.from_dict(json.loads(p.read_text()))
    assert back.to_dict() == doc.to_dict()
    assert chem.analyse(back).smiles == 'C[C@@H](N)C(=O)O'
    i = back.add_atom(9, 9)
    assert i not in doc.atoms or i > max(doc.atoms)        # new ids do not collide with loaded ones


def test_partial_dict_keeps_inner_bonds():
    doc = Doc()
    doc.insert_rdkit(chem.mol_from_text('CCO'))
    atoms = sorted(doc.atoms)[:2]
    d = doc.to_dict([('atom', a) for a in atoms])
    assert len(d['atoms']) == 2 and len(d['bonds']) == 1
