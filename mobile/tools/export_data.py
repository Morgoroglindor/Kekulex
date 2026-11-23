"""Writes www/core/data.js and tests/reference.json from the desktop program, so that both share one set of data and the
phone's results can be compared with the desktop's. Run with the desktop's Python after changing model.py, chem.py,
templates.py or layout2d.py:   ..\\.venv\\Scripts\\python tools\\export_data.py"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MOBILE = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(MOBILE))
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6 import QtWidgets
app = QtWidgets.QApplication([])
from rdkit import Chem

import chem
import layout2d
import model
import templates
from model import Doc

PT = Chem.GetPeriodicTable()
elements = {PT.GetElementSymbol(z): [PT.GetMostCommonIsotopeMass(z), PT.GetAtomicWeight(z)] for z in range(1, 119)}
data = dict(GROUPS=model.GROUPS, NICKNAMES=model.NICKNAMES, VALENCE=model.VALENCE, VALENCE_ELECTRONS=model.VALENCE_ELECTRONS,
            ELEMENTS=elements, ELECTRON=chem.ELECTRON, ADDUCTS=chem.ADDUCTS_POS + chem.ADDUCTS_NEG,
            TEMPLATES=dict(templates.TEMPLATES),
            LAYOUT=[[n, s, c] for n, s, c in layout2d.TEMPLATES], SAME_SHAPE=layout2d.SAME_SHAPE)
import canvas  # noqa: E402  (the desktop's ChemDraw-style hotkeys: the web version uses the same tables)
import render  # noqa: E402
data['KEYS'] = dict(ATOM=canvas.ATOM_KEYS, SPROUT={k: list(v) for k, v in canvas.SPROUT_KEYS.items()},
                    FUSE={k: list(v) for k, v in canvas.FUSE_KEYS.items()}, BOND_STYLE=canvas.BOND_STYLE_KEYS,
                    TOOL=canvas.TOOL_KEYS)
data['STYLE_SHEETS'] = {name: {k: v for k, v in s.items() if isinstance(v, (int, float, str, bool))}
                        for name, s in render.STYLE_SHEETS.items()}
data['DEFAULT_STYLE'] = {k: v for k, v in render.DEFAULT_STYLE.items() if isinstance(v, (int, float, str, bool))}
with open(os.path.join(MOBILE, 'www', 'core', 'data.js'), 'w', encoding='utf-8') as f:
    f.write('// Made by tools/export_data.py from the desktop program: do not edit.\n')
    for k, v in data.items():
        f.write(f'export const {k} = {json.dumps(v, ensure_ascii=False)};\n')

SMILES = ['CCO', 'CC(=O)O', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
          'C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'C/C=C/C(=O)O', 'C/C=C\\C(=O)O', 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
          'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C',
          'CC(C)CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
          'CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C',
          'C[C@]12CC[C@@H](C([C@@H]1CC[C@@]3([C@@H]2CC=C4[C@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)C(=O)O)C)C)(C)C)O',
          'O=C([O-])[C@@H]1Cc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c(O)cc2[N+]1=CC=C1C[C@@H](C(=O)O)NC(C(=O)O)=C1',
          'CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5',
          'COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O',
          'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C', 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N',
          'FC(F)(F)c1ccc(Br)cc1Cl', 'CS(C)(=O)=O', 'C[NH3+]', 'CC(=O)[O-].[Na+]', 'O=[N+]([O-])c1ccccc1',
          'CCCCCCCCCCCCCCCC(=O)O', 'C#N', '[Cl-]', 'O', 'P(=O)(O)(O)O']
rows = []
for s in SMILES:
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(s))
    info = chem.analyse(d)
    rows.append(dict(smiles=s, canonical=Chem.MolToSmiles(Chem.MolFromSmiles(s)), formula=info.formula, charge=info.charge,
                     exact=info.exact, average=info.average, analysis=info.analysis, inchikey=info.inchikey,
                     ions=[[i.label, i.z, i.mz] for i in info.ions], sed=d.to_dict()))
PUB = 'C1[C@H](N=C(C=C1/C=C/N2[C@@H](CC3=CC(=C(C=C32)O)O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)C(=O)O)C(=O)O)C(=O)O'
usual = [dict(smiles=s, usual=layout2d.usual_form(s)) for s in (PUB, 'CCO', 'c1ccncc1')]
import condensed  # noqa: E402
CONDENSED = ['SO3K', 'SO3Na', 'CO2Na', 'CH2CH3', 'OC2H5', 'COCH3', 'C(O)CH3', 'OCOCH3', 'NHCOCH3', 'CONH2', 'CH2CO2H',
             'N(CH3)2', 'CH(CH3)2', '(CH2)3CH3', 'B(OH)2', 'P(O)(OH)2', 'PO(OEt)2', 'SO2NH2', 'SO2Cl', 'SO3Me', 'PO3H2',
             'OCF3', 'CN', 'NCO', 'SCN', 'CH=CH2', 'C6H5', 'SO3-', 'CO2-', 'NMe3+', 'N(CH3)3+', 'CH2Ph', 'CONHMe',
             'CO2tBu', 'OTBS', 'CCl3', 'CH2OMe']
condensed_rows = [[t, Chem.MolToSmiles(Chem.MolFromSmiles(condensed.smiles(t)))] for t in CONDENSED]
LABELS = ['N', 'OH', 'NH2', 'Na+', 'O-', 'OMe', 'Boc', 'cl', 'Xx', 'ome', 'so3k', 'ch2cl', 'nh2', 'cf3', 'no2', 'boc',
          'CH2', 'CH3•', 'C••', 'O•+', 'ch2.', 'N•', 'CH3+', 'SO3K', 'N(CH3)2']
labels = [[t, list(model.parse_label(t) or []), list(model.label_extras(t))] for t in LABELS]
radical_rows = []
for s, tweak in (('[CH3]', None), ('[CH2]', None), ('C[CH2]', None), ('CC(C)=O', 'O+.'), ('CC#[O+]', None)):
    d = Doc()
    d.insert_rdkit(Chem.MolFromSmiles(s))
    if tweak:
        o = next(a for a, at in d.atoms.items() if at.el == 'O')
        d.atoms[o].charge, d.atoms[o].radical = 1, 1
    info = chem.analyse(d)
    radical_rows.append(dict(sed=d.to_dict(), formula=info.formula, ions=[[i.label, i.mz] for i in info.ions[:3]]))
with open(os.path.join(MOBILE, 'tests', 'reference.json'), 'w', encoding='utf-8') as f:
    json.dump(dict(rows=rows, usual=usual, condensed=condensed_rows, labels=labels, radicals=radical_rows,
                   parse=[[t, list(model.parse_label(t) or [])] for t in ('N', 'OH', 'NH2', 'Na+', 'O-', 'OMe', 'Boc', 'cl', 'Xx', 'ome')]), f)
print(len(rows), 'reference molecules')
