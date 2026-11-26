// The phone's chemistry against what the desktop computes (tests/reference.json, from tools/export_data.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { test } from 'node:test';

import { TEMPLATES } from '../www/core/data.js';
import { analyse, insertMolblock, insertSmiles, labelExtras, liveCaption, parseLabel, refreshLiveCaptions, setNickname } from '../www/core/chem.js';
import { condensedSmiles, setRD } from '../www/core/condensed.js';
import { Doc } from '../www/core/doc.js';
import { findRings, layoutSmiles, usualForm } from '../www/core/layout.js';

const require = createRequire(import.meta.url);
const RD = await require('@rdkit/rdkit')();
const REF = JSON.parse(readFileSync(new URL('./reference.json', import.meta.url), 'utf8'));
setRD(RD);
const strip = s => s.replace(/[/\\]/g, '');
const canon = s => { const m = RD.get_mol(s); const c = m.get_smiles(); m.delete(); return c; };

for (const row of REF.rows) {
  test(`analysis as on the desktop: ${row.smiles.slice(0, 30)}`, () => {
    const d = new Doc();
    assert.ok(insertSmiles(RD, d, row.smiles));
    const info = analyse(RD, d);
    assert.equal(info.formula, row.formula);
    assert.equal(info.charge, row.charge);
    assert.ok(Math.abs(info.exact - row.exact) < 1e-6, `${info.exact} ${row.exact}`);
    assert.ok(Math.abs(info.average - row.average) < 1e-3);
    assert.equal(info.analysis, row.analysis);
    assert.deepEqual(info.ions.map(i => [i.label, i.z]), row.ions.map(i => [i[0], i[1]]));
    info.ions.forEach((ion, k) => assert.ok(Math.abs(ion.mz - row.ions[k][2]) < 1e-6));
    assert.equal(info.inchikey, row.inchikey);
  });

  test(`the drawing is the same molecule: ${row.smiles.slice(0, 30)}`, () => {
    const d = new Doc();
    const r = insertSmiles(RD, d, row.smiles);
    assert.ok(r.stereoOk, 'stereo checked');
    const [block] = d.toMolblock();
    const m = RD.get_mol(block);
    m.remove_hs_in_place();
    const back = m.get_smiles();
    m.delete();
    if (/[/\\]/.test(row.canonical)) assert.equal(back, canon(row.smiles));
    else assert.equal(strip(back), strip(canon(row.smiles)));
  });

  test(`a desktop .sed reads the same: ${row.smiles.slice(0, 30)}`, () => {
    const d = Doc.fromDict(row.sed);
    const info = analyse(RD, d);
    assert.equal(info.formula, row.formula);
    assert.equal(info.inchikey, row.inchikey);
    assert.deepEqual(JSON.parse(JSON.stringify(d.toDict().atoms)), row.sed.atoms);
    assert.deepEqual(JSON.parse(JSON.stringify(d.toDict().bonds)), row.sed.bonds);
  });
}

const lengths = g => g.bonds.map(b => Math.hypot(g.atoms[b.a].x - g.atoms[b.b].x, g.atoms[b.a].y - g.atoms[b.b].y));
const row = name => REF.rows.find(r => r.smiles.startsWith(name));

test('usual drawing: lupeol with its template, even bonds, ring-fusion H and no wedges in rings', () => {
  const g = layoutSmiles(RD, row('CC(=C)[C@@H]1').smiles);
  assert.equal(g.template, 'lupane');
  assert.equal(g.clashes, 0);
  const L = lengths(g);
  assert.ok(Math.max(...L) / Math.min(...L) < 1.12);
  assert.equal(g.atoms.filter(a => a.el === 'H').length, 4);
  const inRing = new Set();
  for (const r of findRings(g)) r.forEach((a, k) => inRing.add(`${a},${r[(k + 1) % r.length]}`));
  for (const b of g.bonds) if (b.stereo) assert.ok(!inRing.has(`${b.a},${b.b}`) && !inRing.has(`${b.b},${b.a}`), 'wedge on a ring bond');
});

test('usual drawing: oleanane rises from the oxygenated ring A (y is down)', () => {
  const g = layoutSmiles(RD, row('C[C@]12CC[C@@H](C(').smiles);
  assert.equal(g.template, 'oleanane');
  const o = g.atoms.find(a => a.el === 'O' && g.bonds.some(b => (b.a === g.atoms.indexOf(a) || b.b === g.atoms.indexOf(a)) && b.order === 1));
  const cx = g.atoms.reduce((s, a) => s + a.x, 0) / g.atoms.length;
  assert.ok(g.atoms.some(a => a.el === 'O' && a.x < cx && a.y > 0), 'an OH at the lower left');
  assert.ok(o);
});

test('usual drawing: betanin has its ring NH at the bottom and its carboxylate right of N+', () => {
  const g = layoutSmiles(RD, row('O=C([O-])[C@@H]1').smiles);
  assert.equal(g.template, 'betalamic');
  assert.equal(g.clashes, 0);
  const L = lengths(g);
  assert.ok(Math.max(...L) / Math.min(...L) < 1.15, `${Math.max(...L)} ${Math.min(...L)}`);
  const plus = g.atoms.findIndex(a => a.el === 'N' && a.charge === 1), ringN = g.atoms.findIndex(a => a.el === 'N' && !a.charge);
  assert.ok(g.atoms[ringN].y > g.atoms[plus].y, 'ring NH below the iminium');
  const nb = g.bonds.flatMap(b => (b.a === plus ? [b.b] : b.b === plus ? [b.a] : []));
  const acidSide = nb.find(n => g.bonds.some(b => (b.a === n || b.b === n) && g.atoms[b.a === n ? b.b : b.a].el === 'C' && g.bonds.some(c => (c.a === (b.a === n ? b.b : b.a) || c.b === (b.a === n ? b.b : b.a)) && c.order === 2 && g.atoms[c.a === (b.a === n ? b.b : b.a) ? c.b : c.a].el === 'O')));
  assert.ok(g.atoms[acidSide].x > g.atoms[plus].x);
});

test('name lookup form of betanin as the desktop makes it', () => {
  for (const { smiles, usual } of REF.usual) assert.equal(strip(canon(usualForm(RD, smiles))), strip(canon(usual)));
});

test('labels as on the desktop', () => {
  for (const [text, want] of REF.parse) assert.deepEqual(parseLabel(text) || [], want, text);
});

test('nicknames count in full', () => {
  const d = new Doc();
  const a = d.addAtom(0, 0), b = d.addAtom(1, 0);
  d.addBond(a, b);
  setNickname(RD, d, b, 'OMe');
  assert.equal(analyse(RD, d).formula, 'C2H6O');
  assert.equal(d.atoms.size, 2, 'the drawing itself keeps the label');
});

test('drawing tools: benzene, valence error, merge, delete', () => {
  const d = new Doc();
  d.addRing(d.ringPoints(6, 0, 0, Math.PI / 2), true);
  assert.equal(d.atoms.size, 6);
  assert.equal([...d.bonds.values()].filter(b => b.order === 2).length, 3);
  assert.equal(analyse(RD, d).formula, 'C6H6');
  const first = [...d.bonds.keys()][0];
  d.addRing(d.ringOnBond(first, 6), true);
  assert.equal(d.atoms.size, 10);
  assert.equal(analyse(RD, d).formula, 'C10H8');
  const e = new Doc();
  const c = e.addAtom(0, 0);
  for (let k = 0; k < 5; k++) e.addBond(c, e.addAtom(Math.cos(k), Math.sin(k)));
  assert.ok(analyse(RD, e).errors[0].includes('takes at most 4'));
  e.removeAtom(c);
  assert.equal(e.bonds.size, 0);
});

test('every template of the desktop goes in', () => {
  for (const [cat, items] of Object.entries(TEMPLATES)) {
    for (const [name, smiles] of items) {
      const d = new Doc();
      const r = insertSmiles(RD, d, smiles);
      assert.ok(r, `${cat} ${name}`);
      assert.ok(r.stereoOk, `stereo of ${cat} ${name}`);
    }
  }
});

test('betalains stay out of the prepacked molecules', () => {
  const names = Object.values(TEMPLATES).flat().map(t => t[0].toLowerCase()).join(' ');
  for (const word of ['betal', 'betan', 'xanthin']) assert.ok(!names.includes(word));
});

test('a MOL file keeps its coordinates and wedges', () => {
  const d = new Doc();
  insertSmiles(RD, d, 'C[C@H](N)C(=O)O');
  const [block] = d.toMolblock();
  const e = new Doc();
  insertMolblock(e, block);
  assert.equal(analyse(RD, e).inchikey, analyse(RD, d).inchikey);
});

test('kept kinds survive a round trip', () => {
  const dict = { atoms: [{ id: 3, x: 0, y: 0, el: 'C', charge: 0, label: '', color: '', highlight: '', charge_pos: null }], shapes: [{ id: 9, x1: 0, y1: 0, x2: 1, y2: 1, kind: 'plasmid', params: { bp: 100 } }] };
  const d = Doc.fromDict(dict);
  assert.deepEqual(d.toDict().shapes, dict.shapes);
  assert.equal(d.otherCount(), 1);
  assert.ok(d.addAtom(1, 1) > 9, 'new ids stay clear of kept ones');
});

// ---------------------------------------------------------------------------------------------- condensed, radicals, groups
for (const [text, smiles] of REF.condensed) {
  test(`condensed label as on the desktop: ${text}`, () => {
    const mine = condensedSmiles(text);
    assert.ok(mine, `${text} not read`);
    assert.equal(canon(mine), smiles);
  });
}

for (const [text, parsed, extras] of REF.labels) {
  test(`typed label as on the desktop: ${text}`, () => {
    assert.deepEqual(parseLabel(text) || [], parsed);
    assert.deepEqual(labelExtras(text), extras);
  });
}

for (const row of REF.radicals) {
  test(`radicals and radical ions as on the desktop: ${row.formula}`, () => {
    const d = Doc.fromDict(row.sed);
    const info = analyse(RD, d);
    assert.equal(info.formula, row.formula);
    row.ions.forEach(([label, mz], k) => { assert.equal(info.ions[k].label, label); assert.ok(Math.abs(info.ions[k].mz - mz) < 1e-6); });
  });
}

test('a nickname written as a condensed formula counts in full', () => {
  const d = new Doc();
  insertSmiles(RD, d, 'c1ccccc1');
  const a = [...d.atoms.keys()][0];
  const end = d.addAtom(d.atoms.get(a).x + 1, d.atoms.get(a).y);
  d.addBond(a, end);
  setNickname(RD, d, end, 'SO3K');
  assert.equal(analyse(RD, d).formula, 'C6H5KO3S');
});

test('groups are kept and pick up their whole group; captions follow their structure', () => {
  const d = new Doc();
  const ids = insertSmiles(RD, d, 'CCO').ids;
  const t = d.addText(0, 3, 'x');
  d.texts.get(t).live = ['formula'];
  d.other.groups = [{ id: 99, members: [...ids.map(i => ['atom', i]), ['text', t]] }];
  const back = Doc.fromDict(JSON.parse(JSON.stringify(d.toDict())));
  assert.equal(back.other.groups.length, 1);
  assert.ok(back.withGroups(new Set([`atom:${ids[0]}`])).has(`text:${t}`));
  back.atoms.get(ids[2]).el = 'N';
  refreshLiveCaptions(RD, back);
  assert.equal(back.texts.get(t).text, 'Chemical Formula: C2H7N');
  assert.equal(liveCaption(RD, back, new Set(ids), ['name']), null);   // names stay as they were
});

test('tapping an element again cycles its hydrogens, as ChemDraw does', () => {
  const d = new Doc();
  const a = d.addAtom(0, 0, 'N'), b = d.addAtom(1, 0, 'C');
  d.addBond(a, b);
  const seen = [];
  for (let k = 0; k < 4; k++) { d.cycleHydrogens(a); seen.push([d.implicitH(a), Boolean(d.valenceProblem(a))]); }
  assert.deepEqual(seen, [[1, true], [0, true], [2, false], [1, true]]);
  const c = [];
  for (let k = 0; k < 5; k++) { d.cycleHydrogens(b); c.push(d.implicitH(b)); }
  assert.deepEqual(c, [3, 2, 1, 0, 3]);
});

test('a formula text: digits after a letter or bracket are subscripts, as the desktop render.runs', async () => {
  const { runs } = await import('../www/core/draw.js');
  const subs = (line) => runs(line, true).filter(([, s]) => s).map(([t]) => t);
  assert.deepEqual(subs('Pd(PPh3)4'), ['3', '4']);
  assert.deepEqual(subs('K2CO3, toluene/H2O'), ['2', '3', '2']);
  assert.deepEqual(subs('100 °C, 12 h, 1H'), []);
  assert.deepEqual(runs('H2O', false), [['H2O', false]]);
  assert.equal(runs('SnCl2·2H2O', true).map(([t]) => t).join(''), 'SnCl2·2H2O');
});

// ------------------------------------------------------------------------------------------ the desktop's keys (keys.js)
const smilesOf = d => { const [block] = d.toMolblock(); const m = RD.get_mol(block); const s = m.get_smiles(); m.delete(); return s; };

test('hotkeys over an atom, as the desktop: elements, nicknames, the same key again, bonds sprouted', async () => {
  const { hotkey } = await import('../www/core/keys.js');
  const d = new Doc();
  insertSmiles(RD, d, 'CC');
  const [a, b] = [...d.atoms.keys()];
  assert.equal(hotkey(RD, d, 'o', ['atom', a]), 'edited');
  assert.equal(canon(smilesOf(d)), canon('CO'));                     // o: oxygen
  hotkey(RD, d, 'o', ['atom', a]);                                    // again: one H fewer (O with no H: a radical-free O-)
  assert.equal(d.atoms.get(a).hcount, 0);
  assert.equal(hotkey(RD, d, 'C', ['atom', b]), 'edited');           // Shift+C: chlorine
  assert.equal(d.atoms.get(b).el, 'Cl');
  const e = new Doc(); insertSmiles(RD, e, 'CC');
  const c = [...e.atoms.keys()][0];
  hotkey(RD, e, '2', ['atom', c]);                                    // 2: C=O sprouted
  assert.equal(canon(smilesOf(e)), canon('CC=O'));
  hotkey(RD, e, 'O', ['atom', [...e.atoms.keys()][1]]);              // Shift+O: OMe, a nickname
  assert.equal(e.atoms.get([...e.atoms.keys()][1]).label, 'OMe');
  hotkey(RD, e, '+', ['atom', c]);
  assert.equal(e.atoms.get(c).charge, 1);
  assert.equal(hotkey(RD, e, 't', ['atom', c]), 'label');
});

test('hotkeys over a bond: order, wedge and its turn, a fused ring; tool keys elsewhere', async () => {
  const { hotkey } = await import('../www/core/keys.js');
  const d = new Doc();
  insertSmiles(RD, d, 'CCC');
  const bond = [...d.bonds.keys()][0];
  hotkey(RD, d, '2', ['bond', bond]);
  assert.equal(d.bonds.get(bond).order, 2);
  hotkey(RD, d, 'w', ['bond', bond]);
  const bd = d.bonds.get(bond), first = bd.a;
  assert.deepEqual([bd.order, bd.stereo], [1, 'wedge']);
  hotkey(RD, d, 'w', ['bond', bond]);                                 // again: the other way round
  assert.notEqual(d.bonds.get(bond).a, first);
  const atoms = d.atoms.size;
  hotkey(RD, d, '6', ['bond', bond]);                                 // a six-membered ring fused on it
  assert.equal(d.atoms.size, atoms + 4);
  assert.deepEqual(hotkey(RD, d, ' ', null), ['select']);
  assert.deepEqual(hotkey(RD, d, 'j', null), ['ring', 'benzene']);
  assert.equal(hotkey(RD, d, 'q', null), null);                       // nothing there
});

test('missing charges, as the desktop: N with four bonds is N+, O typed with no H on one bond is O−', async () => {
  const { checkStructure, fixCharges, neededCharge } = await import('../www/core/keys.js');
  const d = new Doc();
  insertSmiles(RD, d, 'C[N+](C)(C)C');
  const n = [...d.atoms.entries()].find(([, a]) => a.el === 'N')[0];
  d.atoms.get(n).charge = 0;                                          // drawn without its charge
  assert.equal(neededCharge(d, n), 1);
  assert.ok(checkStructure(d).some(p => p.includes('needs a charge')));
  assert.deepEqual(fixCharges(d), { [n]: 1 });
  assert.equal(neededCharge(d, [...d.atoms.keys()][0]), null);         // carbon: the chemist's call
});

test('align and arrange as a reaction scheme', async () => {
  const { align, arrangeScheme, objects } = await import('../www/core/keys.js');
  const d = new Doc();
  insertSmiles(RD, d, 'CCO', 0, 0);
  insertSmiles(RD, d, 'CC(=O)O', 10, 3);
  d.addArrow(4, -2, 6, -2, 'forward');
  assert.ok(align(d, null, 'vcenter'));
  const mids = objects(d).map(o => (o.box[1] + o.box[3]) / 2);
  assert.ok(Math.max(...mids) - Math.min(...mids) < 1e-9);
  assert.ok(arrangeScheme(d, null));
  const objs = objects(d).sort((p, q) => p.box[0] - q.box[0]);
  assert.deepEqual(objs.map(o => o.kind), ['mol', 'arrow', 'mol']);
  for (let k = 1; k < objs.length; k++) assert.ok(objs[k].box[0] > objs[k - 1].box[2]);   // one row, in order, apart
  assert.ok(Math.abs(d.arrows.values().next().value.x2 - d.arrows.values().next().value.x1) >= 3);
});

// ------------------------------------------------------------------------------------------- MS fragments
import { predict as predictFragments } from '../www/core/fragments.js';

for (const row of REF.fragments) {
  test(`MS fragments as on the desktop: ${row.smiles.slice(0, 28)} ${row.mode}`, () => {
    const d = Doc.fromDict(row.sed);
    const { pre, frags } = predictFragments(RD, d, null, row.mode, 2);
    assert.ok(pre, 'a precursor');
    assert.equal(pre.ion, row.pre[0]);
    assert.equal(pre.label, row.pre[2]);
    assert.ok(Math.abs(pre.mz - row.pre[1]) < 1e-6);
    const got = frags.map(f => [f.ion, f.lost, f.how, f.step]);
    assert.deepEqual(got, row.frags.map(f => [f[1], f[2], f[3], f[5]]));
    frags.forEach((f, k) => {
      assert.ok(Math.abs(f.mz - row.frags[k][0]) < 1e-6, `${f.ion} m/z`);
      assert.ok(Math.abs(f.score - row.frags[k][4]) < 1e-9, `${f.ion} ${f.how} score ${f.score} ${row.frags[k][4]}`);
      assert.deepEqual([...f.atoms].sort((p, q) => p - q), row.frags[k][6]);
      assert.ok(f.mech, `${f.how}: a mechanism`);
    });
  });
}
