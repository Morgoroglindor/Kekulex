// chem.js - formula, masses, adducts, SMILES and InChIKey, as the desktop's chem.py, and getting structures into the drawing.
import { ADDUCTS, ELECTRON, ELEMENTS, GROUPS, NICKNAMES } from './data.js';
import { parseBlock } from './doc.js';
import { caseGuess, condensedSmiles } from './condensed.js';
import { layoutSmiles } from './layout.js';

const mono = el => ELEMENTS[el][0];
const avg = el => ELEMENTS[el][1];

export function hill(counts) {
  const els = Object.keys(counts).sort();
  if (!counts.C) return els;
  return ['C', ...(counts.H ? ['H'] : []), ...els.filter(e => e !== 'C' && e !== 'H')];
}
export const formulaText = counts => hill(counts).filter(e => counts[e] > 0).map(e => e + (counts[e] > 1 ? counts[e] : '')).join('');
export const formulaHtml = (counts, charge = 0, radicals = 0) => hill(counts).filter(e => counts[e] > 0)
  .map(e => e + (counts[e] > 1 ? `<sub>${counts[e]}</sub>` : '')).join('') +
  (charge || radicals ? `<sup>${'•'.repeat(radicals)}${chargeText(charge)}</sup>` : '');
export const chargeText = z => (!z ? '' : (Math.abs(z) === 1 ? '' : Math.abs(z)) + (z > 0 ? '+' : '−'));

function monoOf(f) { return Object.entries(f).reduce((s, [e, n]) => s + mono(e) * n, 0); }

/** The SMILES a nickname or a condensed formula (SO3K, N(CH3)2 ...) stands for, its first atom attaching; or null. */
export const nicknameSmiles = label => NICKNAMES[label] || condensedSmiles(label);

// CH2, CH3•, C••, O•+, N•: a symbol, typed H, a sign, radical dots (. or * too), a sign after the dots
const ELEMENT_LABEL = /^([A-Za-z][a-z]?)(H(\d*))?\s*(\d*)\s*([+-]*)\s*([•·.*]*)\s*(?:(\d*)\s*([+-]*))?$/;
const isElement = s => !!ELEMENTS[s];

/** What a typed atom label means: ['element', symbol, charge], ['nickname', label] or null (model.parse_label):
 * an element (OH, NH2, Na+, CH3•), a nickname, a condensed formula (SO3K, CH2CH3), also typed in lower case. */
export function parseLabel(raw) {
  const text = raw.trim().replace('−', '-');
  if (NICKNAMES[text]) return ['nickname', text];
  const fixed = Object.keys(NICKNAMES).find(k => k.toLowerCase() === text.toLowerCase());
  if (fixed && text !== text.toUpperCase()) return ['nickname', fixed];   // boc is Boc, not B-O-C
  const guess = caseGuess(text, isElement);                // nh2 is NH2, not nihonium
  if (guess && guess !== text) { const r = parseLabel(guess); if (r) return r; }
  if (condensedSmiles(text)) return ['nickname', text];
  const m = ELEMENT_LABEL.exec(text);
  if (m) {
    const sym = m[1][0].toUpperCase() + m[1].slice(1).toLowerCase();
    if (ELEMENTS[sym]) {
      let signs = m[5], n = m[4] ? parseInt(m[4], 10) : signs.length;
      if (!signs && m[8]) { signs = m[8]; n = m[7] ? parseInt(m[7], 10) : signs.length; }
      return ['element', sym, !signs ? 0 : signs[0] === '+' ? n : -n];
    }
  }
  return fixed ? ['nickname', fixed] : null;
}

/** [typed H count or null, unpaired electrons] of an element label: 'CH2' -> [2, 0], 'CH3•' -> [3, 1], 'C••' -> [null, 2]. */
export function labelExtras(raw) {
  let text = raw.trim().replace('−', '-');
  const guess = caseGuess(text, isElement);
  if (guess) text = guess;
  const m = ELEMENT_LABEL.exec(text);
  if (!m || !ELEMENTS[m[1][0].toUpperCase() + m[1].slice(1).toLowerCase()]) return [null, 0];
  return [m[2] === undefined ? null : parseInt(m[3] || '1', 10), Math.min(m[6].length, 2)];
}

/** A group (SMILES) on atom a, its atoms laid out away from a's bonds (Doc.put_group). */
export function putGroup(RD, doc, a, smiles) {
  const m = RD.get_mol(smiles);
  if (!m) throw new Error(`Bad group ${smiles}`);
  m.set_new_coords(false);
  const g = parseBlock(m.get_molblock());
  m.delete();
  const at = doc.atoms.get(a);
  at.el = g.atoms[0].el; at.charge = g.atoms[0].charge; at.label = '';
  if (g.atoms.length === 1) return [a];
  const lens = g.bonds.map(b => Math.hypot(g.atoms[b.a].x - g.atoms[b.b].x, g.atoms[b.a].y - g.atoms[b.b].y)).sort((p, q) => p - q);
  const L = lens[Math.floor(lens.length / 2)] || 1;
  const pts = g.atoms.map(p => [p.x - g.atoms[0].x, p.y - g.atoms[0].y]);
  const cx = pts.slice(1).reduce((s, p) => s + p[0], 0) / (pts.length - 1), cy = pts.slice(1).reduce((s, p) => s + p[1], 0) / (pts.length - 1);
  const [fx, fy] = doc.neighbors(a).length ? doc.freeDirection(a) : [1, 0];
  const rot = Math.atan2(fy, fx) - Math.atan2(cy, cx), c = Math.cos(rot), s = Math.sin(rot);
  const ids = [a];
  for (let k = 1; k < pts.length; k++) {
    const x = pts[k][0] / L, y = pts[k][1] / L;
    ids.push(doc.addAtom(at.x + x * c - y * s, at.y + x * s + y * c, g.atoms[k].el, g.atoms[k].charge));
  }
  for (const b of g.bonds) doc.addBond(ids[b.a], ids[b.b], b.order);
  return ids;
}

export function setNickname(RD, doc, a, label) {
  const m = RD.get_mol(nicknameSmiles(label));
  const first = parseBlock(m.get_molblock()).atoms[0];
  m.delete();
  const at = doc.atoms.get(a);
  at.label = label; at.el = first.el; at.charge = 0; at.radical = 0; at.hcount = null;
}

/** A copy with every nickname written out: [doc, Map(atom id -> ids it became)]. */
export function expanded(RD, doc) {
  const d = doc.copy(), grown = new Map();
  for (const [a, at] of [...d.atoms]) {
    if (!at.label) continue;
    const smiles = nicknameSmiles(at.label);
    if (!smiles) { at.label = ''; continue; }
    const before = new Set(d.atoms.keys());
    putGroup(RD, d, a, smiles);
    grown.set(a, [...d.atoms.keys()].filter(i => !before.has(i)));
  }
  return [d, grown];
}

function masses(info, counts, charge, m, w) {
  Object.assign(info, { counts, charge, exact: m, average: w, formula: formulaText(counts), formulaHtml: formulaHtml(counts, charge, info.radicals || 0) });
  if (w) info.analysis = hill(counts).map(e => `${e}, ${(100 * avg(e) * counts[e] / w).toFixed(2)}`).join('; ');
  const odd = n => (n * (info.radicals || 0)) % 2 === 1;   // an odd-electron ion: [M]•+
  const label = (base, z, n = 1) => `[${base}]${odd(n) ? '•' : ''}${chargeText(z)}`;
  if (charge) info.ions.push({ label: label('M', charge), z: charge, mz: (m - charge * ELECTRON) / Math.abs(charge) });
  for (const [name, n, gain, loss, dz] of ADDUCTS) {
    const z = n * charge + dz;
    if (z === 0 || (charge && n > 1)) continue;
    if (charge && (z > 0) !== (charge > 0)) continue;
    const ion = {};
    for (const [e, c] of Object.entries(counts)) ion[e] = n * c;
    for (const [e, c] of Object.entries(gain)) ion[e] = (ion[e] || 0) + c;
    for (const [e, c] of Object.entries(loss)) ion[e] = (ion[e] || 0) - c;
    if (Object.values(ion).some(c => c < 0)) continue;
    info.ions.push({ label: label(name, z, n), z, mz: (n * m + monoOf(gain) - monoOf(loss) - z * ELECTRON) / Math.abs(z) });
  }
}

/** Formula, masses, ions and identifiers of the atoms in ids (all when null), as chem.analyse. */
export function analyse(RD, doc, ids = null) {
  const info = { atoms: 0, formula: '', formulaHtml: '', counts: {}, charge: 0, exact: 0, average: 0, analysis: '', ions: [], smiles: '', inchikey: '', errors: [], fragments: 0, radicals: 0 };
  let list = [...doc.atoms.keys()].filter(a => !ids || ids.has(a));
  info.atoms = list.length;
  if (!list.length) return info;
  info.fragments = doc.fragments(list).length;
  for (const a of list) { const p = doc.valenceProblem(a); if (p) info.errors.push(p); }
  const [full, grown] = expanded(RD, doc);
  const keep = new Set(list);
  for (const a of list) for (const g of grown.get(a) || []) keep.add(g);
  list = [...full.atoms.keys()].filter(a => keep.has(a));
  const counts = {};
  let m = 0, w = 0, charge = 0;
  const unknown = [];
  for (const a of list) {
    const at = full.atoms.get(a);
    if (!ELEMENTS[at.el]) { unknown.push(at.el); continue; }
    const hs = full.implicitH(a);
    counts[at.el] = (counts[at.el] || 0) + 1;
    if (hs) counts.H = (counts.H || 0) + hs;
    m += mono(at.el) + hs * mono('H');
    w += avg(at.el) + hs * avg('H');
    charge += at.charge;
    info.radicals += at.radical || 0;
  }
  if (unknown.length) info.errors.push(`Unknown element: ${[...new Set(unknown)].sort().join(', ')}`);
  masses(info, counts, charge, m, w);
  const [block] = full.toMolblock(keep);
  const fixedH = [...keep].some(a => full.atoms.has(a) && !full.atoms.get(a).label &&
    ((full.atoms.get(a).radical || 0) || (full.atoms.get(a).hcount !== null && full.atoms.get(a).hcount !== undefined)));
  const mol = RD.get_mol(fixedH ? withHydrogens(full, keep) : block);
  if (mol) {
    mol.remove_hs_in_place();                              // the drawn ring-fusion H do not belong in the SMILES
    info.smiles = mol.get_smiles();
    try { const inchi = mol.get_inchi(); if (inchi) info.inchikey = RD.get_inchikey_for_inchi(inchi); } catch { /* no InChI for this one */ }
    mol.delete();
  } else if (!info.errors.length) info.errors.push('The structure could not be read.');
  return info;
}

/** Adds a laid-out structure centred on (cx, cy), its bond length 1; returns the new atom ids. */
export function insertLayout(doc, layout, cx = 0, cy = 0) {
  const n = layout.atoms.length || 1;
  const mx = layout.atoms.reduce((s, a) => s + a.x, 0) / n, my = layout.atoms.reduce((s, a) => s + a.y, 0) / n;
  const lens = layout.bonds.map(b => Math.hypot(layout.atoms[b.a].x - layout.atoms[b.b].x, layout.atoms[b.a].y - layout.atoms[b.b].y)).sort((p, q) => p - q);
  const L = lens[Math.floor(lens.length / 2)] || 1;
  const ids = layout.atoms.map(a => doc.addAtom(cx + (a.x - mx) / L, cy + (a.y - my) / L, a.el, a.charge));
  layout.atoms.forEach((a, k) => { if (a.radical) doc.atoms.get(ids[k]).radical = a.radical; });
  for (const b of layout.bonds) doc.addBond(ids[b.a], ids[b.b], b.order, b.stereo || '');
  return ids;
}

export function insertSmiles(RD, doc, smiles, cx = 0, cy = 0) {
  const layout = layoutSmiles(RD, smiles);
  if (!layout) return null;
  return { ids: insertLayout(doc, layout, cx, cy), stereoOk: layout.stereoOk };
}

/** A MOL block's structure, keeping its coordinates and wedges. */
export function insertMolblock(doc, block, cx = 0, cy = 0) {
  return insertLayout(doc, parseBlock(block), cx, cy);
}

/** A SMILES with every atom's H written (for radicals and typed H counts, which a MOL block cannot say). */
function withHydrogens(doc, keep) {
  const order = [...doc.atoms.keys()].filter(a => keep.has(a));
  const idx = new Map(order.map((a, k) => [a, k]));
  const atoms = order.map(a => {
    const at = doc.atoms.get(a), h = doc.implicitH(a);
    const q = at.charge ? (at.charge > 0 ? '+' : '-') + (Math.abs(at.charge) > 1 ? Math.abs(at.charge) : '') : '';
    return `[${at.el}${h ? 'H' + (h > 1 ? h : '') : ''}${q}]`;
  });
  const bonds = [...doc.bonds.values()].filter(b => idx.has(b.a) && idx.has(b.b));
  // a SMILES by depth-first walk with ring closures
  const adj = order.map(() => []);
  bonds.forEach((b, k) => { adj[idx.get(b.a)].push([idx.get(b.b), b.order, k]); adj[idx.get(b.b)].push([idx.get(b.a), b.order, k]); });
  const seen = new Set(), tree = new Set(), closures = new Map();
  let ring = 1;
  const walk = (u, from) => {
    seen.add(u);
    for (const [v, , k] of adj[u]) {
      if (k === from) continue;
      if (seen.has(v)) { if (!tree.has(k) && !closures.has(k)) closures.set(k, ring++); continue; }
      tree.add(k); walk(v, k);
    }
  };
  const sym = o => (o === 2 ? '=' : o === 3 ? '#' : '');
  const done = new Set();
  const emit = u => {
    done.add(u);
    let s = atoms[u];
    for (const [, o, k] of adj[u]) if (closures.has(k)) { const r = closures.get(k); s += sym(o) + (r > 9 ? '%' + r : r); }
    const kids = adj[u].filter(([v, , k]) => tree.has(k) && !done.has(v));
    kids.forEach(([v, o], i) => { if (done.has(v)) return; const sub = sym(o) + emit(v); s += i < kids.length - 1 ? `(${sub})` : sub; });
    return s;
  };
  const parts = [];
  for (let u = 0; u < order.length; u++) if (!seen.has(u)) { walk(u, -1); parts.push(emit(u)); }
  return parts.join('.');
}

const ANALYSIS_LINES = {
  formula: i => `Chemical Formula: ${i.formula}${chargeText(i.charge)}`,
  exact: i => `Exact Mass: ${i.exact.toFixed(2)}`,
  mw: i => `Molecular Weight: ${i.average.toFixed(2)}`,
  analysis: i => `Elemental Analysis: ${i.analysis}`,
};

/** The text of a desktop "live" caption (Text.live) for these atoms, or null to leave it (keys the phone cannot work
 * out - names, properties, isotope patterns - keep their last text). */
export function liveCaption(RD, doc, ids, keys) {
  if (!keys.every(k => ANALYSIS_LINES[k] || k === 'ms')) return null;
  const info = analyse(RD, doc, ids);
  if (!info.atoms) return null;
  return keys.map(k => {
    if (k !== 'ms') return ANALYSIS_LINES[k](info);
    const dots = '•'.repeat(info.radicals % 3);
    if (!info.charge) return `${info.formula}${dots}\n${info.exact.toFixed(4)} (neutral)`;
    return `${info.formula}${dots}${chargeText(info.charge)}\nm/z ${((info.exact - info.charge * ELECTRON) / Math.abs(info.charge)).toFixed(4)}`;
  }).join('\n');
}

/** Rewrites the live captions of the structures they are grouped with (after an edit on the phone). */
export function refreshLiveCaptions(RD, doc) {
  for (const [i, t] of doc.texts) {
    if (!t.live || !t.live.length) continue;
    const g = doc.groupOf('text', i);
    if (!g) continue;
    const ids = new Set(g.members.filter(([k, j]) => k === 'atom' && doc.atoms.has(j)).map(([, j]) => j));
    if (!ids.size) continue;
    const text = liveCaption(RD, doc, ids, t.live);
    if (text !== null && text !== t.text) t.text = text;
  }
}

export { GROUPS };
