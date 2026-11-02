// chem.js - formula, masses, adducts, SMILES and InChIKey, as the desktop's chem.py, and getting structures into the drawing.
import { ADDUCTS, ELECTRON, ELEMENTS, GROUPS, NICKNAMES } from './data.js';
import { parseBlock } from './doc.js';
import { layoutSmiles } from './layout.js';

const mono = el => ELEMENTS[el][0];
const avg = el => ELEMENTS[el][1];

export function hill(counts) {
  const els = Object.keys(counts).sort();
  if (!counts.C) return els;
  return ['C', ...(counts.H ? ['H'] : []), ...els.filter(e => e !== 'C' && e !== 'H')];
}
export const formulaText = counts => hill(counts).filter(e => counts[e] > 0).map(e => e + (counts[e] > 1 ? counts[e] : '')).join('');
export const formulaHtml = (counts, charge = 0) => hill(counts).filter(e => counts[e] > 0)
  .map(e => e + (counts[e] > 1 ? `<sub>${counts[e]}</sub>` : '')).join('') + (charge ? `<sup>${chargeText(charge)}</sup>` : '');
export const chargeText = z => (!z ? '' : (Math.abs(z) === 1 ? '' : Math.abs(z)) + (z > 0 ? '+' : '−'));

function monoOf(f) { return Object.entries(f).reduce((s, [e, n]) => s + mono(e) * n, 0); }

/** What a typed atom label means: ['element', symbol, charge], ['nickname', label] or null (model.parse_label). */
export function parseLabel(raw) {
  const text = raw.trim().replace('−', '-');
  if (NICKNAMES[text]) return ['nickname', text];
  const m = /^([A-Za-z][a-z]?)(?:H\d*)?\s*(\d*)\s*([+-]*)$/.exec(text);
  if (m) {
    const sym = m[1][0].toUpperCase() + m[1].slice(1).toLowerCase();
    if (ELEMENTS[sym]) {
      const n = m[2] ? parseInt(m[2], 10) : m[3].length;
      return ['element', sym, !m[3] ? 0 : m[3][0] === '+' ? n : -n];
    }
  }
  const low = Object.keys(NICKNAMES).find(k => k.toLowerCase() === text.toLowerCase());
  return low ? ['nickname', low] : null;
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
  const m = RD.get_mol(NICKNAMES[label]);
  const first = parseBlock(m.get_molblock()).atoms[0];
  m.delete();
  const at = doc.atoms.get(a);
  at.label = label; at.el = first.el; at.charge = 0;
}

/** A copy with every nickname written out: [doc, Map(atom id -> ids it became)]. */
export function expanded(RD, doc) {
  const d = doc.copy(), grown = new Map();
  for (const [a, at] of [...d.atoms]) {
    if (!at.label) continue;
    const smiles = NICKNAMES[at.label];
    if (!smiles) { at.label = ''; continue; }
    const before = new Set(d.atoms.keys());
    putGroup(RD, d, a, smiles);
    grown.set(a, [...d.atoms.keys()].filter(i => !before.has(i)));
  }
  return [d, grown];
}

function masses(info, counts, charge, m, w) {
  Object.assign(info, { counts, charge, exact: m, average: w, formula: formulaText(counts), formulaHtml: formulaHtml(counts, charge) });
  if (w) info.analysis = hill(counts).map(e => `${e}, ${(100 * avg(e) * counts[e] / w).toFixed(2)}`).join('; ');
  const label = (base, z) => `[${base}]${chargeText(z)}`;
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
    info.ions.push({ label: label(name, z), z, mz: (n * m + monoOf(gain) - monoOf(loss) - z * ELECTRON) / Math.abs(z) });
  }
}

/** Formula, masses, ions and identifiers of the atoms in ids (all when null), as chem.analyse. */
export function analyse(RD, doc, ids = null) {
  const info = { atoms: 0, formula: '', formulaHtml: '', counts: {}, charge: 0, exact: 0, average: 0, analysis: '', ions: [], smiles: '', inchikey: '', errors: [], fragments: 0 };
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
  }
  if (unknown.length) info.errors.push(`Unknown element: ${[...new Set(unknown)].sort().join(', ')}`);
  masses(info, counts, charge, m, w);
  const [block] = full.toMolblock(keep);
  const mol = RD.get_mol(block);
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

export { GROUPS };
