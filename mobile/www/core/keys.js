// keys.js - the desktop's ChemDraw-style hotkeys and structure commands, for a keyboard and mouse (the web version on a
// PC): the same tables as canvas.py (KEYS, exported by tools/export_data.py) and the same rules as canvas.hotkey,
// model.needed_charge / fix_charges and canvas.align. Kept apart from app.js so the tests can run them.
import { putGroup, setNickname } from './chem.js';
import { GROUPS, KEYS, NICKNAMES, VALENCE, VALENCE_ELECTRONS } from './data.js';

/** A new bond from atom a, in its next free direction, to a group (SMILES; 'C' a carbon) by order and stereo: returns the
 * new atom (model.Doc.sprout). */
export function sprout(RD, doc, a, smiles = 'C', order = 1, stereo = '') {
  const at = doc.atoms.get(a), ang = doc.nextBondAngle(a);
  const b = doc.atomAtOrNew(at.x + Math.cos(ang), at.y + Math.sin(ang), 'C', [a]);
  const id = doc.addBond(a, b, order, stereo);
  if (id !== null && stereo) {
    const bd = doc.bonds.get(id);
    if (bd.a !== a) { const t = bd.a; bd.a = bd.b; bd.b = t; }
  }
  if (smiles !== 'C' && doc.neighbors(b).length <= 1) putGroup(RD, doc, b, smiles);
  return b;
}

/** A hotkey ch over the hovered item (['atom', id] / ['bond', id] / null). Returns what it did: 'edited' (the
 * drawing changed: the caller took a snapshot first), 'label' (type the atom's label), a tool name (no hover: a
 * tool key), or null (not a hotkey here). canvas.Canvas.hotkey, key for key. */
export function hotkey(RD, doc, ch, hover) {
  if (!ch) return null;
  if (hover && hover[0] === 'atom') {
    const a = hover[1], at = doc.atoms.get(a);
    if (!at) return null;
    if (ch === 't') return 'label';
    if (ch === '.') { at.radical = ((at.radical || 0) + 1) % 3; return 'edited'; }
    if (['+', '=', '-', '_'].includes(ch)) { at.charge += ['+', '='].includes(ch) ? 1 : -1; return 'edited'; }
    if (ch in KEYS.ATOM) {
      const name = KEYS.ATOM[ch];
      if (NICKNAMES[name]) { setNickname(RD, doc, a, name); at.hcount = null; at.radical = 0; }
      else if (at.el === name && !at.label) doc.cycleHydrogens(a);           // the same key again: one H fewer
      else { at.el = name; at.label = ''; at.hcount = null; at.iso = 0; }      // another element: all its H
      return 'edited';
    }
    if (ch in KEYS.SPROUT) {
      const [group, order, stereo] = KEYS.SPROUT[ch];
      sprout(RD, doc, a, GROUPS[group] || group, order, stereo);
      return 'edited';
    }
    if (ch === '9') { sprout(RD, doc, a); sprout(RD, doc, a); return 'edited'; }
    return null;
  }
  if (hover && hover[0] === 'bond') {
    const bd = doc.bonds.get(hover[1]);
    if (!bd) return null;
    if (['1', '2', '3'].includes(ch)) { Object.assign(bd, { order: Number(ch), stereo: '', style: '', side: null }); return 'edited'; }
    if (['l', 'r', 'c'].includes(ch) && bd.order === 2) { bd.side = { l: -1, r: 1, c: 0 }[ch]; return 'edited'; }
    if (ch in KEYS.FUSE) {
      const [n, aromatic] = KEYS.FUSE[ch];
      doc.addRing(doc.ringOnBond(hover[1], n), aromatic);
      return 'edited';
    }
    if (['w', 'h', 'W'].includes(ch)) {
      const stereo = ch === 'w' ? 'wedge' : 'hash';
      if (bd.stereo === stereo) { const t = bd.a; bd.a = bd.b; bd.b = t; }
      Object.assign(bd, { stereo, order: 1, style: '' });
      return 'edited';
    }
    if (ch in KEYS.BOND_STYLE) {
      const style = KEYS.BOND_STYLE[ch];
      if (style === 'dashed' && bd.order === 2) bd.style = 'dashed';
      else Object.assign(bd, { style, order: 1, stereo: '' });
      return 'edited';
    }
    return null;
  }
  if (!hover && ch in KEYS.TOOL) return WEB_TOOL[KEYS.TOOL[ch]] || null;
  return null;
}

/** The desktop's tool names as the web version's tools ([tool, its choice]); tools it does not have are left out. */
export const WEB_TOOL = { select: ['select'], bond1: ['bond', 'single'], chain: ['chain'], arrow: ['arrow'], text: ['text'],
  benzene: ['ring', 'benzene'], cyclopentadiene: ['ring', 'r5'] };

// ------------------------------------------------------------------------------------------ charges (model.py)
/** The formal charge an atom's bonds and H call for when it has none (N or P with four bonds: +1; O with one bond and
 * no H: -1 ...), or null: it is fine, or carbon (a cation or an anion is the chemist's call), or no charge fits. */
export function neededCharge(doc, a) {
  const at = doc.atoms.get(a);
  if (at.label || at.charge || at.radical || at.el === 'C' || at.el === 'H' || !(at.el in VALENCE_ELECTRONS)) return null;
  const typed = at.hcount !== null && at.hcount !== undefined;
  if (!doc.valenceProblem(a) && !typed) return null;
  const total = doc.bondSum(a) + (at.hcount || 0);
  if (!doc.valenceProblem(a) && typed && total === doc.bondSum(a) + doc.usualH(a)) return null;
  for (const z of [1, -1, 2, -2]) {
    at.charge = z;
    try {
      const vals = doc.valences(a);
      if (vals && vals.includes(total)) return z;
    } finally { at.charge = 0; }
  }
  return null;
}

/** Each atom given the charge neededCharge calls for: {atom: charge}. */
export function fixCharges(doc, ids = null) {
  const done = {};
  for (const a of [...doc.atoms.keys()]) {
    if (ids && !ids.has(a)) continue;
    const z = neededCharge(doc, a);
    if (z === null) continue;
    const at = doc.atoms.get(a), hs = doc.implicitH(a);
    at.charge = z;
    if (at.hcount !== null && at.hcount !== undefined) {   // NH3 typed on one bond: NH3+, whose H are the usual now
      const keep = at.hcount;
      at.hcount = null;
      if (doc.implicitH(a) !== hs) at.hcount = keep;
    }
    done[a] = z;
  }
  return done;
}

/** What is wrong with the structures (all, or those with atoms in ids): valence problems, unknown labels, charges
 * they need (canvas.check_structure, without the desktop-only stereo and curved-arrow checks). */
export function checkStructure(doc, ids = null) {
  const atoms = ids && ids.size ? ids : new Set(doc.atoms.keys());
  const out = [];
  for (const a of atoms) {
    const at = doc.atoms.get(a);
    const name = at.label || at.el;
    const p = doc.valenceProblem(a);
    if (p) out.push(`${name}: ${p}`);
    if (at.label && !NICKNAMES[at.label] && !VALENCE[at.label]) out.push(`Label "${at.label}" is not a known nickname: it counts as C.`);
    const z = neededCharge(doc, a);
    if (z !== null) out.push(`${at.el} looks like it needs a charge (${Math.abs(z) === 1 ? '' : Math.abs(z)}${z > 0 ? '+' : '−'}): Add missing charges.`);
  }
  return [...new Set(out)];
}

// ------------------------------------------------------------------------------------------ objects (canvas.py)
/** The objects of a selection (or of everything): each structure (its atoms), arrow and caption, with its box. */
export function objects(doc, keys) {
  const take = keys ? new Set(keys) : null, has = (k, i) => !take || take.has(`${k}:${i}`);
  const out = [];
  for (const frag of doc.fragments()) {
    const atoms = frag.filter(a => has('atom', a));
    if (!atoms.length) continue;
    const all = doc.fragmentOf(atoms[0]);
    const xs = all.map(a => doc.atoms.get(a).x), ys = all.map(a => doc.atoms.get(a).y);
    out.push({ kind: 'mol', ids: all, box: [Math.min(...xs) - 0.3, Math.min(...ys) - 0.3, Math.max(...xs) + 0.3, Math.max(...ys) + 0.3] });
  }
  for (const [i, ar] of doc.arrows) if (has('arrow', i)) out.push({ kind: 'arrow', ids: [i], box: [Math.min(ar.x1, ar.x2), Math.min(ar.y1, ar.y2), Math.max(ar.x1, ar.x2), Math.max(ar.y1, ar.y2)] });
  for (const [i, t] of doc.texts) if (has('text', i)) out.push({ kind: 'text', ids: [i], box: [t.x, t.y, t.x + String(t.text).length * 0.4, t.y + 0.8] });
  return out;
}

/** Moves one object by (dx, dy). */
export function moveObject(doc, o, dx, dy) {
  if (o.kind === 'mol') for (const a of o.ids) { const at = doc.atoms.get(a); at.x += dx; at.y += dy; }
  if (o.kind === 'arrow') for (const i of o.ids) { const ar = doc.arrows.get(i); ar.x1 += dx; ar.x2 += dx; ar.y1 += dy; ar.y2 += dy; }
  if (o.kind === 'text') for (const i of o.ids) { const t = doc.texts.get(i); t.x += dx; t.y += dy; }
}

/** The objects lined up: 'left' | 'hcenter' | 'right' | 'top' | 'vcenter' | 'bottom' (canvas.align). False: fewer than two. */
export function align(doc, keys, how) {
  const objs = objects(doc, keys);
  if (objs.length < 2) return false;
  const at = (b, h) => ({ left: b[0], right: b[2], top: b[1], bottom: b[3], hcenter: (b[0] + b[2]) / 2, vcenter: (b[1] + b[3]) / 2 })[h];
  const vals = objs.map(o => at(o.box, how));
  const target = how === 'left' || how === 'top' ? Math.min(...vals) : how === 'right' || how === 'bottom' ? Math.max(...vals)
    : vals.reduce((s, v) => s + v, 0) / vals.length;
  const horizontal = ['left', 'right', 'hcenter'].includes(how);
  objs.forEach((o, k) => moveObject(doc, o, horizontal ? target - vals[k] : 0, horizontal ? 0 : target - vals[k]));
  return true;
}

/** Structures and straight arrows in one row, left to right as they are, centred on one line with even gaps; an arrow
 * three bond lengths long at least, its captions kept with it (canvas.arrange_scheme). False: fewer than two. */
export function arrangeScheme(doc, keys) {
  const objs = objects(doc, keys).filter(o => o.kind === 'mol' || (o.kind === 'arrow' && !String(doc.arrows.get(o.ids[0]).kind || '').startsWith('curved')));
  if (objs.length < 2) return false;
  objs.sort((p, q) => (p.box[0] + p.box[2]) - (q.box[0] + q.box[2]));
  const mols = objs.filter(o => o.kind === 'mol');
  const line = mols.length ? mols.reduce((s, o) => s + (o.box[1] + o.box[3]) / 2, 0) / mols.length : objs[0].box[1];
  const gap = 0.9;
  let cursor = objs[0].box[0];
  for (const o of objs) {
    if (o.kind === 'arrow') {
      const ar = doc.arrows.get(o.ids[0]);
      const len = Math.max(3, Math.abs(ar.x2 - ar.x1));
      const captions = [...doc.texts.values()].filter(t => t.x >= Math.min(ar.x1, ar.x2) - 0.5 && t.x <= Math.max(ar.x1, ar.x2) + 0.5
        && Math.abs(t.y - (ar.y1 + ar.y2) / 2) < 1.6);
      const ox = (ar.x1 + ar.x2) / 2, oy = (ar.y1 + ar.y2) / 2;
      const forward = ar.x2 >= ar.x1;
      ar.x1 = forward ? cursor : cursor + len; ar.x2 = forward ? cursor + len : cursor; ar.y1 = ar.y2 = line;
      for (const t of captions) { t.x += (cursor + len / 2) - ox; t.y += line - oy; }
      cursor += len + gap;
    } else {
      moveObject(doc, o, cursor - o.box[0], line - (o.box[1] + o.box[3]) / 2);
      cursor += o.box[2] - o.box[0] + gap;
    }
  }
  return true;
}
