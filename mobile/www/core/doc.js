// doc.js - the drawing, the same as the desktop's model.py: atoms and bonds with the geometry the editor needs, in bond
// lengths with y down. The file is the desktop's .sed JSON; what the phone cannot show (shapes, grids, sequences, TLC
// plates, ring fills) is kept as it is and saved back.
import { NICKNAMES, VALENCE, VALENCE_ELECTRONS } from './data.js';

export const MERGE = 0.3;
const MOL_BOND = 1.5;
const OTHER_KINDS = ['shapes', 'grids', 'fills', 'seqs', 'tlcs', 'groups'];
const TAU = 2 * Math.PI;
const mod = (a, n) => ((a % n) + n) % n;
export const snapAngle = (a, deg = 15) => { const s = deg * Math.PI / 180; return Math.round(a / s) * s; };

export function segDist(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1, dy = y2 - y1, L2 = dx * dx + dy * dy;
  const t = L2 === 0 ? 0 : Math.max(0, Math.min(1, ((px - x1) * dx + (py - y1) * dy) / L2));
  return Math.hypot(px - (x1 + t * dx), py - (y1 + t * dy));
}

export class Doc {
  constructor() {
    this.atoms = new Map(); this.bonds = new Map(); this.arrows = new Map(); this.texts = new Map();
    this.other = {};                                       // the kinds the phone only keeps
    this.style = {}; this.next = 1; this.header = {};
  }

  id() { return this.next++; }
  addAtom(x, y, el = 'C', charge = 0) {
    const i = this.id();
    this.atoms.set(i, { x, y, el, charge, label: '', color: '', highlight: '', charge_pos: null });
    return i;
  }
  atomAtOrNew(x, y, el = 'C', exclude = []) {
    const hit = this.findAtom(x, y, MERGE, exclude);
    return hit !== null ? hit : this.addAtom(x, y, el);
  }
  bondBetween(a, b) {
    for (const [i, bd] of this.bonds) if ((bd.a === a && bd.b === b) || (bd.a === b && bd.b === a)) return i;
    return null;
  }
  addBond(a, b, order = 1, stereo = '') {
    if (a === b) return null;
    const old = this.bondBetween(a, b);
    if (old !== null) return old;
    const i = this.id();
    this.bonds.set(i, { a, b, order, stereo, side: null, style: '', color: '', highlight: '' });
    return i;
  }
  addArrow(x1, y1, x2, y2, kind = 'forward', bend = 0) {
    const i = this.id();
    this.arrows.set(i, { x1, y1, x2, y2, kind, bend, color: '' });
    return i;
  }
  addText(x, y, text) {
    const i = this.id();
    this.texts.set(i, { x, y, text, align: 'left', size: 0, bold: false, italic: false, formula: false, color: '' });
    return i;
  }
  neighbors(a) {
    const out = [];
    for (const [i, bd] of this.bonds) {
      if (bd.a === a) out.push([bd.b, i]);
      else if (bd.b === a) out.push([bd.a, i]);
    }
    return out;
  }
  removeBond(i, prune = true) {
    const bd = this.bonds.get(i);
    this.bonds.delete(i);
    if (bd && prune) {
      for (const a of [bd.a, bd.b]) {
        const at = this.atoms.get(a);
        if (at && at.el === 'C' && !at.charge && !this.neighbors(a).length) this.atoms.delete(a);
      }
    }
  }
  removeAtom(a) {
    for (const [, b] of this.neighbors(a)) this.removeBond(b, false);
    this.atoms.delete(a);
    const fills = this.other.fills || [];
    this.other.fills = fills.filter(f => !f.atoms.includes(a));
  }
  stripAtom(a) {                                           // a labelled atom turns back into a plain carbon first
    const at = this.atoms.get(a);
    if (at.el === 'C' && !at.label && !at.charge) return false;
    at.el = 'C'; at.label = ''; at.charge = 0;
    return true;
  }
  mergeAtom(src, dst) {
    for (const [n, b] of this.neighbors(src)) {
      const bd = this.bonds.get(b);
      if (n === dst || this.bondBetween(n, dst) !== null) { this.bonds.delete(b); continue; }
      if (bd.a === src) bd.a = dst; else bd.b = dst;
    }
    this.atoms.delete(src);
  }
  isEmpty() {
    return !this.atoms.size && !this.arrows.size && !this.texts.size && !this.otherCount();
  }
  otherCount() { return OTHER_KINDS.filter(k => k !== 'fills').reduce((n, k) => n + (this.other[k] || []).length, 0); }

  // ---------------------------------------------------------------------------------------- hit testing
  findAtom(x, y, r, exclude = []) {
    let best = null, bd = r;
    for (const [i, a] of this.atoms) {
      if (exclude.includes(i)) continue;
      const d = Math.hypot(a.x - x, a.y - y);
      if (d <= bd) { best = i; bd = d; }
    }
    return best;
  }
  findBond(x, y, r) {
    let best = null, bd = r;
    for (const [i, b] of this.bonds) {
      const p = this.atoms.get(b.a), q = this.atoms.get(b.b);
      const d = segDist(x, y, p.x, p.y, q.x, q.y);
      if (d <= bd) { best = i; bd = d; }
    }
    return best;
  }
  findArrow(x, y, r) {
    let best = null, bd = r;
    for (const [i, a] of this.arrows) {
      const d = segDist(x, y, a.x1, a.y1, a.x2, a.y2);
      if (d <= bd) { best = i; bd = d; }
    }
    return best;
  }

  // ---------------------------------------------------------------------------------------- chemistry
  bondSum(a) { return this.neighbors(a).reduce((s, [, b]) => s + this.bonds.get(b).order, 0); }
  valences(a) {
    const at = this.atoms.get(a), base = VALENCE[at.el], r = at.radical || 0;
    if (!base) return null;
    if (!at.charge) return r ? base.map(v => v - r).filter(v => v >= 0) : base;
    if (at.el === 'H') return [0];
    let e = VALENCE_ELECTRONS[at.el];
    if (e === undefined) return base.map(v => v - Math.abs(at.charge)).filter(v => v >= 0);
    e -= at.charge;
    return [Math.max(0, Math.min(e, 8 - e) - r)];        // each unpaired electron (•) takes a bond away
  }
  implicitH(a) {
    const at = this.atoms.get(a);
    if (at.label) return 0;
    if (at.hcount !== null && at.hcount !== undefined) return at.hcount;   // typed: CH2, NH
    const vals = this.valences(a);
    if (!vals || !vals.length) return 0;
    const s = this.bondSum(a);
    for (const v of vals) if (v >= s) return v - s;
    return 0;
  }
  valenceProblem(a) {
    const at = this.atoms.get(a), bonds = this.bondSum(a);
    if (at.label) return bonds > 1 && NICKNAMES[at.label] ? `${at.label} attaches by one bond, but has ${bonds}.` : '';
    const vals = this.valences(a);
    if (at.hcount !== null && at.hcount !== undefined) {
      return vals && ['C', 'N', 'O', 'F', 'B'].includes(at.el) && bonds + at.hcount > 4
        ? `${at.el} with ${at.hcount} H has room for ${Math.max(0, 4 - at.hcount)} bonds, but has ${bonds}.` : '';
    }
    if (vals && vals.length && bonds > Math.max(...vals)) {
      const m = Math.max(...vals);
      const q = at.charge ? ` ${at.charge > 0 ? '+' : '−'}${Math.abs(at.charge)}` : '';
      return `${at.el}${q} takes at most ${m} bond${m !== 1 ? 's' : ''}, but has ${bonds}.`;
    }
    return '';
  }
  fragments(ids = null) {
    const todo = new Set(ids ? [...ids].filter(a => this.atoms.has(a)) : this.atoms.keys());
    const out = [];
    while (todo.size) {
      const seed = Math.min(...todo), comp = new Set([seed]), stack = [seed];
      while (stack.length) {
        for (const [n] of this.neighbors(stack.pop())) if (todo.has(n) && !comp.has(n)) { comp.add(n); stack.push(n); }
      }
      for (const c of comp) todo.delete(c);
      out.push([...comp].sort((a, b) => a - b));
    }
    return out;
  }
  fragmentOf(a) { return this.fragments().find(f => f.includes(a)) || []; }

  /** Ends of curved arrows tied to an atom or bond on the desktop (tail_on / head_on: [kind, id, dx, dy]) follow it. */
  updateAttachedArrows() {
    for (const ar of this.arrows.values()) {
      for (const end of ['tail', 'head']) {
        const on = ar[end + '_on'];
        if (!on) continue;
        let p = null;
        if (on[0] === 'atom' && this.atoms.has(on[1])) { const at = this.atoms.get(on[1]); p = [at.x, at.y]; }
        if (on[0] === 'bond' && this.bonds.has(on[1])) {
          const bd = this.bonds.get(on[1]), u = this.atoms.get(bd.a), v = this.atoms.get(bd.b);
          p = [(u.x + v.x) / 2, (u.y + v.y) / 2];
        }
        if (!p) { ar[end + '_on'] = null; continue; }
        if (end === 'tail') { ar.x1 = p[0] + on[2]; ar.y1 = p[1] + on[3]; } else { ar.x2 = p[0] + on[2]; ar.y2 = p[1] + on[3]; }
      }
    }
  }

  // ---------------------------------------------------------------------------------------- groups (desktop Object > Group)
  groupOf(kind, i) {
    if (kind === 'bond') { const bd = this.bonds.get(i); if (!bd) return null; kind = 'atom'; i = bd.a; }
    for (const g of this.other.groups || []) if (g.members.some(([k, j]) => k === kind && j === i)) return g;
    return null;
  }
  /** Keys 'kind:id' of everything in the groups the given items belong to, plus the items. */
  withGroups(keys) {
    const out = new Set(keys);
    for (const key of keys) {
      const [kind, id] = key.split(':');
      const g = this.groupOf(kind, Number(id));
      if (!g) continue;
      const atoms = new Set();
      for (const [k, j] of g.members) { out.add(`${k}:${j}`); if (k === 'atom') atoms.add(j); }
      for (const [b, bd] of this.bonds) if (atoms.has(bd.a) && atoms.has(bd.b)) out.add(`bond:${b}`);
    }
    return out;
  }
  /** Forgets deleted members (atoms bonded to a grouped structure join it); drops groups left with under two objects. */
  pruneGroups() {
    const groups = this.other.groups || [];
    if (!groups.length) return;
    const table = { atom: this.atoms, bond: this.bonds, arrow: this.arrows, text: this.texts };
    const exists = (k, j) => (table[k] ? table[k].has(j) : (this.other[k + 's'] || []).some(o => o.id === j));
    const frags = this.fragments();
    this.other.groups = groups.filter(g => {
      const members = g.members.filter(([k, j]) => exists(k, j));
      const atoms = new Set(members.filter(([k]) => k === 'atom').map(([, j]) => j));
      for (const f of frags) if (f.some(a => atoms.has(a))) for (const a of f) if (!atoms.has(a)) { atoms.add(a); members.push(['atom', a]); }
      const pieces = frags.filter(f => f.some(a => atoms.has(a))).length + members.filter(([k]) => k !== 'atom').length;
      g.members = members;
      return pieces >= 2;
    });
  }

  // ---------------------------------------------------------------------------------------- geometry
  bondAngles(a) {
    const at = this.atoms.get(a);
    return this.neighbors(a).map(([n]) => { const o = this.atoms.get(n); return Math.atan2(o.y - at.y, o.x - at.x); });
  }
  nextBondAngle(a) {
    const nb = this.neighbors(a), at = this.atoms.get(a);
    if (!nb.length) return -30 * Math.PI / 180;
    const angles = this.bondAngles(a);
    if (nb.length === 1) {
      const base = angles[0];
      if (nb.some(([, b]) => this.bonds.get(b).order === 3)) return base + Math.PI;
      const cands = [base + 2 * Math.PI / 3, base - 2 * Math.PI / 3];
      const room = ang => {
        const x = at.x + Math.cos(ang), y = at.y + Math.sin(ang);
        let m = 9;
        for (const [i, o] of this.atoms) if (i !== a && i !== nb[0][0]) m = Math.min(m, Math.hypot(o.x - x, o.y - y));
        return m;
      };
      const r = cands.map(room);
      if (Math.abs(r[0] - r[1]) > 1e-6) return r[0] > r[1] ? cands[0] : cands[1];
      return Math.abs(Math.sin(cands[0])) <= Math.abs(Math.sin(cands[1])) ? cands[0] : cands[1];
    }
    angles.sort((p, q) => p - q);
    let gap = -1, start = 0;
    angles.forEach((ang, k) => {
      const g = mod(angles[(k + 1) % angles.length] - ang, TAU) || TAU;
      if (g > gap) { gap = g; start = ang; }
    });
    return start + gap / 2;
  }
  freeDirection(a) {
    const at = this.atoms.get(a);
    let sx = 0, sy = 0;
    for (const [n] of this.neighbors(a)) {
      const o = this.atoms.get(n), d = Math.hypot(o.x - at.x, o.y - at.y) || 1;
      sx += (o.x - at.x) / d; sy += (o.y - at.y) / d;
    }
    const d = Math.hypot(sx, sy);
    if (d < 1e-6) return this.neighbors(a).length ? [0, -1] : [1, 0];
    return [-sx / d, -sy / d];
  }
  ringOfBond(i, limit = 8) {
    const bd = this.bonds.get(i), prev = new Map([[bd.a, null]]);
    let frontier = [bd.a];
    for (let depth = 0; depth < limit - 1; depth++) {
      const nxt = [];
      for (const u of frontier) {
        for (const [n, b] of this.neighbors(u)) {
          if (b === i || prev.has(n)) continue;
          prev.set(n, u);
          if (n === bd.b) {
            const path = [n];
            while (prev.get(path[path.length - 1]) !== null) path.push(prev.get(path[path.length - 1]));
            return path;
          }
          nxt.push(n);
        }
      }
      frontier = nxt;
    }
    return null;
  }
  doubleSide(i) {
    const bd = this.bonds.get(i);
    if (bd.side !== null && bd.side !== undefined) return bd.side;
    const p = this.atoms.get(bd.a), q = this.atoms.get(bd.b);
    const nx = -(q.y - p.y), ny = q.x - p.x;
    const sideOf = (x, y) => { const v = (x - p.x) * nx + (y - p.y) * ny; return Math.abs(v) < 1e-6 ? 0 : v > 0 ? 1 : -1; };
    const ring = this.ringOfBond(i);
    if (ring) {
      const cx = ring.reduce((s, k) => s + this.atoms.get(k).x, 0) / ring.length;
      const cy = ring.reduce((s, k) => s + this.atoms.get(k).y, 0) / ring.length;
      const s = sideOf(cx, cy);
      if (s) return s;
    }
    const na = this.neighbors(bd.a).map(([n]) => n).filter(n => n !== bd.b);
    const nb = this.neighbors(bd.b).map(([n]) => n).filter(n => n !== bd.a);
    if (!na.length && !nb.length) return 0;
    if ((!na.length && nb.length !== 1) || (!nb.length && na.length !== 1)) return 0;
    const total = [...na, ...nb].reduce((s, n) => s + sideOf(this.atoms.get(n).x, this.atoms.get(n).y), 0);
    return total === 0 ? 0 : total > 0 ? 1 : -1;
  }
  ringPoints(n, cx, cy, start) {
    const R = 1 / (2 * Math.sin(Math.PI / n));
    return Array.from({ length: n }, (_, k) => [cx + R * Math.cos(start + TAU * k / n), cy + R * Math.sin(start + TAU * k / n)]);
  }
  ringOnBond(i, n) {
    const bd = this.bonds.get(i), p = this.atoms.get(bd.a), q = this.atoms.get(bd.b);
    const mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2, L = Math.hypot(q.x - p.x, q.y - p.y) || 1;
    const nx = -(q.y - p.y) / L, ny = (q.x - p.x) / L;
    let score = 0;
    for (const a of [bd.a, bd.b]) {
      for (const [o] of this.neighbors(a)) {
        if (o === bd.a || o === bd.b) continue;
        const at = this.atoms.get(o);
        score += (at.x - mx) * nx + (at.y - my) * ny;
      }
    }
    const s = score > 0 ? -1 : 1, ap = 1 / (2 * Math.tan(Math.PI / n));
    const cx = mx + s * nx * ap, cy = my + s * ny * ap;
    return this.ringPoints(n, cx, cy, Math.atan2(p.y - cy, p.x - cx));
  }
  ringOnAtom(a, n) {
    const at = this.atoms.get(a);
    const ang = this.neighbors(a).length ? this.nextBondAngle(a) : Math.PI / 2;
    const R = 1 / (2 * Math.sin(Math.PI / n));
    return this.ringPoints(n, at.x + R * Math.cos(ang), at.y + R * Math.sin(ang), ang + Math.PI);
  }
  addRing(pts, aromatic = false) {
    const ids = pts.map(([x, y]) => this.atomAtOrNew(x, y)), n = ids.length;
    const edges = ids.map((u, k) => [u, ids[(k + 1) % n]]);
    const bonds = edges.map(([u, v]) => this.addBond(u, v));
    if (aromatic) {
      const doubles = par => { const out = []; for (let k = par; k < n - n % 2; k += 2) out.push(k); return out; };
      const hasDouble = (atom, skip) => this.neighbors(atom).some(([, b]) => !skip.includes(b) && this.bonds.get(b).order >= 2);
      const conflicts = par => {
        let c = 0;
        const ds = doubles(par);
        for (const k of ds) {
          const [u, v] = edges[k];
          if (this.bonds.get(bonds[k]).order === 1 && (hasDouble(u, bonds) || hasDouble(v, bonds))) c++;
        }
        for (let k = 0; k < n; k++) if (!ds.includes(k) && this.bonds.get(bonds[k]).order >= 2) c++;
        return c;
      };
      const par = conflicts(1) < conflicts(0) ? 1 : 0;
      for (const k of doubles(par)) {
        const [u, v] = edges[k], b = this.bonds.get(bonds[k]);
        if (b.order === 1 && !hasDouble(u, bonds) && !hasDouble(v, bonds)) b.order = 2;
      }
    }
    return ids;
  }
  chainPoints(x0, y0, x1, y1) {
    const ang = snapAngle(Math.atan2(y1 - y0, x1 - x0), 30), dist = Math.hypot(x1 - x0, y1 - y0);
    const step = Math.sqrt(3) / 2, n = Math.max(1, Math.round(dist / step));
    const ux = Math.cos(ang), uy = Math.sin(ang), vx = -uy, vy = ux;
    let sgn = vy > 0 ? -1 : 1;
    if (Math.abs(uy) > 0.99) sgn = vx > 0 ? 1 : -1;
    return Array.from({ length: n }, (_, j) => {
      const k = j + 1;
      return [x0 + k * step * ux + (k % 2) * 0.5 * vx * sgn, y0 + k * step * uy + (k % 2) * 0.5 * vy * sgn];
    });
  }

  // ---------------------------------------------------------------------------------------- files
  toDict() {
    const items = m => [...m].map(([id, v]) => ({ id, ...v }));
    const out = { format: 'structure-editor', version: 2, style: { ...this.style }, atoms: items(this.atoms), bonds: items(this.bonds),
      arrows: items(this.arrows), texts: items(this.texts) };
    for (const k of OTHER_KINDS) out[k] = this.other[k] || [];
    return out;
  }
  static fromDict(d) {
    const doc = new Doc();
    doc.style = { ...(d.style || {}) };
    for (const [key, map] of [['atoms', doc.atoms], ['bonds', doc.bonds], ['arrows', doc.arrows], ['texts', doc.texts]]) {
      for (const { id, ...v } of d[key] || []) { map.set(id, v); doc.next = Math.max(doc.next, id + 1); }
    }
    for (const k of OTHER_KINDS) {
      doc.other[k] = d[k] || [];
      for (const item of doc.other[k]) if (item.id) doc.next = Math.max(doc.next, item.id + 1);
    }
    return doc;
  }
  copy() { return Doc.fromDict(JSON.parse(JSON.stringify(this.toDict()))); }

  /** V2000 MOL text of the atoms in ids (all when null), as model.to_molblock: [text, atom ids in file order]. */
  toMolblock(ids = null, title = '') {
    const order = [...this.atoms.keys()].filter(a => !ids || ids.has(a));
    const index = new Map(order.map((a, k) => [a, k + 1]));
    const bonds = [...this.bonds.values()].filter(b => index.has(b.a) && index.has(b.b));
    return [buildBlock(order.map(a => this.atoms.get(a)), bonds.map(b => ({ a: index.get(b.a) - 1, b: index.get(b.b) - 1, order: b.order, stereo: b.stereo })), title), order];
  }
}

/** A V2000 block from plain atoms {x, y (down, bond lengths), el, charge} and bonds {a, b (0-based), order, stereo}. */
export function buildBlock(atoms, bonds, title = '') {
  const f = (v, w) => v.toFixed(4).padStart(w);
  const n3 = v => String(v).padStart(3);
  const lines = [title, '  Kekulex 2D', '', `${n3(atoms.length)}${n3(bonds.length)}  0  0  0  0  0  0  0  0999 V2000`];
  for (const at of atoms) lines.push(`${f(at.x * MOL_BOND, 10)}${f(-at.y * MOL_BOND, 10)}${f(0, 10)} ${at.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  const radicals = atoms.map((a, k) => [k + 1, a.radical === 1 ? 2 : a.radical ? 3 : 0]).filter(r => r[1]);
  for (const b of bonds) lines.push(`${n3(b.a + 1)}${n3(b.b + 1)}${n3(b.order)}${n3({ wedge: 1, hash: 6 }[b.stereo] || 0)}`);
  const charged = atoms.map((a, k) => [k + 1, a.charge]).filter(c => c[1]);
  for (let k = 0; k < charged.length; k += 8) {
    const chunk = charged.slice(k, k + 8);
    lines.push(`M  CHG${n3(chunk.length)}` + chunk.map(([i, c]) => String(i).padStart(4) + String(c).padStart(4)).join(''));
  }
  for (let k = 0; k < radicals.length; k += 8) {            // 2 doublet (•), 3 triplet (carbene ••)
    const chunk = radicals.slice(k, k + 8);
    lines.push(`M  RAD${n3(chunk.length)}` + chunk.map(([i, c]) => String(i).padStart(4) + String(c).padStart(4)).join(''));
  }
  lines.push('M  END');
  return lines.join('\n') + '\n';
}

/** Atoms {el, charge, x, y (down)} and bonds {a, b, order, stereo} of a V2000 block. */
export function parseBlock(block) {
  const lines = block.split(/\r?\n/);
  const n = parseInt(lines[3].slice(0, 3), 10), m = parseInt(lines[3].slice(3, 6), 10);
  if (!(n >= 0) || !lines[3].includes('V2000')) throw new Error('Only V2000 MOL files can be read here.');
  const atoms = [], bonds = [];
  for (let k = 0; k < n; k++) {
    const l = lines[4 + k];
    atoms.push({ x: parseFloat(l.slice(0, 10)), y: -parseFloat(l.slice(10, 20)), el: l.slice(31, 34).trim(), charge: 0, radical: 0 });
  }
  for (let k = 0; k < m; k++) {
    const l = lines[4 + n + k];
    const code = parseInt(l.slice(9, 12), 10) || 0;
    bonds.push({ a: parseInt(l.slice(0, 3), 10) - 1, b: parseInt(l.slice(3, 6), 10) - 1, order: parseInt(l.slice(6, 9), 10),
      stereo: code === 1 ? 'wedge' : code === 6 ? 'hash' : '' });
  }
  for (const l of lines.slice(4 + n + m)) {
    if (l.startsWith('M  CHG')) {
      const count = parseInt(l.slice(6, 9), 10);
      for (let k = 0; k < count; k++) atoms[parseInt(l.slice(9 + 8 * k, 13 + 8 * k), 10) - 1].charge = parseInt(l.slice(13 + 8 * k, 17 + 8 * k), 10);
    }
    if (l.startsWith('M  RAD')) {
      const count = parseInt(l.slice(6, 9), 10);
      for (let k = 0; k < count; k++) {
        const v = parseInt(l.slice(13 + 8 * k, 17 + 8 * k), 10);
        atoms[parseInt(l.slice(9 + 8 * k, 13 + 8 * k), 10) - 1].radical = v === 2 ? 1 : v ? 2 : 0;
      }
    }
    if (l.startsWith('M  END')) break;
  }
  return { atoms, bonds };
}
