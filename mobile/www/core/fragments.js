// fragments.js - the mass-spectrum fragments a structure can logically give: the desktop's fragments.py (predict),
// rule for rule, on RDKit.js. The rule tables, mechanism texts and literature come from the desktop (data.js FRAG).
// Drawing a fragment or its mechanism (draw_doc, mechanism, the Lewis clean-up) stays on the desktop.
import { ELECTRON, ELEMENTS, FRAG } from './data.js';
import { chargeText, expanded, formulaText } from './chem.js';

export const MODES = FRAG.MODES;
export const LIKELY = FRAG.LIKELY, POSSIBLE = FRAG.POSSIBLE;
const SYMBOLS = Object.keys(ELEMENTS);                     // by atomic number: H = 1
const massOf = el => (ELEMENTS[el] ? ELEMENTS[el][0] : 0);

export const ionText = (counts, z, odd) => formulaText(counts) + (odd ? '•' : '') + chargeText(z);

function sub(a, b) {
  const out = { ...a };
  for (const [k, v] of Object.entries(b)) out[k] = (out[k] || 0) - v;
  return out;
}
function clean(c) {
  const out = {};
  for (const [k, v] of Object.entries(c)) if (v) out[k] = v;
  return out;
}
const anyValue = o => Object.values(o).some(v => v);
const sorted = s => [...s].sort((p, q) => p - q);
const sameSet = (a, b) => a.size === b.size && [...a].every(x => b.has(x));
const pairKey = (i, j) => (i < j ? `${i},${j}` : `${j},${i}`);
const sameCounts = (a, b) => { const ka = Object.keys(a), kb = Object.keys(b); return ka.length === kb.length && ka.every(k => a[k] === b[k]); };

export class Fragment {
  constructor(mz, z, odd, counts, ion, lost, how, score, step = 1) {
    Object.assign(this, { mz, z, odd, counts, ion, lost, how, score, step });
    this.atoms = new Set(); this.lostAtoms = new Set(); this.also = [];
    this.S = null; this.marks = null; this.mech = ''; this.ref = ''; this.parent = null;
    this._idx = null; this._imarks = null;
  }
  get likelihood() { return this.score >= LIKELY ? 'likely' : this.score >= POSSIBLE ? 'possible' : 'minor'; }
  get formula() { return formulaText(this.counts); }
}

/** The fragments from the precursor's first product to frag, in order. */
export function chainOf(frag) {
  const out = [];
  for (let f = frag; f; f = f.parent) out.unshift(f);
  return out;
}
export const stepLoss = f => (f.lost ? f.lost.split(' + ').pop() : '');

// ------------------------------------------------------------------------------------------- the molecule
class Mol {
  constructor(RD, doc, atomIds) {
    const [full, grown] = expanded(RD, doc);
    this.owner = new Map([...doc.atoms.keys()].map(a => [a, a]));
    for (const [a, more] of grown) for (const g of more) this.owner.set(g, a);
    const ids = [...doc.atoms.keys()].filter(a => !atomIds || atomIds.has(a));
    const keep = new Set(ids);
    for (const a of ids) for (const g of grown.get(a) || []) keep.add(g);
    const sel = new Set([...keep].filter(a => full.atoms.has(a)));
    const [block, order] = full.toMolblock(sel);
    this.order = order;
    this.mol = order.length ? RD.get_mol(block, JSON.stringify({ removeHs: false })) : null;
    if (!this.mol || !this.mol.is_valid()) { this.mol = null; return; }
    const j = JSON.parse(this.mol.get_json());
    const mol = j.molecules[0], da = j.defaults.atom, db = j.defaults.bond;
    const ext = (mol.extensions || []).find(e => e.name === 'rdkitRepresentation') || {};
    const atoms = mol.atoms.map(a => ({ ...da, ...a })), bonds = (mol.bonds || []).map(b => ({ ...db, ...b }));
    this.n = atoms.length;
    this.el = atoms.map(a => SYMBOLS[a.z - 1]);
    this.key = this.el;
    this.hs = atoms.map(a => a.impHs);
    this.q = atoms.map(a => a.chg);
    this.rad = atoms.map(a => a.nRad);
    this.nb = atoms.map(() => []);
    this.bo = new Map();
    for (const b of bonds) {
      const [a, c] = b.atoms;
      this.nb[a].push(c); this.nb[c].push(a);
      this.bo.set(pairKey(a, c), b.bo);
    }
    this.rings = ext.atomRings || [];
    this.ringBond = new Set();
    for (const r of this.rings) for (let k = 0; k < r.length; k++) this.ringBond.add(pairKey(r[k], r[(k + 1) % r.length]));
    const aroAtom = new Set(ext.aromaticAtoms || []);
    // aromatic as a chemist sees it: rings of aromatic atoms without an exocyclic C=O
    this.arobond = new Set();
    const aro = new Set();
    for (const ring of this.rings) {
      if (!ring.every(i => aroAtom.has(i))) continue;
      if (ring.some(i => this.nb[i].some(k => ['O', 'S', 'N'].includes(this.el[k]) && !ring.includes(k) && this.bo.get(pairKey(i, k)) === 2))) continue;
      for (const i of ring) aro.add(i);
      for (let k = 0; k < ring.length; k++) this.arobond.add(pairKey(ring[k], ring[(k + 1) % ring.length]));
    }
    this.aro = this.el.map((_e, i) => aro.has(i));
    this.sides = new Map();
  }

  dispose() { if (this.mol) this.mol.delete(); this.mol = null; }

  matches(smarts) {
    let q = null;
    try {
      q = this.RD.get_qmol(smarts);
      if (!q || !q.is_valid()) return [];
      const r = JSON.parse(this.mol.get_substruct_matches(q) || '[]');
      return Array.isArray(r) ? r.map(x => x.atoms) : [];
    } catch { return []; } finally { if (q) q.delete(); }
  }

  hasBond(i, j) { return this.bo.has(pairKey(i, j)); }
  aromaticBond(i, j) { return this.arobond.has(pairKey(i, j)); }
  single(i, j, ring = false) {
    const k = pairKey(i, j);
    return this.bo.has(k) && this.bo.get(k) === 1 && (ring || !this.ringBond.has(k)) && !this.arobond.has(k);
  }
  double(i, j) { const k = pairKey(i, j); return this.bo.get(k) === 2 && !this.arobond.has(k); }
  bondOrder(i, j) {
    if (this.aromaticBond(i, j)) return null;
    const o = this.bo.get(pairKey(i, j));
    return o === 1 || o === 2 || o === 3 ? o : null;
  }
  h(i) { return this.hs[i] + this.nb[i].filter(j => this.el[j] === 'H').length; }
  heavy(i) { return this.nb[i].filter(j => this.el[j] !== 'H'); }

  /** The atoms reached from i without crossing bond i-j (and the bonds in cut); null when j is reached. */
  side(i, j, cut = []) {
    const key = `${i}|${j}|${cut.map(c => pairKey(c[0], c[1])).sort().join(';')}`;
    if (this.sides.has(key)) return this.sides.get(key);
    const blocked = new Set([pairKey(i, j), ...cut.map(c => pairKey(c[0], c[1]))]);
    const seen = new Set([i]), todo = [i];
    while (todo.length) {
      const x = todo.pop();
      for (const y of this.nb[x]) if (!seen.has(y) && !blocked.has(pairKey(x, y))) { seen.add(y); todo.push(y); }
    }
    const out = seen.has(j) ? null : seen;
    this.sides.set(key, out);
    return out;
  }

  cutOrders(i, S) {
    let out = 0;
    for (const j of this.nb[i]) {
      if (S.has(j)) continue;
      const o = this.bondOrder(i, j);
      if (o === null) return null;
      out += o;
    }
    return out;
  }
  hcap(S) {
    let total = 0;
    for (const i of S) { const c = this.cutOrders(i, S); if (c === null) return null; total += this.hs[i] + c; }
    return total;
  }
  counts(S, h) {
    const out = {};
    for (const i of sorted(S)) out[this.key[i]] = (out[this.key[i]] || 0) + 1;
    if (h) out.H = (out.H || 0) + h;
    return out;
  }

  carbonylO(c) {
    if (this.el[c] !== 'C' || this.aro[c]) return null;
    const o = this.nb[c].find(j => ['O', 'S'].includes(this.el[j]) && this.double(c, j));
    return o === undefined ? null : o;
  }
  saturated(c) { return this.el[c] === 'C' && !this.aro[c] && this.nb[c].every(j => this.single(c, j, true)); }
  vinyl(c) { return this.nb[c].some(j => this.el[j] === 'C' && this.double(c, j) && !this.aro[j]); }
  anomeric(c, x) {
    if (this.el[c] !== 'C' || this.aro[c]) return false;
    for (const o of this.nb[c]) {
      if (o !== x && this.el[o] === 'O' && this.ringBond.has(pairKey(c, o))) {
        for (const r of this.rings) if (r.includes(c) && r.includes(o) && !r.includes(x) && (r.length === 5 || r.length === 6)) return true;
      }
    }
    return false;
  }
  cation(c, S) {
    if (this.el[c] !== 'C') return 0.3;
    if (this.aro[c]) return 0.15;
    if (this.carbonylO(c) !== null) return 0.85;
    if (!this.saturated(c)) return 0.1;
    const near = this.nb[c].filter(j => S.has(j) && this.el[j] !== 'H');
    if (near.some(j => ['N', 'O', 'S'].includes(this.el[j]) && this.single(c, j, true) && this.carbonylO(j) === null)) return 0.7;
    if (near.some(j => this.aro[j])) return 0.75;
    if (near.some(j => this.vinyl(j)) || this.vinyl(c)) return 0.55;
    return { 0: 0.05, 1: 0.15, 2: 0.35 }[near.length] ?? 0.6;
  }
}

// ------------------------------------------------------------------------------------------- the rules
const marksOf = obj => new Map(obj);                       // [[atom, [dq, rad, dh]], ...]
const sumMarks = (marks, k, S = null) => { let s = 0; for (const [a, v] of marks || []) if (!S || S.has(a)) s += v[k]; return s; };

class Engine {
  constructor(m, mode) {
    this.m = m;
    this.mode = mode;
    this.M = m.counts(new Set(m.el.map((_e, i) => i)), m.hs.reduce((s, v) => s + v, 0));
    const drawnZ = m.q.reduce((s, v) => s + v, 0), drawnRad = m.rad.reduce((s, v) => s + v, 0);
    this.charged = !!drawnZ;
    if (this.charged) {
      this.z = drawnZ; this.odd = !!(drawnRad % 2); this.pre = { ...this.M };
      this.label = 'M' + (this.odd ? '•' : '') + chargeText(drawnZ);
      this.mode = drawnZ > 0 ? 'ESI+' : 'ESI-';
      if (this.odd) this.mode = 'EI';
    } else if (mode === 'EI') {
      this.z = 1; this.odd = !(drawnRad % 2); this.pre = { ...this.M }; this.label = '[M]•+';
    } else if (mode === 'ESI+') {
      this.z = 1; this.odd = !!(drawnRad % 2); this.label = '[M+H]+';
      this.pre = { ...this.M, H: (this.M.H || 0) + 1 };
    } else {
      this.z = -1; this.odd = !!(drawnRad % 2); this.label = '[M−H]−';
      this.pre = { ...this.M, H: (this.M.H || 0) - 1 };
    }
    this.out = [];
    this.all = new Set(m.el.map((_e, i) => i));
  }

  mz(counts, z) {
    let mass = 0;
    for (const [k, n] of Object.entries(counts)) mass += massOf(k) * n;
    return (mass - z * ELECTRON) / Math.abs(z);
  }

  place(S, base, site, dq, rad, dh) {
    const m = this.m;
    const marks = new Map();
    for (const [k, v] of base || []) if (S.has(k)) marks.set(k, [...v]);
    if (!marks.has(site)) marks.set(site, [0, 0, 0]);
    marks.get(site)[0] += dq;
    const room = i => { const c = m.cutOrders(i, S); return c === null ? null : m.hs[i] + c + (marks.has(i) ? marks.get(i)[2] : 0); };
    if (rad && dh < 0) {
      const j = m.nb[site].find(j => S.has(j) && !marks.has(j) && m.el[j] === 'C' && room(j));
      if (j !== undefined) { marks.set(j, [0, rad, -1]); rad = 0; dh += 1; }
    }
    marks.get(site)[1] += rad;
    const r = room(site);
    if (r === null) return null;
    if (dh >= 0 || r + dh >= 0) marks.get(site)[2] += dh;
    else {
      marks.get(site)[2] -= r;
      let need = -(r + dh);
      const seen = new Set([site]), todo = [site];
      while (need && todo.length) {
        const x = todo.shift();
        for (const y of m.nb[x]) {
          if (S.has(y) && !seen.has(y)) {
            seen.add(y); todo.push(y);
            const ry = room(y);
            if (ry && !m.aro[y]) {
              const take = Math.min(ry, need, 1);
              if (!marks.has(y)) marks.set(y, [0, 0, 0]);
              marks.get(y)[2] -= take;
              need -= take;
              if (!need) break;
            }
          }
        }
      }
      if (need) return null;
    }
    for (const [k, v] of marks) if (m.aro[k] && v[2] < 0 && !v[0] && !v[1]) return null;
    const out = new Map();
    for (const [k, v] of marks) if (v.some(x => x)) out.set(k, v);
    return out;
  }

  acidSite(S, avoid = []) {
    const m = this.m;
    let best = null;
    for (const o of sorted(S)) {
      if (avoid.includes(o) || m.el[o] !== 'O' || !m.hs[o] || m.heavy(o).length !== 1) continue;
      const c = m.heavy(o)[0];
      const rank = m.carbonylO(c) !== null ? 0 : m.aro[c] ? 1 : 2;
      if (best === null || rank < best[0]) best = [rank, o];
    }
    return best ? best[1] : null;
  }

  baseSite(S, avoid = []) {
    const m = this.m;
    let best = null;
    for (const i of sorted(S)) {
      if (avoid.includes(i) || m.q[i] || !['N', 'O'].includes(m.el[i])) continue;
      const hv = m.heavy(i);
      let rank;
      if (m.el[i] === 'N') {
        if (m.aro[i]) rank = hv.length === 2 && !m.hs[i] ? 1 : 9;
        else if (hv.some(j => m.carbonylO(j) !== null) || !m.nb[i].every(j => m.single(i, j, true))) rank = 9;
        else rank = hv.some(j => m.aro[j]) ? 2 : 0;
      } else if (m.nb[i].some(j => m.double(i, j))) rank = 3;
      else rank = 4;
      if (rank < 9 && (best === null || rank < best[0])) best = [rank, i];
    }
    return best ? best[1] : null;
  }

  make(S, marks, odd, how, score, step = 1, lostName = null, parent = null, mech = null) {
    const m = this.m;
    S = new Set(S);
    if (!S.size || (S.size === m.n && (!marks || !marks.size))) return null;
    const hc = m.hcap(S);
    if (hc === null || !marks) return null;
    const h = hc + sumMarks(marks, 2);
    for (const [i, v] of marks) if (!S.has(i) || m.hs[i] + m.cutOrders(i, S) + v[2] < 0) return null;
    let z = 0;
    for (const i of S) z += m.q[i];
    z += sumMarks(marks, 0);
    if (z !== this.z || h < 0) return null;
    const counts = clean(m.counts(S, h));
    if (parent && (!counts.C || Object.entries(counts).filter(([k]) => k !== 'H').reduce((s, [, v]) => s + v, 0) < 2)) return null;
    const start = parent ? parent.counts : this.pre;
    let lost = sub(start, counts);
    if (Object.values(lost).some(v => v < 0) || !anyValue(lost)) return null;
    lost = clean(lost);
    const wasOdd = parent ? parent.odd : this.odd;
    const name = lostName || (formulaText(lost) + (wasOdd !== odd ? '•' : ''));
    const total = parent ? `${parent.lost} + ${name}` : name;
    return this.frag(counts, odd, total, how, score, step, S, marks, mech, parent);
  }

  formulaOnly(lost, radical, how, score, highlight, parent = null, step = 1, mech = null) {
    const start = parent ? parent.counts : this.pre;
    let counts = sub(start, lost);
    if (Object.values(counts).some(v => v < 0)) return null;
    counts = clean(counts);
    if (!counts.C && !counts.N) return null;
    const odd = (parent ? parent.odd : this.odd) !== radical;
    const name = formulaText(lost) + (radical ? '•' : '');
    const total = parent ? `${parent.lost} + ${name}` : name;
    const f = this.frag(counts, odd, total, how, score, step, null, null, mech, parent);
    f.atoms = new Set(highlight);
    return f;
  }

  fixedIon(counts, odd, how, score, highlight, mech = null) {
    const lost = sub(this.pre, counts);
    if (Object.values(lost).some(v => v < 0) || !anyValue(lost)) return null;
    const name = formulaText(clean(lost)) + (odd !== this.odd ? '•' : '');
    const f = this.frag({ ...counts }, odd, name, how, score, 1, null, null, mech, null);
    f.atoms = new Set(highlight);
    return f;
  }

  frag(counts, odd, lost, how, score, step, S, marks, mech = null, parent = null) {
    const m = this.m;
    const f = new Fragment(this.mz(counts, this.z), this.z, odd, counts, ionText(counts, this.z, odd), lost, how,
      Math.max(0, Math.min(1, score)), step);
    f.parent = parent;
    if (mech) { f.mech = mech[0]; f.ref = FRAG.REFS[mech[1]] ?? mech[1] ?? ''; }
    if (S) {
      f.S = new Set([...S].map(i => m.order[i]));
      f.marks = new Map([...marks].map(([i, v]) => [m.order[i], v]));
      f._idx = S; f._imarks = marks;
      f.atoms = new Set([...S].map(i => m.owner.get(m.order[i])));
    }
    return f;
  }

  add(f) { if (f) this.out.push(f); return f; }

  betaH(c, S, avoid = []) {
    const j = this.m.heavy(c).find(j => S.has(j) && !avoid.includes(j) && this.m.h(j));
    return j === undefined ? null : j;
  }

  run(steps = 2) {
    if (this.charged) this.chargeRemote();
    else if (this.mode === 'EI') { this.eiCleavages(); this.mclafferty(); } else this.esiCleavages();
    if (!this.charged) this.rda();
    this.groupLosses();
    let frontier = this.out.filter(f => f._idx);
    for (let step = 2; step <= steps; step++) {
      frontier = [...frontier].sort((a, b) => b.score - a.score).slice(0, 60);
      const fresh = [];
      for (const f of frontier) fresh.push(...this.nextLosses(f, step));
      this.out.push(...fresh);
      frontier = fresh.filter(f => f._idx);
    }
    return this.merged();
  }

  eiCleavages() {
    const m = this.m;
    const alpha = { N: 0.85, S: 0.7, O: 0.65 };
    for (let y = 0; y < m.n; y++) {
      if (!(m.el[y] in alpha) || m.aro[y] || m.q[y] > 0) continue;
      const acyl = m.nb[y].some(j => m.carbonylO(j) !== null);
      const Y = m.el[y], kind = { N: 'iminium', O: 'oxonium', S: 'sulfonium' }[Y];
      for (const ca of m.nb[y]) {
        if (!m.single(ca, y, true) || !m.saturated(ca)) continue;
        if (m.h(ca)) {
          this.add(this.make(this.all, marksOf([[ca, [1, 0, -1]]]), false, `loss of H• α to ${Y} (${kind} ion)`,
            (Y === 'N' ? 0.45 : 0.3) - (acyl ? 0.2 : 0), 1, null, null,
            [`Ionisation at ${Y} (${Y}•+). One electron of an α C–H bond pairs with the unpaired electron on ${Y} (C=${Y}+ forms), the other leaves with the hydrogen atom: [M−H]+.`, 'ei']));
        }
        for (const r of m.heavy(ca)) {
          if (r === y || !m.single(ca, r)) continue;
          const S = m.side(ca, r);
          if (!S) continue;
          const score = alpha[Y] - (acyl ? 0.3 : 0) + (m.side(r, ca).size >= 2 ? 0.05 : 0);
          this.add(this.make(S, marksOf([[ca, [1, 0, -1]]]), false, `α-cleavage at ${Y} (${kind} ion)`, score, 1, null, null,
            [`Ionisation removes an electron from a lone pair of ${Y} (${Y}•+). One electron of the Cα–${m.el[r]} bond pairs with the unpaired electron on ${Y}, making the C=${Y}+ double bond of the ${kind} ion; the other electron leaves with that side as a radical. The largest radical is lost most easily (Stevenson's rule).`, 'ei']));
        }
      }
    }
    for (let c = 0; c < m.n; c++) {
      const o = m.carbonylO(c);
      if (o === null) continue;
      if (m.hs[c]) {
        this.add(this.make(this.all, marksOf([[c, [1, 0, -1]]]), false, 'loss of H• from the aldehyde (acylium ion)', 0.55, 1, null, null,
          ["Ionisation at the C=O oxygen. The aldehyde C–H bond gives one electron to the oxygen's unpaired electron (C≡O+ forms) and the hydrogen atom leaves: [M−H]+, strong for aldehydes.", 'ei']));
      }
      for (const r of m.heavy(c)) {
        if (r === o || !m.single(c, r)) continue;
        const S = m.side(c, r);
        if (!S) continue;
        const score = m.el[r] === 'C' ? 0.85 : m.el[r] === 'O' ? (m.h(r) ? 0.45 : 0.75) : m.el[r] === 'N' ? 0.55 : 0.4;
        this.add(this.make(S, marksOf([[c, [1, 0, -1]]]), false, 'α-cleavage at C=O (acylium ion)', score, 1, null, null,
          [`Ionisation removes a lone-pair electron of the C=O oxygen (O•+). One electron of the C(=O)–${m.el[r]} bond pairs with it, making the C≡O+ triple bond of the acylium ion; the other electron leaves with the group as a radical.`, 'ei']));
        if (m.el[r] === 'C') {
          const S2 = m.side(r, c);
          this.add(this.make(S2, marksOf([[r, [1, 0, -1]]]), false, 'α-cleavage at C=O (charge on the other side)', 0.6 * m.cation(r, S2), 1, null, null,
            ['Inductive cleavage: the charged C=O oxygen pulls the C–C bond pair toward the carbonyl carbon; the acyl radical leaves and the alkyl / aryl cation remains.', 'ei']));
        }
      }
    }
    for (let b = 0; b < m.n; b++) {
      if (!m.saturated(b)) continue;
      const aro = m.nb[b].filter(j => m.aro[j]), allyl = m.nb[b].filter(j => m.vinyl(j));
      if (!aro.length && !allyl.length) continue;
      for (const r of m.heavy(b)) {
        if (aro.includes(r) || allyl.includes(r) || !m.single(b, r)) continue;
        const S = m.side(b, r);
        if (!S) continue;
        const f = this.make(S, marksOf([[b, [1, 0, -1]]]), false, aro.length ? 'benzylic cleavage' : 'allylic cleavage',
          aro.length ? 0.9 : 0.6, 1, null, null, [aro.length
            ? "Ionisation of the aromatic π system (ring•+). One electron of the benzylic C–R bond pairs with the ring's unpaired electron, the other leaves with R•. The benzyl cation formed rearranges to the seven-membered aromatic tropylium ion."
            : 'Ionisation of the C=C π bond. One electron of the allylic C–R bond pairs with the π radical, the other leaves with R•: the resonance-stabilised allyl cation remains.', 'ei']);
        if (f && sameCounts(f.counts, { C: 7, H: 7 })) f.how = 'benzylic cleavage (tropylium ion C7H7+)';
        this.add(f);
      }
    }
    const base = { I: 0.75, Br: 0.7, Cl: 0.45, F: 0.15 };
    for (let x = 0; x < m.n; x++) {
      for (const c of m.nb[x]) {
        if (m.el[c] !== 'C' || !m.single(c, x)) continue;
        const S = m.side(c, x);
        if (!S) continue;
        const X = m.el[x];
        const mech = [`Ionisation at ${X} (${X}•+). The C–${X} bond pair moves to ${X} (inductive cleavage, i): ${X} leaves as a radical and the carbon keeps the positive charge.`, 'ei'];
        if (X in base) {
          const s = base[X] * (0.5 + 0.5 * m.cation(c, S)) * (m.aro[c] ? 0.6 : 1);
          this.add(this.make(S, marksOf([[c, [1, 0, -1]]]), false, `inductive cleavage at C–${X}`, s, 1, null, null, mech));
        } else if ((X === 'O' || X === 'S') && !m.h(x) && m.saturated(c)) {
          const s = 0.3 * (0.5 + m.cation(c, S));
          this.add(this.make(S, marksOf([[c, [1, 0, -1]]]), false, `inductive cleavage at C–${X}`, s, 1, null, null, mech));
        }
      }
    }
    for (let a = 0; a < m.n; a++) {
      if (!m.saturated(a)) continue;
      for (const b of m.nb[a]) {
        if (b < a || !m.saturated(b) || !m.single(a, b)) continue;
        for (const [c, other] of [[a, b], [b, a]]) {
          const S = m.side(c, other);
          if (!S || S.size < 2) continue;
          const s = m.cation(c, S);
          if (s >= 0.35) {
            this.add(this.make(S, marksOf([[c, [1, 0, -1]]]), false, 'σ-cleavage of C–C (alkyl cation)', 0.1 + 0.4 * s, 1, null, null,
              ['σ-Ionisation of the C–C bond (no lone pair or π bond nearby). The one-electron bond breaks; the more substituted carbon keeps the charge (the most stable carbocation) and the other side leaves as a radical - branch points break most.', 'ei']));
          }
        }
      }
    }
    for (let x = 0; x < m.n; x++) {
      if (!['O', 'N'].includes(m.el[x]) || m.aro[x] || m.q[x]) continue;
      for (const c of m.heavy(x)) {
        const a = m.heavy(x).find(j => j !== c && (m.aro[j] || m.carbonylO(j) !== null));
        if (a === undefined || !m.saturated(c) || !m.single(c, x)) continue;
        const Sx = m.side(x, c), Sc = m.side(c, x);
        if (!Sx) continue;
        const bh = this.betaH(c, Sc, [x]);
        if (bh === null) continue;
        const ester = m.carbonylO(a) !== null;
        this.add(this.make(Sx, marksOf([[x, [1, 1, 0]]]), true, 'alkene loss, H to ' + (ester ? 'the ester O' : m.el[x] + ' (four-centre)'),
          ester ? 0.4 : 0.6, 1, null, null,
          ['A β-hydrogen of the alkyl group moves to the heteroatom as the C–' + m.el[x] + ' bond breaks (four-membered transition state); the alkene leaves and the ' +
            (ester ? 'acid' : 'phenol / amine') + ' radical cation remains (phenetole → phenol•+, m/z 94).', 'ei']));
        if (ester && [...Sc].filter(j => m.el[j] === 'C').length >= 2) {
          this.add(this.make(Sx, marksOf([[x, [1, 0, 1]]]), false, 'double hydrogen rearrangement (McLafferty+1)', 0.5, 1, null, null,
            ['Two hydrogens of the alkoxy group move to the ester oxygens (two consecutive H transfers) and an allyl / alkenyl radical leaves: the protonated acid [RC(OH)2]+ (ethyl acetate → m/z 61).', 'mclafferty']));
        }
      }
    }
  }

  mclafferty() {
    const m = this.m;
    const done = new Set();
    for (let c1 = 0; c1 < m.n; c1++) {
      const ys = m.nb[c1].filter(y => m.double(c1, y) && ['O', 'N', 'S', 'C'].includes(m.el[y]) && !m.aro[y]);
      const ortho = m.aro[c1] ? m.nb[c1].filter(y => m.aro[y] && m.el[y] === 'C') : [];
      if (m.el[c1] !== 'C' || (!ys.length && !ortho.length)) continue;
      for (const xa of m.heavy(c1)) {
        if (ys.includes(xa) || m.aro[xa] || !['C', 'O', 'N'].includes(m.el[xa]) || !m.single(c1, xa, true)) continue;
        if (m.el[xa] === 'C' && !m.saturated(xa)) continue;
        for (const cb of m.heavy(xa)) {
          if (cb === c1 || !m.saturated(cb) || !m.single(xa, cb)) continue;
          const S = m.side(xa, cb);
          if (!S) continue;
          const gs = m.nb[cb].filter(g => !S.has(g) && m.el[g] === 'C' && m.h(g) && !m.aro[g]);
          if (!gs.length) continue;
          for (const y of ys) {
            const key = `${y},${xa},${cb}`;
            if (done.has(key)) continue;
            done.add(key);
            const s = { O: 0.85, S: 0.6, N: 0.6, C: 0.35 }[m.el[y]] - (m.el[xa] === 'O' ? 0.15 : 0);
            const g = gs[0], Y = m.el[y];
            const enol = Y !== 'C' ? marksOf([[y, [1, 0, 1]], [xa, [0, 1, -1]]]) : marksOf([[y, [1, 1, 0]]]);
            this.add(this.make(S, enol, true, 'McLafferty rearrangement', s, 1, null, null,
              [`Ionisation at ${Y} (${Y}•+). Through a six-membered ring the unpaired electron on ${Y} takes a γ-hydrogen; the Cα–Cβ bond then breaks: one of its electrons pairs with the one left on Cγ (the neutral alkene Cβ=Cγ leaves), the other makes the C=C of the enol radical cation.`, 'mclafferty']));
            const S2 = m.side(cb, xa);
            this.add(this.make(S2, marksOf([[cb, [1, 0, -1]], [g, [0, 1, -1]]]), true, 'McLafferty rearrangement (charge on the alkene)', 0.3, 1, null, null,
              ['The same γ-hydrogen transfer and Cα–Cβ cleavage, but the charge stays on the alkene side (when its ionisation energy is the lower).', 'mclafferty']));
          }
          if (ortho.length && !ys.length && m.el[xa] === 'C') {
            this.add(this.make(S, marksOf([[ortho[0], [1, 1, 0]]]), true, 'McLafferty-type H shift to the ring', 0.55, 1, null, null,
              ['The ring radical cation takes a γ-hydrogen at its ortho carbon (six-membered transition state) and the benzylic Cα–Cβ bond breaks: an alkene leaves and the methylene-cyclohexadiene radical cation remains (propylbenzene → m/z 92).', 'mclafferty']));
          }
        }
      }
    }
  }

  rda() {
    const m = this.m;
    const seen = new Set();
    for (const ring of m.rings) {
      if (ring.length !== 6) continue;
      for (let k = 0; k < 6; k++) {
        const [p, a, b, q, u, w] = [-1, 0, 1, 2, 3, 4].map(d => ring[(k + d + 6) % 6]);
        if (!(m.double(a, b) || m.aromaticBond(a, b))) continue;
        if (!(m.single(p, a, true) && m.single(b, q, true) && m.single(q, u, true) && m.single(p, w, true))) continue;
        if (m.aromaticBond(u, w)) continue;
        const key = [pairKey(q, u), pairKey(p, w)].sort().join('|');
        if (seen.has(key)) continue;
        seen.add(key);
        const D = m.side(a, -1, [[q, u], [p, w]]);
        if (!D || D.has(u) || D.has(w)) continue;
        const E = new Set([...this.all].filter(i => !D.has(i)));
        const text = 'The six-membered ring opens at the two bonds across from its double (or fused aromatic) bond - a retro-Diels–Alder reaction: the three electron pairs of the ring shift together, giving a diene and an alkene (the dienophile). ';
        if (this.mode === 'EI') {
          const site = m.el[p] === 'C' ? p : q, other = site === p ? q : p;
          const mech = [text + 'In the radical cation the charge usually stays on the diene.', 'ei'];
          this.add(this.make(D, marksOf([[site, [1, 0, -1]], [other, [0, 1, -1]]]), true, 'retro-Diels–Alder (diene)', 0.7, 1, null, null, mech));
          this.add(this.make(E, marksOf([[w, [1, 0, -1]], [u, [0, 1, -1]]]), true, 'retro-Diels–Alder (dienophile)', 0.3, 1, null, null, mech));
        } else {
          const sign = this.z > 0 ? 1 : -1;
          for (const [piece, ends, s, part] of [[D, [p, q], 0.65, '1,3A-type diene'], [E, [w, u], 0.35, 'dienophile']]) {
            let marks;
            if (sign > 0) {
              const c = [...ends].sort((i, j) => {
                const ki = [m.el[i] !== 'C', m.carbonylO(i) === null], kj = [m.el[j] !== 'C', m.carbonylO(j) === null];
                return (ki[0] - kj[0]) || (ki[1] - kj[1]);
              })[0];
              marks = this.place(piece, null, c, 1, 0, -1);
            } else {
              const acid = this.acidSite(piece, ends);
              if (acid !== null) {
                marks = this.place(piece, null, acid, -1, 0, -1);
                for (const e of ends) marks = marks && marks.size ? this.place(piece, marks, e, 0, 0, -1) : null;
              } else marks = null;
            }
            this.add(this.make(piece, marks, false, `retro-Diels–Alder (${part})`, s, 1, null, null,
              [text + 'The protonated (or deprotonated) molecule splits the same way; for flavonoids the A-ring piece is the 1,3A ion and the B-ring piece the 1,3B ion.', 'rda']));
          }
        }
      }
    }
  }

  esiCleavages() {
    const m = this.m;
    const pos = this.z > 0;
    for (let c = 0; c < m.n; c++) {
      if (m.el[c] !== 'C') continue;
      for (const x of m.nb[c]) {
        if (!['N', 'O', 'S'].includes(m.el[x]) || !m.single(c, x)) continue;
        const Sc = m.side(c, x), Sx = m.side(x, c);
        if (!Sc) continue;
        const X = m.el[x];
        const acyl = m.carbonylO(c) !== null;
        const glyco = m.anomeric(c, x) && ['O', 'N', 'S'].includes(X);
        const pendant = m.heavy(x).length === 1;
        const xAcyl = m.heavy(x).some(j => j !== c && m.carbonylO(j) !== null);
        const xAryl = m.heavy(x).some(j => j !== c && m.aro[j]);
        const betaH = m.heavy(c).some(j => Sc.has(j) && j !== x && m.h(j)) && !m.aro[c];
        if (pos && !pendant) {
          let s, how;
          if (glyco) [s, how] = [0.55, 'glycosidic cleavage (sugar oxonium ion, B)'];
          else if (acyl) [s, how] = X === 'N' ? [0.7, 'amide cleavage (acylium, b-type)'] : [0.8, 'ester cleavage (acylium ion)'];
          else [s, how] = [0.9 * m.cation(c, Sc), `C–${X} cleavage (carbocation)`];
          if (s >= 0.15) {
            this.add(this.make(Sc, marksOf([[c, [1, 0, -1]]]), false, how, s, 1, null, null,
              [`Protonation of ${X} turns it into a good leaving group. The C–${X} bond pair goes to ${X}, which leaves as a neutral molecule, and the ` +
                (acyl ? 'acylium (b) ion' : glyco ? 'oxocarbenium (B) ion' : 'carbocation') + ' remains (charge-migration fragmentation).',
              acyl && X === 'N' ? 'peptide' : glyco ? 'glyco' : 'esi']));
          }
        }
        if (!pos && glyco && (xAryl || m.heavy(x).some(j => j !== c && m.vinyl(j)))) {
          const acid = this.acidSite(Sx, [x]);
          if (acid !== null) {
            const flavonol = m.heavy(x).some(j => j !== c && m.vinyl(j));
            this.add(this.make(Sx, marksOf([[acid, [-1, 0, -1]], [x, [0, 1, -1]]]), true, 'homolytic glycosidic cleavage (radical aglycone [Y0−H]•−)',
              flavonol ? 0.75 : 0.4, 1, null, null,
              ['The glycosidic O–C bond breaks homolytically: the aglycone keeps the oxygen as a radical next to the deprotonated phenol ([Y0−H]•−) and the sugar leaves as a radical. Strong for flavonol 3-O-glycosides, where the radical is conjugated with the 2,3-double bond and C-4 carbonyl.', 'radical']));
          }
        }
        if (!pos && acyl && X === 'O' && !pendant) {
          const acid = this.acidSite(Sc);
          const alpha = this.betaH(c, Sc, [x]);
          if (acid !== null && alpha !== null) {
            let marks = this.place(Sc, marksOf([[acid, [-1, 0, -1]]]), c, 0, 0, -1);
            marks = marks && marks.size ? this.place(Sc, marks, alpha, 0, 0, -1) : null;
            this.add(this.make(Sc, marks, false, 'ester cleavage (acyl side as the anion, ketene)', 0.45, 1, null, null,
              ['The deprotonated phenol / acid on the acyl side stays charged while the ester O–C(=O) bond breaks with an α-hydrogen moving to the alkoxy oxygen: the alcohol leaves as a neutral molecule and the acyl part remains as a ketene anion (chlorogenic acids: caffeoyl m/z 161).', 'chlorogenic']));
          }
        }
        if (!betaH || (pendant && !glyco)) continue;
        if (pos) {
          let s, how;
          if (glyco) [s, how] = [0.9, 'glycosidic cleavage (aglycone or sugar loss, Y ion)'];
          else if (acyl && X === 'N') [s, how] = [0.65, 'amide cleavage, charge on the amine (y-type)'];
          else if (acyl) [s, how] = [0.35, 'ester cleavage, charge on the alcohol (ketene lost)'];
          else if (xAcyl && X === 'O') [s, how] = [0.5, 'ester cleavage: [acid+H]+, alkene lost'];
          else if (xAcyl) [s, how] = [0.35, 'C–N cleavage: [amide+H]+, alkene lost'];
          else if (X === 'N') [s, how] = [0.45, 'C–N cleavage: amine, alkene lost'];
          else [s, how] = [0.3, `C–${X} cleavage, alkene lost`];
          this.add(this.make(Sx, marksOf([[x, [1, 0, 1]]]), false, how, s, 1, null, null,
            [`Proton on ${X}. As the C–${X} bond breaks, a hydrogen next to the carbon goes with the electrons into a C=C (or C=C=O) bond and is handed over to ${X} (1,2-elimination, charge retention on the heteroatom side): the neutral ` +
              (glyco ? 'anhydro sugar leaves and the aglycone (Y) ion remains.' : acyl ? 'ketene leaves.' : 'alkene leaves.'),
            glyco ? 'glyco' : acyl && X === 'N' ? 'peptide' : 'esi']));
        } else {
          let s, how;
          if (glyco) [s, how] = [0.85, 'glycosidic cleavage (aglycone or sugar loss, Y⁻ ion)'];
          else if (xAcyl && X === 'O') [s, how] = [0.7, 'ester cleavage (carboxylate ion)'];
          else if (acyl && X === 'O') [s, how] = [0.6, 'ester cleavage (alkoxide ion, ketene lost)'];
          else if (xAryl) [s, how] = [0.5, 'C–O cleavage (phenoxide ion)'];
          else if (acyl) [s, how] = [0.3, 'amide cleavage (amide anion)'];
          else [s, how] = [0.2, `C–${X} cleavage (anion)`];
          this.add(this.make(Sx, marksOf([[x, [-1, 0, -1]]]), false, how, s, 1, null, null,
            [`The negative charge (the deprotonated acidic OH) takes a hydrogen next to the carbon as the C–${X} bond pair goes to ${X}: ${X} leaves with the charge and the carbon side becomes a neutral ` +
              (glyco ? 'anhydro sugar.' : acyl ? 'ketene.' : 'alkene.'),
            xAcyl || acyl ? 'chlorogenic' : glyco ? 'glyco' : 'esi']));
        }
      }
    }
  }

  chargeRemote() {
    const m = this.m;
    for (let i = 0; i < m.n; i++) {
      for (const j of m.heavy(i)) {
        if (j < i || !m.single(i, j)) continue;
        for (const [a, b] of [[i, j], [j, i]]) {
          const S = m.side(a, b);
          if (!S || ![...S].some(k => m.q[k])) continue;
          if (!m.heavy(b).some(k => !S.has(k) && m.h(k)) && !m.h(b)) continue;
          this.add(this.make(S, new Map(), this.odd, 'cleavage, neutral molecule lost (charge-remote)', 0.45, 1, null, null,
            ['The charge stays where it is drawn: a hydrogen moves to the charged side as the bond breaks and the other side leaves as a neutral molecule (charge-remote fragmentation).', 'esi']));
        }
      }
    }
  }

  groupMatches() {
    if (!this.gm) this.gm = FRAG.GROUP_LOSSES.map(row => [row, this.m.matches(row[1])]);
    return this.gm;
  }

  groupLosses() {
    const m = this.m;
    for (const [row, matches] of this.groupMatches()) {
      const [name, , rem, att, lh, radical, scores, what] = row;
      const s = scores[this.mode];
      if (!s) continue;
      for (const match of matches) {
        const gone = new Set(rem.map(k => match[k]));
        for (const g of [...gone]) for (const j of m.nb[g]) if (m.el[j] === 'H') gone.add(j);
        const S = new Set([...this.all].filter(i => !gone.has(i)));
        const attach = match[att];
        if (!S.has(attach)) continue;
        const reach = m.side(attach, -1, [...gone].map(g => [attach, g]));
        if (!reach || !sameSet(reach, S)) continue;
        this.add(this.loss(S, null, attach, radical, name, lh, s, `loss of ${name} (${what})`, null, 1, this.lossMech(row)));
      }
    }
    const everything = new Set([...this.all].map(i => m.owner.get(m.order[i])));
    for (const [name, smarts, lost, radical, scores, what] of FRAG.FORMULA_LOSSES) {
      const s = scores[this.mode];
      if (!s) continue;
      if (m.matches(smarts).length) {
        const ref = what.includes('phospho') ? 'lipid' : what.includes('chromone') ? 'rda' : what.includes('glycoside') ? 'glyco' : 'ei';
        this.add(this.formulaOnly(lost, radical, `loss of ${name} (${what})`, s, everything, null, 1, [FRAG.FORMULA_TEXT[what] || '', ref]));
      }
    }
    for (const [counts, odd, [smarts, most], scores, what, text, ref] of FRAG.DIAGNOSTIC) {
      const s = scores[this.mode];
      if (!s || m.matches(smarts).length < most) continue;
      this.add(this.fixedIon(counts, odd, what, s, everything, [text, ref]));
    }
    if (this.z > 0 && !this.odd) this.immonium();
  }

  immonium() {
    const m = this.m;
    for (let ca = 0; ca < m.n; ca++) {
      if (!m.saturated(ca)) continue;
      const n = m.heavy(ca).find(j => m.el[j] === 'N' && !m.aro[j] && m.single(ca, j, true));
      const co = m.heavy(ca).find(j => m.carbonylO(j) !== null && m.single(ca, j));
      if (n === undefined || co === undefined) continue;
      const prev = m.heavy(n).filter(j => j !== ca && m.carbonylO(j) !== null && m.single(n, j));
      const S = m.side(ca, co, prev.map(j => [n, j]));
      if (!S || prev.some(j => S.has(j))) continue;
      const aromatic = [...S].some(j => m.aro[j]);
      this.add(this.make(S, marksOf([[ca, [1, 0, -1]]]), false, 'immonium ion of the residue', aromatic ? 0.7 : 0.45, 1, null, null,
        ["An internal fragment (b → a ion of one residue): the amide bond before the residue and its Cα–C(=O) bond break, CO leaves, and the iminium ion H2N+=CH–R remains. It tells which amino acid is present (Phe 120.08, Tyr 136.08, Trp 159.09, His 110.07, Leu / Ile 86.10, Pro 70.07).", 'peptide']));
    }
  }

  lossMech(row) {
    const what = row[7];
    const texts = FRAG.LOSS_TEXT[what] || ['', null];
    const text = this.mode === 'EI' && texts[1] ? texts[1] : texts[0];
    const ref = what.includes('silyl') ? 'tms' : what.includes('phosph') ? 'lipid' : this.mode !== 'EI' ? 'esi' : 'ei';
    return [text, ref];
  }

  loss(S, base, attach, radical, name, lh, s, how, parent, step, mech = null) {
    const m = this.m;
    const start = parent || null;
    const odd = (start ? start.odd : this.odd) !== radical;
    const lost = {};
    const gone = start ? [...start._idx].filter(i => !S.has(i)) : [...this.all].filter(i => !S.has(i));
    for (const g of gone) if (m.el[g] !== 'H') lost[m.key[g]] = (lost[m.key[g]] || 0) + 1;
    lost.H = lh;
    const target = sub(start ? start.counts : this.pre, lost);
    const hc = m.hcap(S);
    if (hc === null || Object.values(target).some(v => v < 0)) return null;
    const have = hc + sumMarks(base, 2, S);
    let dh = (target.H || 0) - have, dq, rad;
    if (!start) {
      dq = this.z; for (const i of S) dq -= m.q[i];
      rad = odd ? 1 : 0;
      if (dq < 0 && (lh || radical)) {
        const acid = this.acidSite(S, [attach]);
        if (acid !== null) { base = marksOf([[acid, [-1, 0, -1]]]); dq = 0; dh += 1; }
      } else if (dq > 0 && this.mode !== 'EI' && (radical || !lh)) {
        const site = this.baseSite(S, [attach]);
        if (site === null) return null;
        base = marksOf([[site, [1, 0, 1]]]); dq = 0; dh -= 1;
      }
    } else {
      if (this.z > 0 && lh && !radical) {
        let proton = null;
        for (const [k, v] of base) if (S.has(k) && v[0] > 0 && v[2] > 0 && 'NOS'.includes(m.el[k])) { proton = k; break; }
        if (proton !== null) {
          const v = base.get(proton);
          base = new Map(base);
          base.set(proton, [v[0] - 1, v[1], v[2] - 1]);
          dh += 1;
        }
      }
      dq = this.z; for (const i of S) dq -= m.q[i];
      dq -= sumMarks(base, 0, S);
      rad = (odd ? 1 : 0) - sumMarks(base, 1, S);
      if (rad < 0) {
        base = new Map([...base].map(([k, v]) => [k, [...v]]));
        for (const v of base.values()) if (v[1] > 0 && rad < 0) { v[1] -= 1; rad += 1; }
      }
    }
    const marks = this.place(S, base, attach, dq, rad, dh);
    if (!marks) return null;
    return this.make(S, marks, odd, how, s, step, name + (!radical || name.endsWith('•') ? '' : '•'), parent, mech);
  }

  nextLosses(f, step) {
    const m = this.m;
    const out = [];
    const S0 = f._idx;
    const cutOut = [];
    for (const x of S0) for (const y of m.nb[x]) if (!S0.has(y)) cutOut.push([x, y]);
    for (const [row, matches] of this.groupMatches()) {
      const [name, , rem, att, lh, radical, scores] = row;
      const s = scores[this.mode];
      if (!s || (radical && !f.odd)) continue;
      for (const match of matches) {
        const gone = new Set(rem.map(k => match[k]));
        for (const g of [...gone]) for (const j of m.nb[g]) if (m.el[j] === 'H') gone.add(j);
        const attach = match[att];
        if (![...gone].every(g => S0.has(g)) || !S0.has(attach) || [...gone].some(g => f._imarks.has(g))) continue;
        let charged = false;
        for (const g of [...gone, attach]) for (const j of [...m.nb[g], g]) if (f._imarks.has(j) && f._imarks.get(j)[0]) charged = true;
        if (charged) continue;
        const S = new Set([...S0].filter(i => !gone.has(i)));
        const reach = m.side(attach, -1, [...[...gone].map(g => [attach, g]), ...cutOut]);
        if (!reach || !sameSet(reach, S)) continue;
        const g = this.loss(S, f._imarks, attach, radical, name, lh, f.score * s * 0.85, `${f.how}, then −${name}`, f, step, this.lossMech(row));
        if (g) out.push(g);
      }
    }
    for (const [c, v] of f._imarks) {                       // an acylium ion loses CO
      const o = m.carbonylO(c);
      if (v[0] <= 0 || o === null || !S0.has(o)) continue;
      const rest = m.heavy(c).filter(j => S0.has(j) && j !== o);
      if (rest.length !== 1 || m.el[rest[0]] !== 'C') continue;
      const S = new Set([...S0].filter(i => i !== c && i !== o));
      const base = new Map([...f._imarks].filter(([k]) => k !== c));
      const hc = m.hcap(S);
      if (hc === null) continue;
      const target = sub(f.counts, { C: 1, O: 1 });
      const dh = (target.H || 0) - hc - sumMarks(base, 2, S);
      const marks = this.place(S, base, rest[0], v[0], 0, dh);
      const g = this.make(S, marks, f.odd, `${f.how}, then −CO`, f.score * 0.75, step, 'CO', f,
        ['The acylium ion R–C≡O+ expels carbon monoxide: the R–C bond pair goes with the CO carbon and the carbocation R+ remains (benzoyl 105 → phenyl 77; peptide b → a ions).', 'ei']);
      if (g) out.push(g);
    }
    for (const [x, v] of f._imarks) {                       // a carboxylate ion loses CO2
      let c = null;
      if (v[0] < 0 && m.el[x] === 'O') { const k = m.heavy(x).find(j => S0.has(j)); c = k === undefined ? null : k; }
      const o = c !== null ? m.carbonylO(c) : null;
      if (o === null || !S0.has(o)) continue;
      const rest = m.heavy(c).filter(j => S0.has(j) && j !== o && j !== x);
      if (rest.length !== 1) continue;
      const S = new Set([...S0].filter(i => i !== c && i !== o && i !== x));
      const base = new Map([...f._imarks].filter(([k]) => k !== x));
      const hc = m.hcap(S);
      if (hc === null) continue;
      const target = sub(f.counts, { C: 1, O: 2 });
      const dh = (target.H || 0) - hc - sumMarks(base, 2, S);
      const marks = this.place(S, base, rest[0], v[0], 0, dh);
      const g = this.make(S, marks, f.odd, `${f.how}, then −CO2`, f.score * 0.8, step, 'CO2', f,
        ['The carboxylate COO− expels CO2: the C–C bond pair stays on the carbon, which keeps the negative charge (caffeate 179 → 135).', 'chlorogenic']);
      if (g) out.push(g);
    }
    if (f.z > 0 && !f.odd) {                                // an onium ion loses an alkene, its H moving to N / O
      for (const y of sorted(S0)) {
        if (!['N', 'O'].includes(m.el[y]) || m.aro[y]) continue;
        for (const c of m.heavy(y)) {
          if (!S0.has(c) || !m.saturated(c) || !m.single(c, y)) continue;
          const R = m.side(c, y);
          if (!R || ![...R].every(k => S0.has(k)) || [...R].some(k => f._imarks.has(k))) continue;
          if (!m.heavy(c).some(j => R.has(j) && m.h(j))) continue;
          const S = new Set([...S0].filter(i => !R.has(i)));
          const hc = m.hcap(S);
          if (hc === null) continue;
          const kept = clean(m.counts(S, hc + sumMarks(f._imarks, 2, S)));
          const lost = clean(sub(f.counts, kept));
          const g = this.make(S, new Map([...f._imarks].filter(([k]) => S.has(k))), f.odd,
            `${f.how}, then −${formulaText(lost)} (alkene; H to ${m.el[y]})`, f.score * 0.6, step, null, f,
            [`The onium ion loses an alkene: a hydrogen of the alkyl group moves to ${m.el[y]} (four-membered transition state) as the ${m.el[y]}–C bond breaks (triethylamine 86 → 58).`, 'ei']);
          if (g) out.push(g);
        }
      }
    }
    if ((sameCounts(f.counts, { C: 7, H: 7 }) || sameCounts(f.counts, { C: 6, H: 5 })) && f.z > 0 && !f.odd) {
      const g = this.formulaOnly({ C: 2, H: 2 }, false, `${f.how}, then −C2H2`, f.score * 0.5, f.atoms, f, step,
        ['The tropylium (or phenyl) ion breaks its ring and loses acetylene: 91 → 65, 77 → 51.', 'ei']);
      if (g) out.push(g);
    }
    return out;
  }

  merged() {
    const best = new Map();
    for (const f of this.out) {
      const key = `${f.formula}|${f.z}|${f.odd}`;
      const old = best.get(key);
      if (!old) { best.set(key, f); continue; }
      if (f.score > old.score) {
        f.also = [old.how, ...old.also.filter(h => h !== f.how)];
        best.set(key, f);
      } else if (f.how !== old.how && !old.also.includes(f.how)) old.also.push(f.how);
    }
    const out = [...best.values()];
    for (const f of out) f.score = Math.min(1, f.score + 0.05 * Math.min(f.also.length, 3));
    out.sort((a, b) => b.mz - a.mz);
    return out.slice(0, 500);
  }
}

/** {pre, frags} - the precursor ion and its fragments (highest m/z first) of the structure with these atoms (all
 * when null), as the desktop's fragments.predict; {pre: null, frags: []} when the structure cannot be read. */
export function predict(RD, doc, atomIds = null, mode = 'ESI+', steps = 2) {
  const m = new Mol(RD, doc, atomIds);
  m.RD = RD;
  if (!m.mol || !m.n) { m.dispose(); return { pre: null, frags: [] }; }
  try {
    if ((mode === 'ESI+NH4' || mode === 'ESI+Na') && !m.q.reduce((s, v) => s + v, 0)) return adduct(m, mode, steps);
    const e = new Engine(m, mode);
    if ((e.pre.H || 0) < 0) return { pre: null, frags: [] };
    const atoms = new Set(m.order.map(a => m.owner.get(a)));
    const pre = new Fragment(e.mz(e.pre, e.z), e.z, e.odd, e.pre, ionText(e.pre, e.z, e.odd), '', `precursor ion ${e.label}`, 1, 0);
    pre.atoms = atoms; pre.label = e.label; pre.drawnIon = e.charged;
    const frags = e.run(steps);
    for (const f of frags) f.lostAtoms = f.S ? new Set([...atoms].filter(a => !f.atoms.has(a))) : new Set();
    return { pre, frags };
  } finally { m.dispose(); }
}

function adduct(m, mode, steps) {
  const e = new Engine(m, 'ESI+');
  const frags = e.run(steps);
  const atoms = new Set(m.order.map(a => m.owner.get(a)));
  const plus = mode === 'ESI+NH4' ? { N: 1, H: 4 } : { Na: 1 };
  const counts = { ...e.M };
  for (const [k, v] of Object.entries(plus)) counts[k] = (counts[k] || 0) + v;
  const label = mode === 'ESI+NH4' ? '[M+NH4]+' : '[M+Na]+';
  const pre = new Fragment(e.mz(counts, 1), 1, false, counts, ionText(counts, 1, false), '', `precursor ion ${label}`, 1, 0);
  pre.atoms = atoms; pre.label = label; pre.drawnIon = false;
  let out = [];
  if (mode === 'ESI+NH4') {
    const base = e.baseSite(e.all);
    const mh = new Fragment(e.mz(e.pre, 1), 1, false, { ...e.pre }, ionText(e.pre, 1, false), 'NH3', 'loss of NH3 from the ammonium adduct ([M+H]+)', 0.9, 1);
    mh.atoms = atoms;
    if (base !== null) { mh.S = new Set(m.order); mh.marks = new Map([[m.order[base], [1, 0, 1]]]); }
    mh.mech = 'The ammonium ion is held by hydrogen bonds; on collision it hands its proton to the molecule and ammonia leaves: [M+H]+, which then breaks up as a protonated molecule.';
    mh.ref = FRAG.REFS.esi;
    out.push(mh);
    for (const f of frags) {
      f.lost = 'NH3 + ' + f.lost; f.step += 1; f.score *= 0.9;
      if (!f.parent) f.parent = mh;
      out.push(f);
    }
  } else {
    for (const f of frags) {
      if (f.step !== 1 || !f.marks || !f.marks.size || f.odd) continue;
      if (![...f.marks.values()].some(v => v[0] > 0 && v[2] > 0)) continue;
      let c = { ...f.counts };
      c.H = (c.H || 0) - 1;
      c.Na = (c.Na || 0) + 1;
      c = clean(c);
      const lost = clean(sub(counts, c));
      if (Object.values(lost).some(v => v < 0) || !Object.keys(lost).length) continue;
      const g = new Fragment(e.mz(c, 1), 1, false, c, ionText(c, 1, false), formulaText(lost),
        f.how.split('[acid+H]+').join('[acid+Na]+').split('Y ion').join('Y ion, sodiated') + ' (Na+)', f.score * 0.8, 1);
      g.atoms = f.atoms;
      g.mech = 'With Na+ in place of the proton the same bond breaks, the sodium staying on the oxygen-rich piece (charge-remote cleavage of the sodium adduct). As for [M+H]+: ' + f.mech;
      g.ref = FRAG.REFS.sodium;
      out.push(g);
    }
    const everything = new Set(atoms);
    for (const [name, smarts, lost, score, text] of FRAG.SODIUM_LOSSES) {
      if (!m.matches(smarts).length) continue;
      let c = sub(counts, lost);
      if (Object.values(c).some(v => v < 0)) continue;
      c = clean(c);
      const g = new Fragment(e.mz(c, 1), 1, false, c, ionText(c, 1, false), name, `loss of ${name} (sodium adduct)`, score, 1);
      g.atoms = everything; g.mech = text; g.ref = FRAG.REFS.sodium;
      out.push(g);
    }
  }
  const best = new Map();
  for (const f of out) {
    const key = `${f.formula}|${f.odd}`;
    if (!best.has(key) || f.score > best.get(key).score) best.set(key, f);
  }
  out = [...best.values()].sort((a, b) => b.mz - a.mz);
  for (const f of out) f.lostAtoms = f.S ? new Set([...atoms].filter(a => !f.atoms.has(a))) : new Set();
  return { pre, frags: out };
}
