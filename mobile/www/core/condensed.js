// condensed.js - a typed atom label read as a condensed formula, the same as the desktop's condensed.py: SO3K, CH2CH3,
// COCH3, N(CH3)2, CH(CH3)2, PO3H2, SO2NH2, CONHMe, (CH2)3CH3, CH=CH2, CO2-, NMe3+ ... The first atom is the one on the
// structure. The molecule is built here as a graph and written as a SMILES rooted at that first atom; RDKit.js (setRD)
// reads the nicknames used inside (Ph, Me ...) and checks the result.
import { NICKNAMES } from './data.js';
import { parseBlock } from './doc.js';

const ELECTRONS = { H: 1, B: 3, C: 4, N: 5, O: 6, F: 7, Si: 4, P: 5, S: 6, Cl: 7, Se: 6, Br: 7, I: 7, As: 5, Te: 6, Ge: 4,
  Sn: 4, Li: 1, Na: 1, K: 1, Rb: 1, Cs: 1, Mg: 2, Zn: 2, Al: 3, Ca: 2 };
const HYPER = { P: 5, S: 6, Se: 6, As: 5, Te: 6, I: 7, Cl: 7, Br: 5 };
const SALT_METALS = ['Li', 'Na', 'K', 'Rb', 'Cs'];
const ALIASES = { C6H5: 'Ph', C6H11: 'Cy' };
const ORGANIC = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']);

let RD = null;
const cache = new Map();
export function setRD(rd) { RD = rd; cache.clear(); }

function valence(sym, charge = 0) {
  const e = ELECTRONS[sym] - charge;
  if (charge === 0 && HYPER[sym] !== undefined) return HYPER[sym];
  return e > 0 ? Math.min(e, 8 - e) : 0;
}

const isLower = c => !!c && c >= 'a' && c <= 'z';
const isUpper = c => !!c && c >= 'A' && c <= 'Z';
const isDigit = c => !!c && c >= '0' && c <= '9';

function plainFormula(k) {
  const syms = k.match(/[A-Z][a-z]?/g) || [];
  return /^(?:[A-Z][a-z]?)+$/.test(k) && syms.every(t => ELECTRONS[t] !== undefined);
}

function tokenize(text) {
  const names = Object.keys(NICKNAMES).filter(k => /^[A-Za-z]+$/.test(k) && k.length > 1 && !plainFormula(k))
    .sort((a, b) => b.length - a.length);
  const out = [];
  let i = 0;
  while (i < text.length) {
    const ch = text[i];
    if (/\s/.test(ch)) { i++; continue; }
    if (ch === '(' || ch === '[') { out.push(['open', null, 1]); i++; continue; }
    if (ch === ')' || ch === ']') { out.push(['close', null, 1]); i++; }
    else if (ch === '=' || ch === '#' || ch === '≡') { out.push(['bond', ch === '=' ? 2 : 3, 1]); i++; continue; }
    else if (ch === '+' || ch === '-') {
      const signs = /^[+-]+/.exec(text.slice(i))[0];
      let n = signs.length;
      if (out.length && out[out.length - 1][0] === 'num') n = out.pop()[1];
      out.push(['charge', signs[0] === '+' ? n : -n, 1]); i += signs.length; continue;
    } else {
      const nick = names.find(n => text.startsWith(n, i) && !isLower(text[i + n.length]));
      if (nick) { out.push(['nick', nick, 1]); i += nick.length; }
      else if (isUpper(ch) && isLower(text[i + 1]) && ELECTRONS[text.slice(i, i + 2)] !== undefined) { out.push(['el', text.slice(i, i + 2), 1]); i += 2; }
      else if (isUpper(ch) && ELECTRONS[ch] !== undefined) { out.push(['el', ch, 1]); i++; }
      else if (isDigit(ch)) { const m = /^\d+/.exec(text.slice(i))[0]; out.push(['num', parseInt(m, 10), 1]); i += m.length; continue; }
      else return null;
    }
    const m = /^\d+/.exec(text.slice(i));
    if (m) { out[out.length - 1][2] = parseInt(m[0], 10); i += m[0].length; }
  }
  if (out.some(t => t[0] === 'num')) return null;
  return out;
}

/** {atoms: [{el, charge, aromatic}], bonds: [[a, b, order]]} of a nickname's SMILES (first atom first), via RDKit.js. */
function nickGraph(label) {
  const m = RD && RD.get_mol(NICKNAMES[label]);
  if (!m) throw new Error('nickname');
  try {                                                    // the kekulé form, atoms in the SMILES's order
    const g = parseBlock(m.get_molblock());
    return { atoms: g.atoms.map(a => ({ el: a.el, charge: a.charge })), bonds: g.bonds.map(bd => [bd.a, bd.b, bd.order]) };
  } finally { m.delete(); }
}

class Builder {
  constructor() { this.atoms = []; this.bonds = []; this.free = []; this.hs = []; }
  add(sym, charge = 0) {
    const a = this.atoms.length;
    this.atoms.push({ el: sym, charge });
    this.free.push(valence(sym, charge) - (a === 0 ? 1 : 0));   // the first atom keeps a bond for the structure
    this.hs.push(0);
    return a;
  }
  bond(a, b, order = 1) { this.bonds.push([a, b, order]); this.free[a] -= order; this.free[b] -= order; }
  neighbors(a) { return this.bonds.filter(b => b[0] === a || b[1] === a).map(b => (b[0] === a ? b[1] : b[0])); }
  bondOf(a, b) { return this.bonds.find(x => (x[0] === a && x[1] === b) || (x[0] === b && x[1] === a)); }
}

function parse(text) {
  const toks = tokenize(text);
  if (!toks || !toks.length) return null;
  if (toks[toks.length - 1][0] === 'charge') {               // NMe3+, SO3-: the last atom written, outside brackets
    const charge = toks.pop()[1];
    let last = null, depth = 0;
    for (let k = toks.length - 1; k >= 0; k--) {
      const t = toks[k];
      depth += t[0] === 'close' ? 1 : t[0] === 'open' ? -1 : 0;
      if (depth === 0 && t[0] === 'el' && t[1] !== 'H') { last = t; break; }
    }
    if (!last) return null;
    last.push(charge);
  }
  const b = new Builder();
  const scopes = [[]], entries = [null], pendantOf = new Map();
  let pendingOrder = 1, lastHeavy = null;

  const latestFree = (scope, need = 1) => {
    for (let k = scope.length - 1; k >= 0; k--) {
      const a = scope[k];
      if (b.free[a] >= need) {
        const host = pendantOf.get(a);
        if (host !== undefined) {
          let oxo = 0;
          for (const [p, h] of pendantOf) if (h === host) oxo += b.free[p];
          if (b.free[host] - oxo >= need) return host;
        }
        return a;
      }
    }
    return null;
  };
  const attach = (nw, scope, order = 1, onto = null) => {
    let target = onto !== null ? onto : latestFree(scope, order);
    if (target === null && !scope.length && entries[entries.length - 1] !== null && b.free[entries[entries.length - 1]] >= order) target = entries[entries.length - 1];
    if (target === null && scopes.length > 1) {
      for (let s = scopes.length - 2; s >= 0; s--) { target = latestFree(scopes[s], order); if (target !== null) break; }
    }
    if (target === null) { if (nw !== 0) throw new Error('no free bond'); return; }
    b.bond(target, nw, order);
  };
  const run = tokens => {
    let k = 0;
    while (k < tokens.length) {
      const [kind, value, count] = tokens[k];
      const tokenCharge = tokens[k].length > 3 ? tokens[k][3] : 0;
      if (kind === 'open') {
        let depth = 1, j = k + 1;
        while (j < tokens.length && depth) { depth += tokens[j][0] === 'open' ? 1 : tokens[j][0] === 'close' ? -1 : 0; j++; }
        if (depth) throw new Error('unclosed bracket');
        const inner = tokens.slice(k + 1, j - 1), times = tokens[j - 1][2];
        const branch = latestFree(scopes[scopes.length - 1]);
        if (times === 1 && inner.length === 1 && inner[0][0] === 'el' && (inner[0][1] === 'O' || inner[0][1] === 'S') && inner[0][2] === 1
            && branch !== null && b.free[branch] >= 2) {      // C(O)CH3, P(O)(OH)2: (O) is =O
          const a = b.add(inner[0][1]);
          b.bond(branch, a, 2);
          scopes[scopes.length - 1].push(a);
          k = j; continue;
        }
        let previous = [];
        for (let r = 0; r < times; r++) {
          const chain = latestFree(previous);
          entries.push(chain !== null ? chain : branch);
          scopes.push([]);
          run(inner);
          const done = scopes.pop();
          entries.pop();
          scopes[scopes.length - 1].push(...done);
          previous = done;
        }
        k = j; continue;
      }
      if (kind === 'close') throw new Error('stray bracket');
      if (kind === 'bond') pendingOrder = value;
      else if (kind === 'charge') {
        if (lastHeavy === null) throw new Error('charge on nothing');
        const at = b.atoms[lastHeavy];
        const used = valence(at.el) - b.free[lastHeavy];
        at.charge = value;
        b.free[lastHeavy] = valence(at.el, value) - used;
      } else if (kind === 'el' && value === 'H') {
        for (let r = 0; r < count; r++) {
          let a = scopes[scopes.length - 1].length ? latestFree(scopes[scopes.length - 1]) : null;
          if (a === null) a = latestFree(scopes.flat());
          if (a === null) throw new Error('no room for H');
          b.free[a] -= 1; b.hs[a] += 1;
        }
      } else if (kind === 'el') {
        const scope = scopes[scopes.length - 1];
        const onto = count > 1 ? latestFree(scope) : null;
        for (let r = 0; r < count; r++) {
          const a = b.add(value, r === count - 1 ? tokenCharge : 0);
          const here = onto !== null && b.free[onto] >= pendingOrder ? onto : null;
          attach(a, scope, pendingOrder, here);
          pendingOrder = 1;
          if (here !== null && (value === 'O' || value === 'S')) pendantOf.set(a, here);
          const host = b.neighbors(a);
          const next = tokens[k + 1];
          const carbonyl = count === 1 && (value === 'O' || value === 'S') && host.length && ['C', 'P', 'S'].includes(b.atoms[host[0]].el)
            && b.free[host[0]] >= 2 && next && ['el', 'nick', 'open'].includes(next[0]) && !(next[0] === 'el' && next[1] === 'H');
          if (carbonyl) { b.bondOf(host[0], a)[2] = 2; b.free[host[0]] -= 1; b.free[a] -= 1; }
          scope.push(a);
          lastHeavy = a;
        }
      } else if (kind === 'nick') {
        const sub = nickGraph(value);
        for (let r = 0; r < count; r++) {
          const base = b.atoms.length;
          for (const at of sub.atoms) { b.atoms.push({ el: at.el, charge: at.charge, nick: true }); b.free.push(0); b.hs.push(-1); }
          for (const [p, q, o] of sub.bonds) b.bonds.push([base + p, base + q, o]);
          b.free[base] = 1;
          attach(base, scopes[scopes.length - 1], pendingOrder);
          pendingOrder = 1;
          b.free[base] = 0;
          scopes[scopes.length - 1].push(base);
        }
      }
      k++;
    }
  };
  try { run(toks); } catch { return null; }
  const n = b.atoms.length;
  if (n < 2 || Math.min(...b.free) < 0) return null;
  for (const order of [2, 3]) {                             // leftover free bonds pair up: C=O, S=O, C#N, N=C=O
    let changed = true;
    while (changed) {
      changed = false;
      for (const bd of b.bonds) {
        const [p, q] = bd, sp = b.atoms[p].el, sq = b.atoms[q].el;
        const hyperToC = (HYPER[sp] !== undefined && sq === 'C') || (HYPER[sq] !== undefined && sp === 'C');
        if (bd[2] < order && b.free[p] > 0 && b.free[q] > 0 && b.hs[p] >= 0 && b.hs[q] >= 0 && !hyperToC) {
          bd[2] += 1; b.free[p] -= 1; b.free[q] -= 1; changed = true;
        }
      }
    }
  }
  for (const bd of [...b.bonds]) {                          // SO3K: a potassium sulfonate, O- K+
    const [p, q, o] = bd;
    for (const [metal, other] of [[p, q], [q, p]]) {
      if (SALT_METALS.includes(b.atoms[metal].el) && ['O', 'S', 'N'].includes(b.atoms[other].el) && o === 1) {
        b.bonds.splice(b.bonds.indexOf(bd), 1);
        b.atoms[metal].charge = 1; b.atoms[other].charge -= 1;
        b.free[metal] = 0; b.free[other] = Math.max(0, b.free[other]);
        break;
      }
    }
  }
  // hydrogens: as written; an atom left short takes the rest, except P / S / I ... which take their lower valence
  const hcount = b.atoms.map((at, i) => (b.hs[i] < 0 ? null : b.hs[i] + (b.free[i] > 0 && HYPER[at.el] === undefined ? b.free[i] : 0)));
  const smi = writeSmiles(b.atoms, b.bonds, hcount);
  if (RD) {
    const m = RD.get_mol(smi);
    if (!m) return null;
    m.delete();
  }
  return smi;
}

/** A SMILES rooted at atom 0 (depth-first, ring closures numbered), every atom in brackets with its H count. */
function writeSmiles(atoms, bonds, hcount) {
  const adj = atoms.map(() => []);
  bonds.forEach(([p, q, o], k) => { adj[p].push([q, o, k]); adj[q].push([p, o, k]); });
  const seen = new Set(), used = new Set(), closures = new Map();
  let ring = 1;
  const sym = o => (o === 2 ? '=' : o === 3 ? '#' : '');
  // first pass: find ring-closure bonds (back edges of the depth-first walk)
  const order = [], parentEdge = new Map();
  const walk = (a, from) => {
    seen.add(a); order.push(a);
    for (const [n, , k] of adj[a]) {
      if (k === from) continue;
      if (seen.has(n)) { if (!used.has(k)) { used.add(k); closures.set(k, ring++); } continue; }
      used.add(k); parentEdge.set(n, k); walk(n, k);
    }
  };
  const roots = [];                                        // a salt (SO3K) is two pieces: S... then the K+
  for (let a = 0; a < atoms.length; a++) if (!seen.has(a)) { roots.push(a); walk(a, -1); }
  const atomText = a => {
    const at = atoms[a], h = hcount[a];
    if (h === null && !at.charge && ORGANIC.has(at.el)) return at.el;
    const q = at.charge ? (at.charge > 0 ? '+' : '-') + (Math.abs(at.charge) > 1 ? Math.abs(at.charge) : '') : '';
    return `[${at.el}${h ? 'H' + (h > 1 ? h : '') : ''}${q}]`;
  };
  const done = new Set();
  const emit = a => {
    done.add(a);
    let s = atomText(a);
    for (const [n, o, k] of adj[a]) {                        // ring closures at this atom
      if (closures.has(k)) { const r = closures.get(k); s += sym(o) + (r > 9 ? '%' + r : r); }
    }
    const kids = adj[a].filter(([n, , k]) => parentEdge.get(n) === k && !done.has(n));
    kids.forEach(([n, o], idx) => {
      const sub = sym(o) + emit(n);
      s += idx < kids.length - 1 ? `(${sub})` : sub;
    });
    return s;
  };
  return roots.map(emit).join('.');
}

/** The group a condensed label stands for, as a SMILES whose first atom attaches; null when it is not one. */
export function condensedSmiles(raw) {
  const text = (raw || '').trim().replace('−', '-');
  if (!text || text.length > 40) return null;
  if (ALIASES[text]) return NICKNAMES[ALIASES[text]];
  if (cache.has(text)) return cache.get(text);
  let out = null;
  try { out = parse(text); } catch { out = null; }
  cache.set(text, out);
  return out;
}

// ------------------------------------------------------------------------------------------- lower case, radicals
const LOWER_PAIRS = ['cl', 'br', 'na', 'li', 'mg', 'si', 'me', 'et', 'ph', 'ac', 'bn', 'bz', 'bu', 'pr', 'ts', 'tf', 'cy'];

/** An all-lower-case label in capitals ('so3k' -> 'SO3K', 'ch2cl' -> 'CH2Cl', 'ome' -> 'OMe'), or null (model.case_guess). */
export function caseGuess(text, isElement) {
  if (!text || text !== text.toLowerCase() || !/[a-z]/.test(text)) return null;
  let out = '', k = 0;
  while (k < text.length) {
    const ch = text[k];
    if (!/[a-z]/.test(ch)) { out += ch; k++; continue; }
    const pair = text.slice(k, k + 2), cap = pair[0].toUpperCase() + pair.slice(1);
    if (pair.length === 2 && /^[a-z]{2}$/.test(pair) && LOWER_PAIRS.includes(pair)) { out += cap; k += 2; }
    else if (isElement(ch.toUpperCase())) { out += ch.toUpperCase(); k++; }
    else if (pair.length === 2 && /^[a-z]{2}$/.test(pair) && (isElement(cap) || NICKNAMES[cap])) { out += cap; k += 2; }
    else return null;
  }
  return out;
}
