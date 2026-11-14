// layout.js - 2D coordinates for a structure from a SMILES or a name, as the desktop's layout2d.py: CoordGen lays it out, a
// common skeleton is turned to its usual orientation and takes the template's exact shape, bonds are made even, what hangs
// on a ring goes into the widest gap and ring-fusion stereocentres show their H. RDKit's JavaScript build cannot wedge moved
// coordinates, so each wedge is chosen here and checked against the molecule's stereo.
import { LAYOUT, SAME_SHAPE } from './data.js';
import { buildBlock, parseBlock } from './doc.js';

const TAU = 2 * Math.PI;
const mod = (a, n) => ((a % n) + n) % n;
const median = v => { const s = [...v].sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; };
const noBondMarks = s => s.replace(/[/\\]/g, '');
const KEEP_H = JSON.stringify({ removeHs: false });
const queries = new WeakMap();

function compiled(RD) {
  if (!queries.has(RD)) {
    queries.set(RD, LAYOUT.map(([name, smarts, coords]) => ({ name, query: RD.get_qmol(smarts), coords }))
      .sort((p, q) => q.coords.length - p.coords.length));
  }
  return queries.get(RD);
}

function withMol(RD, text, fn, details) {
  const m = details ? RD.get_mol(text, details) : RD.get_mol(text);
  if (!m) return null;
  try { return fn(m); } finally { m.delete(); }
}

/** The SMILES of a drawn structure, without its drawn hydrogens (a MOL block keeps them unless told otherwise). */
function smilesOf(RD, g) {
  return withMol(RD, buildBlock(g.atoms, g.bonds), m => { m.remove_hs_in_place(); return m.get_smiles(); });
}

function adjacency(g) {
  const nb = g.atoms.map(() => []);
  g.bonds.forEach((b, k) => { nb[b.a].push([b.b, k]); nb[b.b].push([b.a, k]); });
  return nb;
}

/** The smallest ring through each bond, without repeats: for ordinary fused systems this is the SSSR. */
export function findRings(g) {
  const nb = adjacency(g), seen = new Set(), rings = [];
  g.bonds.forEach((bond, k) => {
    const prev = new Map([[bond.a, -1]]);
    let frontier = [bond.a], path = null;
    while (frontier.length && !path) {
      const nxt = [];
      for (const u of frontier) {
        for (const [n, j] of nb[u]) {
          if (j === k || prev.has(n)) continue;
          prev.set(n, u);
          if (n === bond.b) { path = [n]; while (prev.get(path[path.length - 1]) !== -1) path.push(prev.get(path[path.length - 1])); break; }
          nxt.push(n);
        }
        if (path) break;
      }
      frontier = nxt;
    }
    if (path) {
      const key = [...path].sort((p, q) => p - q).join();
      if (!seen.has(key)) { seen.add(key); rings.push(path); }
    }
  });
  return rings;
}

function length(g, b) { return Math.hypot(g.atoms[b.a].x - g.atoms[b.b].x, g.atoms[b.a].y - g.atoms[b.b].y); }

export function clashes(g) {
  const bond = g.bonds.length ? median(g.bonds.map(b => length(g, b))) : 1;
  const bonded = new Set(g.bonds.flatMap(b => [`${b.a},${b.b}`, `${b.b},${b.a}`]));
  let n = 0;
  for (let i = 0; i < g.atoms.length; i++) {
    for (let j = i + 1; j < g.atoms.length; j++) {
      if (Math.hypot(g.atoms[i].x - g.atoms[j].x, g.atoms[i].y - g.atoms[j].y) / (bond || 1) < 0.45 && !bonded.has(`${i},${j}`)) n++;
    }
  }
  return n;
}

/** Turns (and if needed mirrors) the whole layout onto the template; returns the skeleton's RMS distance from it. */
function align(g, fixed) {
  const ids = [...fixed.keys()], n = ids.length;
  const bond = median(g.bonds.map(b => length(g, b))) || 1, scale = 1 / bond;
  const mean = pts => [pts.reduce((s, p) => s + p[0], 0) / n, pts.reduce((s, p) => s + p[1], 0) / n];
  const have = ids.map(a => [g.atoms[a].x, g.atoms[a].y]), want = ids.map(a => fixed.get(a));
  const [hx, hy] = mean(have), [wx, wy] = mean(want);
  const fit = flip => {
    let a = 0, b = 0;
    have.forEach((p, k) => {
      const x = (p[0] - hx) * scale, y = (p[1] - hy) * scale * flip, u = want[k][0] - wx, v = want[k][1] - wy;
      a += x * u + y * v; b += x * v - y * u;
    });
    const t = Math.atan2(b, a), c = Math.cos(t), s = Math.sin(t);
    let e = 0;
    have.forEach((p, k) => {
      const x = (p[0] - hx) * scale, y = (p[1] - hy) * scale * flip;
      e += (x * c - y * s + wx - want[k][0]) ** 2 + (x * s + y * c + wy - want[k][1]) ** 2;
    });
    return { flip, c, s, rms: Math.sqrt(e / n) };
  };
  const a = fit(1), b = fit(-1), best = b.rms < a.rms ? b : a;
  for (const at of g.atoms) {
    const x = (at.x - hx) * scale, y = (at.y - hy) * scale * best.flip;
    at.x = x * best.c - y * best.s + wx;
    at.y = x * best.s + y * best.c + wy;
  }
  return best.rms;
}

function faceRingsLeft(g, rings) {
  const n = g.atoms.length;
  const mean = g.atoms.reduce((s, a) => s + a.x, 0) / n;
  const inRing = new Set(rings.flat());
  let off = 0;
  if (inRing.size) off = [...inRing].reduce((s, i) => s + g.atoms[i].x, 0) / inRing.size - mean;
  else {
    const hetero = g.atoms.filter(a => !['H', 'C'].includes(a.el));
    if (hetero.length) off = mean - hetero.reduce((s, a) => s + a.x, 0) / hetero.length;
  }
  if (off > 0.5) for (const a of g.atoms) a.x = -a.x;
}

function bridged(rings) {
  const out = new Set();
  rings.forEach((r, i) => rings.slice(i + 1).forEach(q => {
    if (r.filter(a => q.includes(a)).length > 2) for (const a of [...r, ...q]) out.add(a);
  }));
  return out;
}

function angleTargets(g, nb, i, ringSize, inRing) {
  const at = g.atoms[i];
  const list = nb[i].map(([n]) => n);
  if (list.length < 2) return [];
  list.sort((p, q) => Math.atan2(g.atoms[p].y - at.y, g.atoms[p].x - at.x) - Math.atan2(g.atoms[q].y - at.y, g.atoms[q].x - at.x));
  const fixed = [], free = [];
  list.forEach((p, k) => {
    const q = list[(k + 1) % list.length], size = ringSize(i, p, q);
    if (size) fixed.push([p, q, (size - 2) * Math.PI / size]); else free.push([p, q]);
  });
  if (inRing.has(i)) return fixed.map(([p, q, t]) => [p, q, 2 * Math.sin(t / 2)]);
  if (list.length === 2 && !fixed.length) {
    const linear = nb[i].some(([, b]) => g.bonds[b].order === 3);
    return [[list[0], list[1], linear ? 2 : 2 * Math.sin(Math.PI / 3)]];
  }
  const used = fixed.reduce((s, f) => s + f[2], 0);
  const share = free.length ? (TAU - used) / free.length : 0;
  return [...fixed, ...free.map(([p, q]) => [p, q, share])].filter(f => f[2] < 2.97).map(([p, q, t]) => [p, q, 2 * Math.sin(t / 2)]);
}

function branch(nb, start, awayFrom) {
  const seen = new Set([start]), stack = [start];
  while (stack.length) for (const [n] of nb[stack.pop()]) if (n !== awayFrom && !seen.has(n)) { seen.add(n); stack.push(n); }
  return [...seen];
}

const CROWD = 0.8;                                           // as layout2d.CROWD

/** Atoms laid almost on top of each other (two OH of a diphosphate): an end atom or a branch is turned about the atom
 * it hangs on, in 15 degree steps, to where nothing crowds it (layout2d.uncrowd). */
function uncrowd(g, fixed = new Set(), rounds = 30) {
  const n = g.atoms.length;
  if (!g.bonds.length) return;
  const nb = adjacency(g);
  const ringBond = new Set();
  for (const r of findRings(g)) r.forEach((a, k) => { const b = r[(k + 1) % r.length]; ringBond.add(`${a},${b}`); ringBond.add(`${b},${a}`); });
  const bond = median(g.bonds.map(b => length(g, b))) || 1;
  const lab = g.atoms.map(a => a.el !== 'C');
  const need = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) =>
    (i === j ? 0 : lab[i] && lab[j] ? 1.25 : lab[i] || lab[j] ? 0.9 : CROWD)));
  for (let i = 0; i < n; i++) for (const [p] of nb[i]) { need[i][p] = need[p][i] = 0; for (const [q] of nb[i]) need[p][q] = 0; }
  const xy = g.atoms.map(a => [a.x, a.y]);
  const score = pts => {
    let s = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (need[i][j]) {
      const d = Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]) / bond;
      if (d < need[i][j]) s += need[i][j] - d;
    }
    return s;
  };
  const squeeze = (pts, i) => {
    if (nb[i].length < 2) return 0;
    const ang = nb[i].map(([m]) => (Math.atan2(pts[m][1] - pts[i][1], pts[m][0] - pts[i][0]) + 2 * Math.PI) % (2 * Math.PI)).sort((a, b) => a - b);
    const gaps = ang.map((a, k) => ((ang[(k + 1) % ang.length] - a) + 2 * Math.PI) % (2 * Math.PI) || 2 * Math.PI);
    const lim = 50 * Math.PI / 180;
    return Math.max(0, lim - Math.min(...gaps)) / lim;
  };
  let total = score(xy), cur = xy;
  const tried = new Set();
  for (let r = 0; r < rounds && total > 1e-6; r++) {
    let worst = null;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) if (need[i][j] && !tried.has(`${i},${j}`)) {
      const d = Math.hypot(cur[i][0] - cur[j][0], cur[i][1] - cur[j][1]) / bond;
      if (d < need[i][j] && (!worst || need[i][j] - d > worst[0])) worst = [need[i][j] - d, i, j];
    }
    if (!worst) break;
    tried.add(`${worst[1]},${worst[2]}`);
    let best = null;
    for (const x of [worst[1], worst[2]]) {
      for (const [pivot] of nb[x]) {
        if (ringBond.has(`${x},${pivot}`)) continue;
        const side = branch(nb, x, pivot);
        if (side.some(a => fixed.has(a)) || side.length > n / 2) continue;
        const base = squeeze(cur, pivot);
        for (let step = 1; step < 24; step++) {
          const turn = (15 * Math.floor((step + 1) / 2) * (step % 2 ? 1 : -1)) * Math.PI / 180, c = Math.cos(turn), s = Math.sin(turn);
          const trial = cur.map(p => [p[0], p[1]]);
          for (const a of side) {
            const dx = cur[a][0] - cur[pivot][0], dy = cur[a][1] - cur[pivot][1];
            trial[a] = [cur[pivot][0] + dx * c - dy * s, cur[pivot][1] + dx * s + dy * c];
          }
          const sc = score(trial) + 2 * Math.max(0, squeeze(trial, pivot) - base);
          if (!best || sc < best[0] - 1e-9) best = [sc, trial];
        }
      }
    }
    if (best && best[0] < total - 1e-3) { total = best[0]; cur = best[1]; tried.clear(); }
  }
  cur.forEach((p, i) => { g.atoms[i].x = p[0]; g.atoms[i].y = p[1]; });
}

function placeSubstituents(g, nb, frozen, rings) {
  if (!rings.length) return;
  const ringBond = new Set();
  for (const r of rings) r.forEach((a, k) => { const b = r[(k + 1) % r.length]; ringBond.add(`${a},${b}`); ringBond.add(`${b},${a}`); });
  const inRing = new Set(rings.flat());
  const start = frozen.size ? [...frozen] : [...rings].sort((p, q) => q.length - p.length)[0];
  const seen = new Set(start), queue = [...start];
  while (queue.length) {
    const i = queue.shift(), at = g.atoms[i], taken = [], children = [];
    for (const [n] of nb[i]) {
      if (ringBond.has(`${i},${n}`) || seen.has(n)) {
        taken.push(n);
        if (!seen.has(n)) { seen.add(n); queue.push(n); }
      } else children.push(n);
    }
    if (!inRing.has(i)) { for (const n of children) { seen.add(n); queue.push(n); } continue; }
    if (!children.length || !taken.length) continue;
    const angle = n => mod(Math.atan2(g.atoms[n].y - at.y, g.atoms[n].x - at.x), TAU);
    const have = taken.map(angle).sort((p, q) => p - q);
    let gap = -1, low = 0;
    have.forEach((h, k) => {
      const w = mod(have[(k + 1) % have.length] - h, TAU) || TAU;
      if (w > gap) { gap = w; low = h; }
    });
    const step = Math.min(gap / (children.length + 1), 70 * Math.PI / 180), middle = low + gap / 2;
    children.sort((p, q) => mod(angle(p) - low, TAU) - mod(angle(q) - low, TAU));
    children.forEach((child, k) => {
      const target = middle + (k - (children.length - 1) / 2) * step, turn = target - angle(child);
      const c = Math.cos(turn), s = Math.sin(turn), ch = g.atoms[child];
      const fx = ch.x - at.x, fy = ch.y - at.y, scale = 1 / (Math.hypot(fx, fy) || 1);
      const nx = at.x + (fx * c - fy * s) * scale, ny = at.y + (fx * s + fy * c) * scale, ox = ch.x, oy = ch.y;
      for (const m of branch(nb, child, i)) {
        const a = g.atoms[m], rx = a.x - ox, ry = a.y - oy;
        a.x = nx + rx * c - ry * s; a.y = ny + rx * s + ry * c;
      }
      seen.add(child); queue.push(child);
    });
  }
}

/** In a bay of fused rings the gaps of two neighbouring atoms point at the same spot (C8 and C14 of a triterpene, each with
 * its methyl): both turn away from each other, as chemists draw them. */
function spreadBays(g, nb, inRing) {
  const parentOf = i => {
    if (inRing.has(i)) return null;
    const towards = nb[i].map(([n]) => n).filter(n => inRing.has(n));
    return towards.length === 1 ? towards[0] : null;
  };
  const turn = (i, parent, angle) => {
    const P = g.atoms[parent], c = Math.cos(angle), s = Math.sin(angle);
    for (const m of branch(nb, i, parent)) {
      const a = g.atoms[m], x = a.x - P.x, y = a.y - P.y;
      a.x = P.x + x * c - y * s; a.y = P.y + x * s + y * c;
    }
  };
  for (let round = 0; round < 3; round++) {
    for (let i = 0; i < g.atoms.length; i++) {
      for (let j = i + 1; j < g.atoms.length; j++) {
        const A = g.atoms[i], B = g.atoms[j];
        if (Math.hypot(A.x - B.x, A.y - B.y) >= 0.75) continue;
        const pa = parentOf(i), pb = parentOf(j);
        if (pa === null || pb === null || pa === pb) continue;
        const away = (k, parent, other) => {          // the sign that takes atom k further from the other parent
          const P = g.atoms[parent], O = g.atoms[other], K = g.atoms[k], t = 25 * Math.PI / 180;
          const at = sign => {
            const x = K.x - P.x, y = K.y - P.y, c = Math.cos(sign * t), s = Math.sin(sign * t);
            return Math.hypot(P.x + x * c - y * s - O.x, P.y + x * s + y * c - O.y);
          };
          return at(1) >= at(-1) ? t : -t;
        };
        turn(i, pa, away(i, pa, pb));
        turn(j, pb, away(j, pb, pa));
      }
    }
  }
}

/** Every bond one length and the usual angles; atoms in fixed and cages stay. Undone if atoms would come to overlap. */
function regularise(g, rings, fixed = new Set(), rounds = 400) {
  if (!g.bonds.length) return;
  const before = clashes(g), original = g.atoms.map(a => [a.x, a.y]);
  if (!fixed.size) {
    const L = median(g.bonds.map(b => length(g, b))) || 1;
    for (const a of g.atoms) { a.x /= L; a.y /= L; }
  }
  const frozen = new Set([...fixed, ...bridged(rings)]);
  const nb = adjacency(g), inRing = new Set(rings.flat());
  const ringSize = (i, p, q) => {
    const sizes = rings.filter(r => r.includes(i) && r.includes(p) && r.includes(q)).map(r => r.length);
    return sizes.length && Math.min(...sizes) <= 8 ? Math.min(...sizes) : 0;
  };
  const pairs = g.bonds.map(b => [b.a, b.b, 1]);
  g.atoms.forEach((_a, i) => pairs.push(...angleTargets(g, nb, i, ringSize, inRing)));
  const move = g.atoms.map((_a, i) => (frozen.has(i) ? 0 : 1));
  for (let r = 0; r < rounds; r++) {
    for (const [p, q, want] of pairs) {
      const total = move[p] + move[q];
      if (!total) continue;
      const A = g.atoms[p], B = g.atoms[q], dx = B.x - A.x, dy = B.y - A.y, now = Math.hypot(dx, dy);
      if (now < 1e-6) continue;
      const k = (now - want) / now * 0.5;
      A.x += dx * k * move[p] * 2 / total; A.y += dy * k * move[p] * 2 / total;
      B.x -= dx * k * move[q] * 2 / total; B.y -= dy * k * move[q] * 2 / total;
    }
  }
  placeSubstituents(g, nb, frozen, rings);
  spreadBays(g, nb, inRing);
  if (clashes(g) > before) g.atoms.forEach((a, i) => { a.x = original[i][0]; a.y = original[i][1]; });
}

// as in betanin (E): the carboxyl next to N+ on the chain's side of the C=N+ bond, unless the SMILES states that bond
function betalainCarboxylateByTheChain(RD, mol, g) {
  const q = RD.get_qmol('[#6;!R]=[#7+;R](-[#6;R]-[#6](~[#8])~[#8])');
  const hit = JSON.parse(mol.get_substruct_match(q));
  q.delete();
  if (!hit.atoms) return;
  const [chain, nitrogen, acid] = hit.atoms, nb = adjacency(g);
  const before = nb[chain].map(([n]) => n).filter(n => n !== nitrogen && g.atoms[n].el !== 'H');
  if (!before.length) return;
  const N = g.atoms[nitrogen], C = g.atoms[chain];
  let ax = N.x - C.x, ay = N.y - C.y;
  const L = Math.hypot(ax, ay) || 1; ax /= L; ay /= L;
  const side = (m, base) => ax * (g.atoms[m].y - base.y) - ay * (g.atoms[m].x - base.x);
  if (side(before[0], C) * side(acid, N) > 0) return;
  for (const m of branch(nb, nitrogen, chain)) {
    const a = g.atoms[m], rx = a.x - N.x, ry = a.y - N.y, d = rx * ax + ry * ay;
    a.x = N.x + 2 * d * ax - rx; a.y = N.y + 2 * d * ay - ry;
  }
}

function betalainTurnedOver(g) {                           // betalains: over the page's vertical axis, ring N still at the bottom
  for (const a of g.atoms) a.x = -a.x;
}

function templateFor(RD, mol, g) {
  for (const { name, query, coords } of compiled(RD)) {
    const matches = JSON.parse(mol.get_substruct_matches(query, JSON.stringify({ uniquify: false, maxMatches: 64 })));
    if (!matches.length) continue;
    const list = matches.map(m => m.atoms);
    let pick = list[0];
    if (name === 'oleanane') {                             // it reads the same from both ends: ring A has C3-O and C4 methyls
      const nb = adjacency(g), rings = findRings(g), inRing = new Set(rings.flat());
      const score = m => 2 * (nb[m[2]].some(([n]) => g.atoms[n].el === 'O') ? 1 : 0)
        + nb[m[3]].filter(([n]) => g.atoms[n].el === 'C' && !inRing.has(n)).length;
      pick = list.reduce((b, m) => (score(m) > score(b) ? m : b), list[0]);
    }
    return { name, fixed: new Map(pick.map((a, k) => [a, [coords[k][0], -coords[k][1]]])) };   // templates have y up
  }
  return null;
}

/** The wedges: each stereocentre gets one on its best bond, whichever of wedge / hash gives the molecule's stereo. */
function chooseWedges(RD, g, targetTags, target) {
  const nb = adjacency(g), rings = findRings(g), ringBond = new Set();
  for (const r of rings) r.forEach((a, k) => { const b = r[(k + 1) % r.length]; ringBond.add(`${a},${b}`); ringBond.add(`${b},${a}`); });
  const chiral = new Set(targetTags.map(t => t[0]));
  for (const b of g.bonds) b.stereo = '';
  const tagsOf = () => withMol(RD, buildBlock(g.atoms, g.bonds), m => new Map(JSON.parse(m.get_stereo_tags()).CIP_atoms), KEEP_H) || new Map();
  const taken = new Set(), placed = [];
  for (const [i, tag] of targetTags) {
    const options = nb[i].filter(([, k]) => !taken.has(k) && g.bonds[k].order === 1).map(([n, k]) => {
      const ring = ringBond.has(`${i},${n}`), atom = g.atoms[n];
      const score = (atom.el === 'H' ? 0 : ring ? 40 : 10) + (chiral.has(n) ? 20 : 0) + nb[n].length;
      return [score, k, n];
    }).sort((p, q) => p[0] - q[0]);
    if (!options.length) continue;
    const [, k, n] = options[0], bond = g.bonds[k];
    bond.a = i; bond.b = n; bond.stereo = 'wedge';
    placed.push([i, tag, bond]);
    taken.add(k);
  }
  // R / S can only be read once every centre has a wedge: flip those that come out wrong until all agree
  for (let pass = 0; pass < 4; pass++) {
    const tags = tagsOf();
    const wrong = placed.filter(([i, tag]) => tags.get(i) !== tag);
    if (!wrong.length) break;
    for (const [, , bond] of wrong) bond.stereo = bond.stereo === 'wedge' ? 'hash' : 'wedge';
  }
  if (!placed.length) return true;                         // no stereocentres: nothing to get wrong
  const same = () => noBondMarks(smilesOf(RD, g) || '') === noBondMarks(target);
  if (same()) return true;
  for (let pass = 0; pass < 2; pass++) {                   // a centre whose R / S hangs on another: try each the other way
    for (const k of taken) {
      const bond = g.bonds[k], was = bond.stereo;
      bond.stereo = was === 'wedge' ? 'hash' : 'wedge';
      if (same()) return true;
      bond.stereo = was;
    }
  }
  return false;
}

function fusionHydrogens(g, targetTags, rings) {
  const nb = adjacency(g), ringBond = new Set();
  for (const r of rings) r.forEach((a, k) => { const b = r[(k + 1) % r.length]; ringBond.add(`${a},${b}`); ringBond.add(`${b},${a}`); });
  for (const [i] of targetTags) {
    const at = g.atoms[i];
    if (at.el === 'C' && !at.charge && nb[i].length === 3 && nb[i].every(([n, k]) => g.bonds[k].order === 1 && ringBond.has(`${i},${n}`))) {
      g.bonds.push({ a: i, b: g.atoms.length, order: 1, stereo: '' });
      g.atoms.push({ x: at.x + 0.2, y: at.y + 0.2, el: 'H', charge: 0 });
    }
  }
}

function attempt(RD, smiles, gentle) {
  const mol = RD.get_mol(smiles);
  if (!mol) return null;
  try {
    const target = mol.get_smiles();
    const tags = JSON.parse(mol.get_stereo_tags()).CIP_atoms.filter(t => !t[1].includes('?'));   // '?': not specified
    RD.prefer_coordgen(true);
    mol.set_new_coords(true);
    const g = parseBlock(mol.get_molblock());
    for (const b of g.bonds) b.stereo = '';
    const rings = findRings(g);
    const found = templateFor(RD, mol, g);
    const heavy = g.atoms.length;
    if (g.bonds.length) {
      if (!found) {
        faceRingsLeft(g, rings);
        const L = median(g.bonds.map(b => length(g, b))) || 1;
        for (const a of g.atoms) { a.x /= L; a.y /= L; }
      } else if (align(g, found.fixed) < SAME_SHAPE) {
        for (const [a, p] of found.fixed) { g.atoms[a].x = p[0]; g.atoms[a].y = p[1]; }
      } else found.fixed.clear();
    }
    fusionHydrogens(g, tags, rings);
    if (!gentle) {
      const rings2 = findRings(g);
      regularise(g, rings2, found ? new Set(found.fixed.keys()) : new Set());
      uncrowd(g, found ? new Set(found.fixed.keys()) : new Set());
      if (found && found.name === 'betalamic') betalainCarboxylateByTheChain(RD, mol, g);
    } else if (g.atoms.length > heavy) {
      const nb = adjacency(g);                             // the hydrogens still need a place
      placeSubstituents(g, nb, new Set(g.atoms.slice(0, heavy).map((_a, i) => i)), findRings(g));
    }
    if (found && found.name === 'betalamic') betalainTurnedOver(g);
    const ok = chooseWedges(RD, g, tags, target);
    const ez = !/[/\\]/.test(target) || smilesOf(RD, g) === target;
    return { atoms: g.atoms, bonds: g.bonds, stereoOk: ok && ez, template: found ? found.name : '', clashes: clashes(g) };
  } finally { mol.delete(); }
}

/** {atoms [{el, charge, x, y}], bonds [{a, b, order, stereo}], stereoOk, template, clashes} of a SMILES, or null. */
export function layoutSmiles(RD, smiles) {
  const first = attempt(RD, smiles, false);
  if (!first || first.stereoOk) return first;
  const second = attempt(RD, smiles, true);                // moving atoms changed an E/Z or R/S: CoordGen's own coordinates
  return second && second.stereoOk ? second : first;
}

/** A betalain as PubChem may give it (ring N=C, chain amine) in its usual form: ring NH, iminium N+ and a carboxylate. */
export function usualForm(RD, smiles) {
  const mol = RD.get_mol(smiles);
  if (!mol) return smiles;
  const q = RD.get_qmol('[#7;R;X2]=[#6;R]-[#6;R]=[#6;R]-[#6;!R]=[#6;!R]-[#7;X3;+0]');
  const acidQ = RD.get_qmol('[#6]-[#6](=O)-[OX2H1]');
  try {
    const hit = JSON.parse(mol.get_substruct_match(q));
    if (!hit.atoms) return smiles;
    const [n1, c2, c3, c4, cc, cb, amine] = hit.atoms;
    const acids = JSON.parse(mol.get_substruct_matches(acidQ)).map(m => m.atoms);
    RD.prefer_coordgen(true);
    mol.set_new_coords(true);
    const g = parseBlock(mol.get_molblock());
    const setOrder = (a, b, order) => {
      const bd = g.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a));
      bd.order = order; bd.stereo = '';
    };
    setOrder(n1, c2, 1); setOrder(c2, c3, 2); setOrder(c3, c4, 1); setOrder(c4, cc, 2); setOrder(cc, cb, 1); setOrder(cb, amine, 2);
    g.atoms[amine].charge = 1;
    const near = new Set(g.bonds.flatMap(b => (b.a === amine ? [b.b] : b.b === amine ? [b.a] : [])));
    const usable = acids.filter(m => m[0] === c2 || near.has(m[0])).sort((p, r) => (p[0] === c2) - (r[0] === c2));
    if (usable.length) g.atoms[usable[0][3]].charge = -1;
    const out = withMol(RD, buildBlock(g.atoms, g.bonds), m => noBondMarks(m.get_smiles()));
    return out || smiles;
  } finally { q.delete(); acidQ.delete(); mol.delete(); }
}
