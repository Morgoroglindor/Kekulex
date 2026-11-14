// app.js - Kekulex on a phone or tablet: the touch tools, the page, the sheets. The chemistry is in core/.
import { analyse, insertLayout, insertMolblock, insertSmiles, labelExtras, parseLabel, putGroup, refreshLiveCaptions, setNickname } from './core/chem.js';
import { setRD } from './core/condensed.js';
import { GROUPS, NICKNAMES, TEMPLATES } from './core/data.js';
import { Doc, snapAngle } from './core/doc.js';
import { bounds, DRAWN_SHAPES, paint, toCanvas } from './core/draw.js';
import { layoutSmiles, usualForm } from './core/layout.js';

const $ = id => document.getElementById(id);
const canvas = $('canvas'), ctx = canvas.getContext('2d');
const STORE = 'kekulex.autosave';
const HIT_PX = 24, DRAG_PX = 9, LONG_MS = 480, MAX_UNDO = 60;
const COMMON = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'H', 'Si', 'B', 'Na', 'K', 'Li', 'Mg'];
const QUICK_GROUPS = ['Me', 'Et', 'iPr', 'tBu', 'Ph', 'Bn', 'OMe', 'OEt', 'OH', 'NH2', 'NO2', 'CN', 'CF3', 'CO2H', 'CO2Me', 'CHO', 'Ac', 'OAc', 'Boc', 'Cbz', 'Fmoc', 'Ts', 'Ms', 'Tf', 'TMS', 'TBS'];

let RD = null;
let doc = new Doc();
let undo = [], redo = [];
let path = 'structure';
const view = { bl: 42, ox: 0, oy: 0, dpr: 1 };
const state = { tool: 'bond', bond: 'single', ring: 'benzene', el: 'N', sel: new Set(), press: null, preview: null, hint: null };
const pointers = new Map();
let pinch = null, longTimer = null, needsPaint = false, infoTimer = null, toastTimer = null, sheetOnClose = null;

// ---------------------------------------------------------------------------------------------- icons
const svg = inner => `<svg viewBox="0 0 32 32">${inner}</svg>`;
const poly = (n, r = 11, start = -90) => Array.from({ length: n }, (_, k) => { const a = (start + 360 * k / n) * Math.PI / 180; return `${(16 + r * Math.cos(a)).toFixed(1)},${(16 + r * Math.sin(a)).toFixed(1)}`; }).join(' ');
const ICONS = {
  select: svg('<rect x="6" y="6" width="20" height="20" rx="2" stroke-dasharray="4 3"/><path d="M13 12l9 6-4 1 2 5-2 1-2-5-3 3z" class="solid"/>'),
  bond: svg('<path d="M7 25L25 7"/>'),
  single: svg('<path d="M7 25L25 7"/>'),
  double: svg('<path d="M5 22L22 5M10 27L27 10"/>'),
  triple: svg('<path d="M7 25L25 7M3 21L21 3M11 29L29 11"/>'),
  wedge: svg('<path class="solid" d="M6 26L24 5l3.5 3z"/>'),
  hash: svg('<path d="M9 23l1.2 1.2M12 19l2.4 2.4M15 15l3.6 3.6M18 11l4.8 4.8M21 7l6 6"/>'),
  wavy: svg('<path d="M5 26c3-1 2-4 4-5s3 1 5-1 0-4 2-5 3 1 5-1 0-4 2-5 3 0 4-2"/>'),
  bold: svg('<path d="M7 25L25 7" style="stroke-width:4.4"/>'),
  dashed: svg('<path d="M7 25L25 7" stroke-dasharray="3.4 3.4"/>'),
  chain: svg('<path d="M3 21l6.5-9 6.5 9 6.5-9L29 21"/>'),
  ring: svg(`<polygon points="${poly(6)}"/>`),
  benzene: svg(`<polygon points="${poly(6)}"/><path d="M16 9l5.7 3.4M21.7 19.3L16 22.7M10.3 19.3v-6.8"/>`),
  r3: svg(`<polygon points="${poly(3)}"/>`), r4: svg(`<polygon points="${poly(4, 11, -45)}"/>`), r5: svg(`<polygon points="${poly(5)}"/>`),
  r6: svg(`<polygon points="${poly(6)}"/>`), r7: svg(`<polygon points="${poly(7)}"/>`), r8: svg(`<polygon points="${poly(8)}"/>`),
  erase: svg('<path d="M12 26h14M6 19l10-12 10 8-9 11h-5z"/><path d="M11 13l10 8"/>'),
  charge: svg('<circle cx="16" cy="16" r="11"/><path d="M16 10v12M10 16h12"/>'),
  arrow: svg('<path d="M4 16h22"/><path class="solid" d="M28 16l-9-5 2 5-2 5z"/>'),
  text: svg('<path d="M7 8h18M16 8v17M12 25h8"/>'),
  plus: svg('<path d="M16 7v18M7 16h18"/>'),
  minus: svg('<path d="M7 16h18"/>'),
  rad: svg('<circle class="solid" cx="16" cy="16" r="3.2"/>'),
  dirad: svg('<circle class="solid" cx="11" cy="16" r="3.2"/><circle class="solid" cx="21" cy="16" r="3.2"/>'),
  'rad+': svg('<circle class="solid" cx="8" cy="16" r="3.2"/><path d="M20 9v14M13 16h14"/>'),
  'rad-': svg('<circle class="solid" cx="8" cy="16" r="3.2"/><path d="M13 16h14"/>'),
};
// The charge tool's choices: [key, name, unpaired electrons, charge] (null: the charge is only added to)
const CHARGES = [['1', 'Plus', null, 1], ['-1', 'Minus', null, -1], ['rad', 'Radical', 1, null], ['dirad', 'Diradical', 2, null],
  ['rad+', 'Radical cation', 1, 1], ['rad-', 'Radical anion', 1, -1]];

const TOOLS = [['select', 'Select and move'], ['bond', 'Bonds'], ['ring', 'Rings'], ['atom', 'Atoms'], ['chain', 'Chain'], ['charge', 'Charge'], ['erase', 'Eraser'], ['arrow', 'Arrow'], ['text', 'Text']];
const BONDS = ['single', 'double', 'triple', 'wedge', 'hash', 'wavy', 'bold', 'dashed'];
const RINGS = ['benzene', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8'];
const ARROWS = [['forward', '→'], ['equilibrium', '⇌'], ['retro', '⇒'], ['resonance', '↔'], ['dashed', '⇢']];

function toast(text, ms = 3200) {
  const t = $('toast');
  t.textContent = text; t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { t.hidden = true; }, ms);
}

/** The charge tool on an atom: + / − add to its charge; a radical choice gives its dots (and charge), or, tapped
 * again on an atom that already has them, takes them away (as the desktop's radical tools do). */
function applyCharge(at, k) {
  const [, , dots, z] = CHARGES.find(c => c[0] === k) || CHARGES[0];
  if (dots === null) { at.charge += z; return; }
  if (at.radical === dots && (z === null || at.charge === z)) { at.radical = 0; if (z !== null) at.charge = 0; }
  else { at.radical = dots; if (z !== null) at.charge = z; }
}

// ---------------------------------------------------------------------------------------------- the tool bars
function renderTools() {
  $('toolRow').innerHTML = TOOLS.map(([k, label]) => `<button class="tool ${state.tool === k ? 'on' : ''}" data-tool="${k}" aria-label="${label}" title="${label}">${k === 'atom' ? (state.el.length > 3 ? 'Atom' : state.el) : k === 'bond' ? ICONS[state.bond] : k === 'ring' ? ICONS[state.ring] : ICONS[k]}</button>`).join('');
  const o = $('options');
  if (state.tool === 'bond') o.innerHTML = BONDS.map(b => `<button class="tool ${state.bond === b ? 'on' : ''}" data-bond="${b}" aria-label="${b} bond">${ICONS[b]}</button>`).join('');
  else if (state.tool === 'ring') o.innerHTML = RINGS.map(r => `<button class="tool ${state.ring === r ? 'on' : ''}" data-ring="${r}" aria-label="${r === 'benzene' ? 'Benzene' : `${r[1]}-membered ring`}">${ICONS[r]}</button>`).join('');
  else if (state.tool === 'atom') o.innerHTML = COMMON.map(e => `<button class="tool ${state.el === e ? 'on' : ''}" data-el="${e}">${e}</button>`).join('') + '<button class="tool" data-el-more="1">More…</button>';
  else if (state.tool === 'charge') o.innerHTML = CHARGES.map(([k, name]) => `<button class="tool ${(state.sign || '1') === k ? 'on' : ''}" data-sign="${k}" aria-label="${name}" title="${name}">${ICONS[k === '1' ? 'plus' : k === '-1' ? 'minus' : k]}</button>`).join('');
  else if (state.tool === 'arrow') o.innerHTML = ARROWS.map(([k, g]) => `<button class="tool ${(state.arrow || 'forward') === k ? 'on' : ''}" data-arrow="${k}" aria-label="${k} arrow" style="font-size:22px">${g}</button>`).join('');
  else o.innerHTML = '';
  requestAnimationFrame(resize);
}

function setTool(k) {
  state.tool = k;
  if (k !== 'select') { state.sel.clear(); renderSelBar(); }
  renderTools(); schedule();
}

// ---------------------------------------------------------------------------------------------- history
function snapshot() {
  undo.push(JSON.stringify(doc.toDict()));
  if (undo.length > MAX_UNDO) undo.shift();
  redo = [];
}
function edited() {
  doc.pruneGroups();
  doc.updateAttachedArrows();                              // curved arrows drawn on the desktop stay on their atoms
  state.sel = new Set([...state.sel].filter(exists));
  renderSelBar(); schedule(); syncHistory();
  clearTimeout(infoTimer);
  infoTimer = setTimeout(() => {
    try { refreshLiveCaptions(RD, doc); schedule(); } catch { /* a caption that cannot be worked out keeps its text */ }
    updateFormula(); autosave();
  }, 120);
}
function exists(key) {
  const [kind, id] = parseKey(key);
  return { atom: doc.atoms, bond: doc.bonds, arrow: doc.arrows, text: doc.texts }[kind].has(id);
}
const key = (kind, id) => `${kind}:${id}`;
const parseKey = k => { const [a, b] = k.split(':'); return [a, Number(b)]; };
function restore(text) { doc = Doc.fromDict(JSON.parse(text)); state.sel.clear(); edited(); }
function doUndo() { if (undo.length) { redo.push(JSON.stringify(doc.toDict())); restore(undo.pop()); } }
function doRedo() { if (redo.length) { undo.push(JSON.stringify(doc.toDict())); restore(redo.pop()); } }
function syncHistory() { $('undoBtn').disabled = !undo.length; $('redoBtn').disabled = !redo.length; }
function autosave() { try { localStorage.setItem(STORE, JSON.stringify(doc.toDict())); } catch { /* full: the file menu still saves */ } }

function updateFormula() {
  const chip = $('formulaChip'), out = $('formulaText');
  if (!RD || !doc.atoms.size) { out.textContent = 'Kekulex'; chip.classList.remove('bad'); return; }
  const ids = selectedAtoms();
  const info = analyse(RD, doc, ids.size ? ids : null);
  chip.classList.toggle('bad', info.errors.length > 0);
  out.innerHTML = `${info.formulaHtml}<span class="mass">${info.exact.toFixed(4)}</span>`;
}

// ---------------------------------------------------------------------------------------------- the view
function resize() {
  const r = canvas.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  const w = Math.round(r.width * dpr), h = Math.round(r.height * dpr);
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; view.dpr = dpr; schedule(); }
}
function schedule() { if (!needsPaint) { needsPaint = true; requestAnimationFrame(render); } }
const toModel = (x, y) => [(x - view.ox) / view.bl, (y - view.oy) / view.bl];
const toPx = (x, y) => [x * view.bl + view.ox, y * view.bl + view.oy];

function fit() {
  const r = canvas.getBoundingClientRect(), b = bounds(doc);
  if (!b) { view.bl = 42; view.ox = r.width / 2; view.oy = r.height / 2; schedule(); return; }
  // labels reach past their atoms: 1.1 bond lengths of room on each side, plus a little for the finger
  const w = Math.max(b[2] - b[0], 1) + 2.2, h = Math.max(b[3] - b[1], 1) + 2.2;
  view.bl = Math.min(64, Math.max(10, Math.min((r.width - 12) / w, (r.height - 24) / h)));
  view.ox = r.width / 2 - (b[0] + b[2]) / 2 * view.bl;
  view.oy = r.height / 2 - (b[1] + b[3]) / 2 * view.bl;
  schedule();
}

function selBox() {
  const pts = [];
  for (const k of state.sel) {
    const [kind, id] = parseKey(k);
    if (kind === 'atom') { const a = doc.atoms.get(id); if (a) pts.push([a.x, a.y]); }
    if (kind === 'arrow') { const a = doc.arrows.get(id); if (a) pts.push([a.x1, a.y1], [a.x2, a.y2]); }
    if (kind === 'text') { const t = doc.texts.get(id); if (t) pts.push([t.x, t.y], [t.x + String(t.text).length * 0.4, t.y + 0.8]); }
  }
  if (!pts.length) return null;
  const xs = pts.map(p => p[0]), ys = pts.map(p => p[1]), m = 0.35;
  return [Math.min(...xs) - m, Math.min(...ys) - m, Math.max(...xs) + m, Math.max(...ys) + m];
}

function render() {
  needsPaint = false;
  const d = window.devicePixelRatio || 1;
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.setTransform(d, 0, 0, d, view.ox * d, view.oy * d);
  const BL = view.bl;
  // the selection behind the drawing
  if (state.sel.size) {
    ctx.fillStyle = '#1b6ca833'; ctx.strokeStyle = '#1b6ca833'; ctx.lineWidth = 0.3 * BL; ctx.lineCap = 'round';
    for (const b of doc.bonds.values()) {
      if (state.sel.has(key('atom', b.a)) && state.sel.has(key('atom', b.b))) { const p = doc.atoms.get(b.a), q = doc.atoms.get(b.b); ctx.beginPath(); ctx.moveTo(p.x * BL, p.y * BL); ctx.lineTo(q.x * BL, q.y * BL); ctx.stroke(); }
    }
    const box = selBox();
    if (box) { ctx.strokeStyle = '#1b6ca8'; ctx.lineWidth = 1; ctx.setLineDash([5, 4]); ctx.strokeRect(box[0] * BL, box[1] * BL, (box[2] - box[0]) * BL, (box[3] - box[1]) * BL); ctx.setLineDash([]); }
  }
  paint(ctx, doc, BL, { errors: true });
  const pv = state.preview;
  if (pv) {
    ctx.strokeStyle = '#1b6ca8'; ctx.lineWidth = 2; ctx.lineCap = 'round';
    if (pv.band) { ctx.setLineDash([5, 4]); const [x0, y0, x1, y1] = pv.band; ctx.strokeRect(x0 * BL, y0 * BL, (x1 - x0) * BL, (y1 - y0) * BL); ctx.setLineDash([]); }
    if (pv.points) { ctx.beginPath(); pv.points.forEach(([x, y], k) => (k ? ctx.lineTo(x * BL, y * BL) : ctx.moveTo(x * BL, y * BL))); ctx.stroke(); }
  }
  if (state.hint) { ctx.fillStyle = '#1b6ca844'; ctx.beginPath(); ctx.arc(state.hint[0] * BL, state.hint[1] * BL, 0.32 * BL, 0, 7); ctx.fill(); }
  if (doc.isEmpty() && !pv) {
    const r = canvas.getBoundingClientRect();
    ctx.setTransform(d, 0, 0, d, 0, 0);
    ctx.fillStyle = '#8894a6'; ctx.font = '15px Roboto, Arial, sans-serif'; ctx.textAlign = 'center';
    ctx.fillText('Tap to start drawing', r.width / 2, r.height / 2 - 12);
    ctx.fillText('Two fingers: move and zoom', r.width / 2, r.height / 2 + 12);
    ctx.textAlign = 'left';
  }
}

// ---------------------------------------------------------------------------------------------- what is under a finger
function itemAt(mx, my) {
  const r = HIT_PX / view.bl;
  const a = doc.findAtom(mx, my, Math.min(r, 0.5));
  if (a !== null) return ['atom', a];
  const b = doc.findBond(mx, my, Math.min(r, 0.34));
  if (b !== null) return ['bond', b];
  const ar = doc.findArrow(mx, my, Math.min(r, 0.4));
  if (ar !== null) return ['arrow', ar];
  for (const [i, t] of doc.texts) {
    if (mx >= t.x - 0.1 && mx <= t.x + String(t.text).length * 0.45 && my >= t.y - 0.1 && my <= t.y + 0.9) return ['text', i];
  }
  return null;
}

function selectedAtoms() {
  const out = new Set();
  for (const k of state.sel) { const [kind, id] = parseKey(k); if (kind === 'atom') out.add(id); }
  return out;
}

function selectMolecule(item) {
  const add = new Set();
  if (item[0] === 'atom' || item[0] === 'bond') {
    const a = item[0] === 'atom' ? item[1] : doc.bonds.get(item[1]).a;
    for (const id of doc.fragmentOf(a)) add.add(key('atom', id));
  } else add.add(key(...item));
  return doc.withGroups(add);                               // a grouped structure comes with its group
}

// ---------------------------------------------------------------------------------------------- touch
function onDown(e) {
  canvas.setPointerCapture(e.pointerId);
  pointers.set(e.pointerId, [e.offsetX, e.offsetY]);
  if (pointers.size === 2) {                                // a second finger: move and zoom, what the first began is off
    cancelPress();
    const [p, q] = [...pointers.values()];
    pinch = { d: Math.hypot(p[0] - q[0], p[1] - q[1]) || 1, cx: (p[0] + q[0]) / 2, cy: (p[1] + q[1]) / 2, bl: view.bl, ox: view.ox, oy: view.oy };
    return;
  }
  if (pointers.size > 2) return;
  const [mx, my] = toModel(e.offsetX, e.offsetY);
  const item = itemAt(mx, my);
  state.press = { x: e.offsetX, y: e.offsetY, mx, my, item, moved: false, snapped: false };
  if (item) { state.hint = item[0] === 'atom' ? [doc.atoms.get(item[1]).x, doc.atoms.get(item[1]).y] : null; schedule(); }
  if (state.tool === 'select') {
    const box = selBox();
    if (item && !state.sel.has(key(...normalItem(item))) && !e.shiftKey) { /* decided on release or drag */ }
    state.press.inBox = !!box && mx >= box[0] && mx <= box[2] && my >= box[1] && my <= box[3];
  }
  clearTimeout(longTimer);
  if (item && item[0] === 'atom') longTimer = setTimeout(() => {
    if (state.press && !state.press.moved) { const a = item[1]; cancelPress(); labelDialog(a); }
  }, LONG_MS);
}

const normalItem = it => (it[0] === 'bond' ? ['atom', doc.bonds.get(it[1]).a] : it);

function cancelPress() {
  clearTimeout(longTimer);
  state.press = null; state.preview = null; state.hint = null; schedule();
}

function onMove(e) {
  if (!pointers.has(e.pointerId)) return;
  pointers.set(e.pointerId, [e.offsetX, e.offsetY]);
  if (pinch && pointers.size >= 2) {
    const [p, q] = [...pointers.values()];
    const d = Math.hypot(p[0] - q[0], p[1] - q[1]) || 1, cx = (p[0] + q[0]) / 2, cy = (p[1] + q[1]) / 2;
    const bl = Math.min(140, Math.max(8, pinch.bl * d / pinch.d));
    view.bl = bl;
    view.ox = cx - (pinch.cx - pinch.ox) * bl / pinch.bl;
    view.oy = cy - (pinch.cy - pinch.oy) * bl / pinch.bl;
    schedule();
    return;
  }
  const pr = state.press;
  if (!pr) return;
  if (!pr.moved && Math.hypot(e.offsetX - pr.x, e.offsetY - pr.y) < DRAG_PX) return;
  pr.moved = true;
  clearTimeout(longTimer);
  const [mx, my] = toModel(e.offsetX, e.offsetY);
  dragTool(pr, mx, my);
  schedule();
}

function bondEnd(pr, mx, my) {
  const from = pr.item && pr.item[0] === 'atom' ? doc.atoms.get(pr.item[1]) : { x: pr.mx, y: pr.my };
  const exclude = pr.item && pr.item[0] === 'atom' ? [pr.item[1]] : [];
  const hit = doc.findAtom(mx, my, 0.4, exclude);
  if (hit !== null) { const a = doc.atoms.get(hit); return [from, a.x, a.y]; }
  const ang = snapAngle(Math.atan2(my - from.y, mx - from.x), 30);
  return [from, from.x + Math.cos(ang), from.y + Math.sin(ang)];
}

function dragTool(pr, mx, my) {
  switch (state.tool) {
    case 'select': {
      if (pr.moving === undefined) {
        const onSel = pr.item && state.sel.has(key(...normalItem(pr.item)));
        pr.moving = !!(onSel || (pr.inBox && !pr.item));
        if (pr.item && !onSel) { state.sel = selectMolecule(pr.item); pr.moving = true; renderSelBar(); }
        if (pr.moving) { snapshot(); pr.lastX = pr.mx; pr.lastY = pr.my; }
      }
      if (pr.moving) { moveSelection(mx - pr.lastX, my - pr.lastY); pr.lastX = mx; pr.lastY = my; }
      else state.preview = { band: [Math.min(pr.mx, mx), Math.min(pr.my, my), Math.max(pr.mx, mx), Math.max(pr.my, my)] };
      break;
    }
    case 'bond': {
      const [from, x, y] = bondEnd(pr, mx, my);
      state.preview = { points: [[from.x, from.y], [x, y]], end: [x, y] };
      break;
    }
    case 'chain': {
      const from = pr.item && pr.item[0] === 'atom' ? doc.atoms.get(pr.item[1]) : { x: pr.mx, y: pr.my };
      state.preview = { points: [[from.x, from.y], ...doc.chainPoints(from.x, from.y, mx, my)] };
      break;
    }
    case 'arrow':
      state.preview = { points: [[pr.mx, pr.my], [mx, my]] };
      break;
    case 'erase': {
      const it = itemAt(mx, my);
      if (it) { if (!pr.snapped) { snapshot(); pr.snapped = true; } eraseItem(it, true); edited(); }
      break;
    }
    default:
      break;
  }
}

function moveSelection(dx, dy) {
  for (const k of state.sel) {
    const [kind, id] = parseKey(k);
    if (kind === 'atom') { const a = doc.atoms.get(id); if (a) { a.x += dx; a.y += dy; } }
    if (kind === 'arrow') { const a = doc.arrows.get(id); if (a) { a.x1 += dx; a.y1 += dy; a.x2 += dx; a.y2 += dy; } }
    if (kind === 'text') { const t = doc.texts.get(id); if (t) { t.x += dx; t.y += dy; } }
  }
}

function onUp(e) {
  if (!pointers.has(e.pointerId)) return;
  pointers.delete(e.pointerId);
  if (pinch) { if (pointers.size < 2) pinch = null; return; }
  clearTimeout(longTimer);
  const pr = state.press;
  state.press = null; state.hint = null;
  if (!pr) { schedule(); return; }
  const pv = state.preview;
  state.preview = null;
  if (e.type === 'pointercancel') { schedule(); return; }
  const [mx, my] = toModel(e.offsetX, e.offsetY);
  if (pr.moved) finishDrag(pr, pv, mx, my); else tap(pr);
  schedule();
}

function finishDrag(pr, pv, mx, my) {
  switch (state.tool) {
    case 'select':
      if (pr.moving) { edited(); break; }
      if (pv && pv.band) {
        const [x0, y0, x1, y1] = pv.band, sel = new Set();
        for (const [i, a] of doc.atoms) if (a.x >= x0 && a.x <= x1 && a.y >= y0 && a.y <= y1) sel.add(key('atom', i));
        for (const [i, a] of doc.arrows) if (Math.min(a.x1, a.x2) >= x0 && Math.max(a.x1, a.x2) <= x1 && Math.min(a.y1, a.y2) >= y0 && Math.max(a.y1, a.y2) <= y1) sel.add(key('arrow', i));
        for (const [i, t] of doc.texts) if (t.x >= x0 && t.x <= x1 && t.y >= y0 && t.y <= y1) sel.add(key('text', i));
        state.sel = sel; renderSelBar(); updateFormula();
      }
      break;
    case 'bond': {
      if (!pv) break;
      snapshot();
      const a = pr.item && pr.item[0] === 'atom' ? pr.item[1] : doc.atomAtOrNew(pr.mx, pr.my);
      const b = doc.atomAtOrNew(pv.end[0], pv.end[1], 'C', [a]);
      makeBond(a, b);
      edited();
      break;
    }
    case 'chain': {
      if (!pv) break;
      snapshot();
      let prev = pr.item && pr.item[0] === 'atom' ? pr.item[1] : doc.atomAtOrNew(pr.mx, pr.my);
      for (const [x, y] of pv.points.slice(1)) { const n = doc.atomAtOrNew(x, y); doc.addBond(prev, n); prev = n; }
      edited();
      break;
    }
    case 'arrow':
      if (Math.hypot(mx - pr.mx, my - pr.my) > 0.4) { snapshot(); doc.addArrow(pr.mx, pr.my, mx, my, state.arrow || 'forward'); edited(); }
      break;
    default:
      break;
  }
}

function makeBond(a, b) {
  const kind = state.bond, id = doc.addBond(a, b);
  if (id === null) return null;
  applyBond(id, kind, a);
  return id;
}

function applyBond(id, kind, narrow = null) {
  const bd = doc.bonds.get(id);
  bd.stereo = ''; bd.style = ''; bd.side = null;
  if (kind === 'double') bd.order = 2;
  else if (kind === 'triple') bd.order = 3;
  else {
    bd.order = 1;
    if (kind === 'wedge' || kind === 'hash') {
      bd.stereo = kind;
      if (narrow !== null && bd.a !== narrow) { const t = bd.a; bd.a = bd.b; bd.b = t; }
    } else if (kind !== 'single') bd.style = kind;
  }
}

function tap(pr) {
  const it = pr.item;
  switch (state.tool) {
    case 'select':
      if (it) {
        const k = key(...normalItem(it));
        if (state.sel.has(k)) state.sel = new Set();    // a second tap lets go of it
        else state.sel = selectMolecule(it);
      } else if (!pr.inBox) state.sel = new Set();
      renderSelBar(); updateFormula();
      break;
    case 'bond': {
      snapshot();
      if (it && it[0] === 'bond') tapBond(it[1]);
      else {
        const a = it && it[0] === 'atom' ? it[1] : doc.addAtom(pr.mx, pr.my);
        const ang = doc.nextBondAngle(a), at = doc.atoms.get(a);
        const b = doc.atomAtOrNew(at.x + Math.cos(ang), at.y + Math.sin(ang), 'C', [a]);
        makeBond(a, b);
      }
      edited();
      break;
    }
    case 'ring': {
      snapshot();
      const aromatic = state.ring === 'benzene', n = aromatic ? 6 : Number(state.ring[1]);
      if (it && it[0] === 'bond') doc.addRing(doc.ringOnBond(it[1], n), aromatic);
      else if (it && it[0] === 'atom') doc.addRing(doc.ringOnAtom(it[1], n), aromatic);
      else doc.addRing(doc.ringPoints(n, pr.mx, pr.my, n % 2 || n === 6 ? -Math.PI / 2 : -Math.PI / 2 + Math.PI / n), aromatic);
      edited();
      break;
    }
    case 'atom':
      snapshot();
      if (it && it[0] === 'atom') {
        const at = doc.atoms.get(it[1]);
        if (at.el === state.el && !at.label) doc.cycleHydrogens(it[1]);   // tapped again: NH2 -> NH -> N -> NH2
        else setElement(it[1], state.el);
      }
      else if (!it) setElement(doc.addAtom(pr.mx, pr.my), state.el);
      edited();
      break;
    case 'chain':
      if (it && it[0] === 'atom') {
        snapshot();
        const a = doc.atoms.get(it[1]), ang = doc.nextBondAngle(it[1]);
        doc.addBond(it[1], doc.atomAtOrNew(a.x + Math.cos(ang), a.y + Math.sin(ang), 'C', [it[1]]));
        edited();
      } else toast('Drag to draw a chain');
      break;
    case 'charge':
      if (it && it[0] === 'atom') { snapshot(); applyCharge(doc.atoms.get(it[1]), state.sign || '1'); edited(); }
      break;
    case 'erase':
      if (it) { snapshot(); eraseItem(it, false); edited(); }
      break;
    case 'arrow':
      if (it && it[0] === 'arrow') {
        snapshot();
        const ar = doc.arrows.get(it[1]), kinds = ARROWS.map(a => a[0]);
        ar.kind = kinds[(kinds.indexOf(ar.kind) + 1) % kinds.length];
        edited();
      } else toast('Drag to draw an arrow');
      break;
    case 'text':
      if (it && it[0] === 'text') textDialog(it[1], pr.mx, pr.my);
      else textDialog(null, pr.mx, pr.my);
      break;
    default:
      break;
  }
}

/** Tapping a bond with the bond tool: single cycles the order; the others put that kind on it, again flips or clears it. */
function tapBond(id) {
  const bd = doc.bonds.get(id), kind = state.bond;
  if (kind === 'single') {
    if (bd.stereo || bd.style) applyBond(id, 'single');
    else bd.order = bd.order === 3 ? 1 : bd.order + 1;
  } else if ((kind === 'wedge' || kind === 'hash') && bd.stereo === kind) { const t = bd.a; bd.a = bd.b; bd.b = t; }
  else if (kind === 'double' && bd.order === 2) {         // the second line: left, right, centred, automatic
    const order = [null, 1, -1, 0];
    bd.side = order[(order.indexOf(bd.side ?? null) + 1) % order.length];
  } else applyBond(id, kind);
}

function setElement(a, text) {
  const parsed = parseLabel(text);
  const at = doc.atoms.get(a);
  if (!parsed) return false;
  if (parsed[0] === 'nickname') setNickname(RD, doc, a, parsed[1]);
  else {
    at.el = parsed[1]; at.charge = parsed[2]; at.label = '';
    const [hs, radical] = labelExtras(text);                // CH2, CH3•, O•+: typed H and unpaired electrons
    at.radical = radical; at.hcount = null;
    if (hs !== null && hs !== doc.implicitH(a)) at.hcount = hs;
  }
  return true;
}

function eraseItem(it, dragging) {
  if (it[0] === 'atom') {
    if (!doc.atoms.has(it[1])) return;
    if (dragging || !doc.stripAtom(it[1])) doc.removeAtom(it[1]);
  } else if (it[0] === 'bond') doc.removeBond(it[1]);
  else if (it[0] === 'arrow') doc.arrows.delete(it[1]);
  else if (it[0] === 'text') doc.texts.delete(it[1]);
}

// ---------------------------------------------------------------------------------------------- the selection's bar
function renderSelBar() {
  const bar = $('selBar');
  if (state.tool !== 'select' || !state.sel.size) { bar.hidden = true; return; }
  bar.hidden = false;
  const atoms = selectedAtoms().size > 0;
  bar.innerHTML = [['del', 'Delete'], ['dup', 'Copy'], ...(atoms ? [['clean', 'Clean up'], ['smiles', 'SMILES']] : []), ['flip', 'Flip'], ['ccw', '⟲'], ['cw', '⟳']]
    .map(([k, t]) => `<button data-sel="${k}">${t}</button>`).join('');
}

function transformSelection(fn) {
  for (const k of state.sel) {
    const [kind, id] = parseKey(k);
    if (kind === 'atom') { const a = doc.atoms.get(id); [a.x, a.y] = fn(a.x, a.y); }
    if (kind === 'arrow') { const a = doc.arrows.get(id); [a.x1, a.y1] = fn(a.x1, a.y1); [a.x2, a.y2] = fn(a.x2, a.y2); }
    if (kind === 'text') { const t = doc.texts.get(id); [t.x, t.y] = fn(t.x, t.y); }
  }
}

function selAction(k) {
  const box = selBox();
  if (!box) return;
  const cx = (box[0] + box[2]) / 2, cy = (box[1] + box[3]) / 2;
  if (k === 'del') {
    snapshot();
    for (const s of [...state.sel]) { const [kind, id] = parseKey(s); eraseItem([kind, id], true); }
    state.sel = new Set();
  } else if (k === 'dup') {
    snapshot();
    const ids = selectedAtoms(), map = new Map(), sel = new Set();
    for (const a of ids) { const at = doc.atoms.get(a), n = doc.addAtom(at.x + 1, at.y + 1, at.el, at.charge); doc.atoms.get(n).label = at.label; map.set(a, n); sel.add(key('atom', n)); }
    for (const b of [...doc.bonds.values()]) if (ids.has(b.a) && ids.has(b.b)) { const i = doc.addBond(map.get(b.a), map.get(b.b), b.order, b.stereo); Object.assign(doc.bonds.get(i), { side: b.side, style: b.style }); }
    state.sel = sel;
  } else if (k === 'flip') {
    // a mirror image: the wedges change over so that the molecule stays the same, as the desktop's Rotate 180
    snapshot();
    transformSelection((x, y) => [2 * cx - x, y]);
    const ids = selectedAtoms();
    for (const b of doc.bonds.values()) if (ids.has(b.a) && ids.has(b.b) && b.stereo) b.stereo = b.stereo === 'wedge' ? 'hash' : 'wedge';
  } else if (k === 'cw' || k === 'ccw') {
    snapshot();
    const t = (k === 'cw' ? 15 : -15) * Math.PI / 180, c = Math.cos(t), s = Math.sin(t);
    transformSelection((x, y) => [cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c]);
  } else if (k === 'clean') cleanUp();
  else if (k === 'smiles') {
    const info = analyse(RD, doc, selectedAtoms());
    copy(info.smiles, 'SMILES copied');
    return;
  }
  edited();
}

function cleanUp() {
  const ids = selectedAtoms();
  const frags = doc.fragments(ids.size ? ids : null).map(f => doc.fragmentOf(f[0]));
  const seen = new Set();
  snapshot();
  let skipped = 0;
  const sel = new Set();
  for (const frag of frags) {
    if (seen.has(frag[0])) continue;
    seen.add(frag[0]);
    if (frag.some(a => doc.atoms.get(a).label)) { skipped++; continue; }
    const [block] = doc.toMolblock(new Set(frag));
    const mol = RD.get_mol(block);
    if (!mol) { skipped++; continue; }
    mol.remove_hs_in_place();
    const smiles = mol.get_smiles();
    mol.delete();
    const layout = layoutSmiles(RD, smiles);
    if (!layout || !layout.stereoOk) { skipped++; continue; }
    const cx = frag.reduce((s, a) => s + doc.atoms.get(a).x, 0) / frag.length, cy = frag.reduce((s, a) => s + doc.atoms.get(a).y, 0) / frag.length;
    for (const a of frag) doc.removeAtom(a);
    for (const n of insertLayout(doc, layout, cx, cy)) sel.add(key('atom', n));
  }
  if (ids.size) state.sel = sel;
  if (skipped) toast(`${skipped} structure${skipped > 1 ? 's' : ''} left as drawn (labels such as OMe, or stereo that could not be kept)`);
}

// ---------------------------------------------------------------------------------------------- sheets
function openSheet(html, onClose = null) {
  const s = $('sheet');
  s.innerHTML = html; s.hidden = false; s.scrollTop = 0; $('scrim').hidden = false;
  sheetOnClose = onClose;
  return s;
}
function closeSheet() {
  if ($('sheet').hidden) return false;
  $('sheet').hidden = true; $('scrim').hidden = true;
  const f = sheetOnClose; sheetOnClose = null;
  if (f) f();
  return true;
}
const esc = t => String(t).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function ask(title, { value = '', placeholder = '', hint = '', ok = 'OK', extra = '' }, done) {
  const s = openSheet(`<h2>${esc(title)}</h2><input class="field" id="askField" value="${esc(value)}" placeholder="${esc(placeholder)}" autocapitalize="off" autocomplete="off" spellcheck="false">
    ${hint ? `<div class="hint">${esc(hint)}</div>` : ''}${extra}
    <div class="actions"><button class="btn" id="askCancel">Cancel</button><button class="btn primary" id="askOk">${esc(ok)}</button></div>`);
  const field = s.querySelector('#askField');
  const submit = () => { const v = field.value.trim(); closeSheet(); if (v) done(v); };
  s.querySelector('#askCancel').onclick = closeSheet;
  s.querySelector('#askOk').onclick = submit;
  field.onkeydown = e => { if (e.key === 'Enter') submit(); };
  field.focus();
  return s;
}

function labelDialog(a) {
  const at = doc.atoms.get(a);
  const groups = QUICK_GROUPS.filter(g => NICKNAMES[g] || GROUPS[g] || parseLabel(g));
  const s = ask('Atom label', { value: at.label || (at.el === 'C' ? '' : at.el), placeholder: 'N, OH, Na+, OMe, Boc …', ok: 'Set',
    extra: `<div class="grid">${groups.map(g => `<button data-g="${g}">${g}</button>`).join('')}</div>` }, text => {
    snapshot();
    if (!setElement(a, text)) { undo.pop(); toast(`"${text}" is not an element or a known group`); return; }
    edited();
  });
  s.querySelectorAll('[data-g]').forEach(b => { b.onclick = () => {
    const g = b.dataset.g;
    closeSheet();
    snapshot();
    if (!setElement(a, g)) { undo.pop(); return; }
    edited();
  }; });
}

function elementDialog() {
  ask('Atom or group', { value: state.el, placeholder: 'Se, Zn, OMe, Boc …', ok: 'Use', hint: 'Then tap an atom. Long-press any atom to label it.' }, text => {
    if (!parseLabel(text)) { toast(`"${text}" is not an element or a known group`); return; }
    state.el = text; renderTools();
  });
}

function textDialog(id, x, y) {
  const t = id !== null ? doc.texts.get(id) : null;
  ask(t ? 'Edit caption' : 'Add caption', { value: t ? t.text : '', ok: t ? 'Save' : 'Add' }, text => {
    snapshot();
    if (t) t.text = text; else doc.addText(x, y, text);
    edited();
  });
}

function showInfo() {
  if (!doc.atoms.size) { toast('Draw a structure first'); return; }
  const ids = selectedAtoms();
  const info = analyse(RD, doc, ids.size ? ids : null);
  const row = (a, b, copyable = '') => `<tr class="${copyable ? 'copy' : ''}" ${copyable ? `data-copy="${esc(copyable)}"` : ''}><td>${a}</td><td>${b}</td></tr>`;
  const pos = info.ions.filter(i => i.z > 0), neg = info.ions.filter(i => i.z < 0);
  const ions = list => list.map(i => row(esc(i.label), i.mz.toFixed(4), i.mz.toFixed(4))).join('');
  const s = openSheet(`<h2>${ids.size ? 'Selection' : info.fragments > 1 ? `Everything drawn (${info.fragments} structures)` : 'Structure'}</h2>
    <div class="big">${info.formulaHtml}</div>
    ${info.errors.map(e => `<div class="problem">${esc(e)}</div>`).join('')}
    <table>${row('Exact mass', info.exact.toFixed(4), info.exact.toFixed(4))}${row('Molecular weight', info.average.toFixed(3), info.average.toFixed(3))}</table>
    ${info.analysis ? `<div class="hint">${esc(info.analysis)}</div>` : ''}
    ${pos.length ? `<h3>Positive ions (m/z)</h3><table>${ions(pos)}</table>` : ''}
    ${neg.length ? `<h3>Negative ions (m/z)</h3><table>${ions(neg)}</table>` : ''}
    <h3>SMILES</h3><div class="mono copy" data-copy="${esc(info.smiles)}">${esc(info.smiles || '—')}</div>
    <h3>InChIKey</h3><div class="mono copy" data-copy="${esc(info.inchikey)}">${esc(info.inchikey || '—')}</div>
    <div class="hint">Tap a value to copy it.</div>
    <div class="actions"><button class="btn primary" id="infoClose">Close</button></div>`);
  s.querySelector('#infoClose').onclick = closeSheet;
  s.querySelectorAll('[data-copy]').forEach(el => { el.onclick = () => copy(el.dataset.copy, 'Copied'); });
}

async function copy(text, done) {
  if (!text) return;
  try { await navigator.clipboard.writeText(text); toast(done); }
  catch { toast('Could not reach the clipboard'); }
}

function showMenu() {
  const items = [['name', 'Structure from a name', 'Looked up online at PubChem: only the name is sent'], ['smiles', 'Paste SMILES or MOL', ''], ['templates', 'Templates', 'Amino acids, sugars, rings, reagents …'],
    ['clean', 'Clean up everything', 'Usual orientation, even bonds'], ['fit', 'Fit to screen', ''],
    ['open', 'Open a file', '.sed from the desktop, .mol, .sdf, .smi'], ['save', 'Save or send the drawing (.sed)', ''], ['mol', 'Send as MOL', ''], ['png', 'Send as a picture (PNG)', ''], ['copySmiles', 'Copy SMILES', ''],
    ['new', 'New drawing', ''], ['about', 'About', '']];
  const s = openSheet(`<h2>Kekulex</h2>${items.map(([k, t, sub]) => `<button class="item ${k === 'new' ? 'danger' : ''}" data-m="${k}">${t}${sub ? `<small>${sub}</small>` : ''}</button>`).join('')}`);
  s.querySelectorAll('[data-m]').forEach(b => { b.onclick = () => { const k = b.dataset.m; closeSheet(); menuAction(k); }; });
}

function menuAction(k) {
  if (k === 'name') nameDialog();
  else if (k === 'smiles') ask('Paste a structure', { placeholder: 'SMILES or MOL text', ok: 'Add' }, addText);
  else if (k === 'templates') templatesSheet();
  else if (k === 'clean') { if (doc.atoms.size) { cleanUp(); edited(); fit(); } }
  else if (k === 'fit') fit();
  else if (k === 'open') $('fileInput').click();
  else if (k === 'save') shareFile(`${path}.sed`, JSON.stringify(doc.toDict(), null, 1), 'application/json');
  else if (k === 'mol') {
    if (doc.atoms.size) shareFile(`${path}.mol`, doc.toMolblock(null, path)[0], 'chemical/x-mdl-molfile');
  } else if (k === 'png') sharePicture();
  else if (k === 'copySmiles') { const info = analyse(RD, doc); copy(info.smiles, 'SMILES copied'); }
  else if (k === 'new') {
    if (!doc.isEmpty()) snapshot();
    doc = new Doc(); state.sel.clear(); path = 'structure'; edited(); fit();
    if (undo.length) toast('New drawing. Undo brings the last one back.');
  } else if (k === 'about') {
    const s = openSheet(`<h2>Kekulex</h2><p>Draw chemical structures and read off formula, exact mass and LC-MS adducts. Drawings open on the desktop program and back.</p>
      <p class="hint">Chemistry by RDKit (BSD licence). Name lookup by PubChem. On the desktop program only: NMR prediction, 3D, MD export, .cdxml, stoichiometry, TLC and bio objects (kept when you save a drawing here).</p>
      <div class="actions"><button class="btn primary" id="aboutClose">Close</button></div>`);
    s.querySelector('#aboutClose').onclick = closeSheet;
  }
}

function placeCentre() {
  const r = canvas.getBoundingClientRect();
  const b = bounds(doc);
  if (!b) return toModel(r.width / 2, r.height / 2);
  return [b[2] + 4, (b[1] + b[3]) / 2];
}

function addSmiles(smiles) {
  const [cx, cy] = placeCentre();
  snapshot();
  const r = insertSmiles(RD, doc, smiles, cx, cy);
  if (!r) { undo.pop(); return null; }
  state.tool === 'select' && (state.sel = new Set(r.ids.map(i => key('atom', i))));
  if (!r.stereoOk) toast('Check the stereo of this drawing: it could not be confirmed');
  edited(); fit();
  return r;
}

function addText(text) {
  if (text.includes('M  END')) {
    const [cx, cy] = placeCentre();
    snapshot();
    try { insertMolblock(doc, text, cx, cy); edited(); fit(); } catch (e) { undo.pop(); toast(e.message); }
    return;
  }
  if (!addSmiles(text)) toast('This is not a SMILES or MOL structure');
}

function nameDialog() {
  ask('Structure from a name', { placeholder: 'caffeine, betanin, lupeol …', ok: 'Look up', hint: 'Looked up online at PubChem: only the name is sent.' }, async name => {
    toast('Looking up…', 15000);
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/property/SMILES,IsomericSMILES/JSON`;
    try {
      const ctl = new AbortController(), timer = setTimeout(() => ctl.abort(), 15000);
      const res = await fetch(url, { signal: ctl.signal });
      clearTimeout(timer);
      if (res.status === 404) { toast(`PubChem has no compound called "${name}"`); return; }
      if (!res.ok) { toast(`PubChem answered ${res.status}`); return; }
      const p = (await res.json()).PropertyTable?.Properties?.[0];
      const smiles = p && (p.SMILES || p.IsomericSMILES || p.CanonicalSMILES);
      if (!smiles) { toast(`PubChem has no compound called "${name}"`); return; }
      $('toast').hidden = true;
      if (!addSmiles(usualForm(RD, smiles))) toast('PubChem\'s structure could not be read');
    } catch { toast('PubChem could not be reached. Are you online?'); }
  });
}

function templatesSheet() {
  const cats = Object.keys(TEMPLATES);
  const s = openSheet(`<h2>Templates</h2>${cats.map((c, i) => `<button class="item" data-c="${i}">${esc(c)}<small>${TEMPLATES[c].length} structures</small></button>`).join('')}`);
  s.querySelectorAll('[data-c]').forEach(b => { b.onclick = () => {
    const items = TEMPLATES[cats[Number(b.dataset.c)]];
    const t = openSheet(`<h2>${esc(cats[Number(b.dataset.c)])}</h2>${items.map(([n], i) => `<button class="item" data-t="${i}">${esc(n)}</button>`).join('')}`);
    t.querySelectorAll('[data-t]').forEach(x => { x.onclick = () => { const sm = items[Number(x.dataset.t)][1]; closeSheet(); addSmiles(sm); }; });
  }; });
}

// ---------------------------------------------------------------------------------------------- files
const plugins = () => window.Capacitor?.Plugins || {};

async function shareFile(name, text, type) {
  if (doc.isEmpty()) { toast('There is nothing to send yet'); return; }
  const { Filesystem, Share } = plugins();
  if (Filesystem && Share) {
    try {
      const { uri } = await Filesystem.writeFile({ path: name, data: text, directory: 'CACHE', encoding: 'utf8' });
      await Share.share({ title: name, files: [uri], dialogTitle: `Send ${name}` });
    } catch (e) { if (!/cancel/i.test(String(e?.message))) toast(`Could not send the file: ${e?.message || e}`); }
    return;
  }
  download(new Blob([text], { type }), name);
}

async function sharePicture() {
  const c = toCanvas(doc, 80);
  if (!c) { toast('There is nothing to send yet'); return; }
  const name = `${path}.png`, { Filesystem, Share } = plugins();
  if (Filesystem && Share) {
    try {
      const data = c.toDataURL('image/png').split(',')[1];
      const { uri } = await Filesystem.writeFile({ path: name, data, directory: 'CACHE' });
      await Share.share({ title: name, files: [uri], dialogTitle: 'Send the picture' });
    } catch (e) { if (!/cancel/i.test(String(e?.message))) toast(`Could not send the picture: ${e?.message || e}`); }
    return;
  }
  c.toBlob(b => download(b, name), 'image/png');
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}

async function openFile(file) {
  const text = await file.text(), low = file.name.toLowerCase();
  const base = file.name.replace(/\.[^.]+$/, '');
  try {
    if (low.endsWith('.sed') || (low.endsWith('.json') && text.includes('"atoms"'))) {
      const d = JSON.parse(text);
      if (!Array.isArray(d.atoms)) throw new Error('This is not a Kekulex drawing.');
      if (!doc.isEmpty()) snapshot();
      doc = Doc.fromDict(d); path = base; state.sel.clear(); edited(); fit();
      const hidden = (doc.other.shapes || []).filter(s => !DRAWN_SHAPES.includes(s.kind)).length + (doc.other.grids || []).length + (doc.other.seqs || []).length + (doc.other.tlcs || []).length;
      if (hidden) toast(`${hidden} object${hidden > 1 ? 's' : ''} (bio, TLC, grids …) can only be seen on the desktop. They are kept when you save.`, 5000);
    } else if (low.endsWith('.mol') || low.endsWith('.sdf') || low.endsWith('.sd')) {
      const blocks = text.split('$$$$').filter(b => b.includes('M  END'));
      if (!blocks.length) throw new Error('No structure found in this file.');
      snapshot();
      let [cx, cy] = placeCentre();
      for (const b of blocks.slice(0, 20)) { insertMolblock(doc, b.replace(/^\s*\n/, m => m), cx, cy); cx += 6; }
      edited(); fit();
      if (blocks.length > 20) toast(`The first 20 of ${blocks.length} structures`);
    } else if (low.endsWith('.smi') || low.endsWith('.smiles') || low.endsWith('.txt')) {
      const lines = text.split(/\r?\n/).map(l => l.trim().split(/\s+/)[0]).filter(Boolean).slice(0, 20);
      let added = 0;
      for (const l of lines) if (addSmiles(l)) added++;
      if (!added) throw new Error('No SMILES found in this file.');
    } else if (low.endsWith('.cdxml') || low.endsWith('.cdx')) throw new Error('ChemDraw files open on the desktop program: save them there as .sed.');
    else throw new Error('This phone version opens .sed, .mol, .sdf and .smi files.');
  } catch (e) { toast(e.message, 5000); }
}

// ---------------------------------------------------------------------------------------------- wiring
function wire() {
  document.addEventListener('contextmenu', e => e.preventDefault());   // a long press is ours, not the browser's
  canvas.addEventListener('pointerdown', onDown);
  canvas.addEventListener('pointermove', onMove);
  canvas.addEventListener('pointerup', onUp);
  canvas.addEventListener('pointercancel', onUp);
  canvas.addEventListener('wheel', e => {                  // a mouse, trying it on a PC
    e.preventDefault();
    const f = Math.exp(-e.deltaY / 400), bl = Math.min(140, Math.max(8, view.bl * f));
    view.ox = e.offsetX - (e.offsetX - view.ox) * bl / view.bl;
    view.oy = e.offsetY - (e.offsetY - view.oy) * bl / view.bl;
    view.bl = bl; schedule();
  }, { passive: false });
  $('toolRow').addEventListener('click', e => { const b = e.target.closest('[data-tool]'); if (b) setTool(b.dataset.tool); });
  $('options').addEventListener('click', e => {
    const b = e.target.closest('button');
    if (!b) return;
    const d = b.dataset;
    if (d.bond) state.bond = d.bond;
    else if (d.ring) state.ring = d.ring;
    else if (d.el) state.el = d.el;
    else if (d.elMore) { elementDialog(); return; }
    else if (d.sign) state.sign = d.sign;
    else if (d.arrow) state.arrow = d.arrow;
    renderTools();
  });
  $('selBar').addEventListener('click', e => { const b = e.target.closest('[data-sel]'); if (b) selAction(b.dataset.sel); });
  $('menuBtn').onclick = showMenu;
  $('formulaChip').onclick = showInfo;
  $('undoBtn').onclick = doUndo;
  $('redoBtn').onclick = doRedo;
  $('scrim').onclick = closeSheet;
  $('fileInput').onchange = e => { const f = e.target.files[0]; e.target.value = ''; if (f) openFile(f); };
  window.addEventListener('resize', resize);
  new ResizeObserver(resize).observe(canvas);
  document.addEventListener('visibilitychange', () => { if (document.hidden) autosave(); });
  plugins().App?.addListener?.('backButton', () => {
    if (closeSheet()) return;
    if (state.sel.size) { state.sel.clear(); renderSelBar(); schedule(); return; }
    plugins().App.exitApp();
  });
}

async function main() {
  const splash = document.createElement('div');
  splash.id = 'splash'; splash.textContent = 'Kekulex';
  document.body.append(splash);
  wire();
  renderTools();
  try {
    RD = await window.initRDKitModule({ locateFile: () => 'rdkit/RDKit_minimal.wasm' });
    RD.prefer_coordgen(true);
    setRD(RD);
  } catch (e) {
    splash.textContent = `The chemistry module did not start: ${e?.message || e}`;
    return;
  }
  try {
    const saved = localStorage.getItem(STORE);
    if (saved) doc = Doc.fromDict(JSON.parse(saved));
  } catch { doc = new Doc(); }
  splash.remove();
  resize();
  fit();
  syncHistory();
  updateFormula();
  window.kekulex = { get doc() { return doc; }, get state() { return state; }, get view() { return view; }, get RD() { return RD; }, addSmiles, setTool, analyse: () => analyse(RD, doc), usualForm: s => usualForm(RD, s), putGroup };
}

main();
