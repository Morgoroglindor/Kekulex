// draw.js - paints the drawing on a canvas, as the desktop's render.py does: labels with their H and charge, bonds cut back
// at labels, double bonds on the ring's side, wedges and hashes, arrows, captions and simple shapes. Proportions are ACS 1996
// (bond 14.4 pt), so a bond of BL px scales everything.
const K = { font: 10 / 14.4, line: 0.6 / 14.4, bold: 2 / 14.4, hash: 2.5 / 14.4, gap: 0.18, margin: 1.6 / 14.4 };
export const FONT = 'Arial, Roboto, "Noto Sans", sans-serif';
const SUB = 0.75;

import { NICKNAMES, VALENCE } from './data.js';

function chargeGlyph(z) { return (Math.abs(z) === 1 ? '' : Math.abs(z)) + (z > 0 ? '+' : '−'); }
/** Radical dots and charge beside an atom: •, ••, +, •+, 2−. */
function signGlyph(at) { return '•'.repeat(at.radical || 0) + (at.charge ? chargeGlyph(at.charge) : ''); }

const TOKEN = /^[A-Z][a-z]*\d*$/;
const tokenOk = t => { const l = t.replace(/\d+$/, ''); return !!VALENCE[l] || !!NICKNAMES[l] || ['Mg', 'Li', 'Na', 'K', 'Zn', 'Cu'].includes(l); };

/** The label's tokens in reading order for an atom bonded on its right (CO2Me -> MeO2C, N(CH3)2 -> (H3C)2N) and the
 * index of the token on the atom (model.label_order). */
function tokens(label, reverse) {
  const sign = /[+\-−]+$/.exec(label);
  if (reverse && sign && sign.index > 0) {                 // NMe3+ -> Me3N+: the charge stays at the end
    const [stem, anchor] = tokens(label.slice(0, sign.index), true);
    if (stem.length > 1) return [[...stem, sign[0]], anchor];
  }
  const toks = label.match(/[A-Z][a-z]*\d*|[^A-Z]+/g) || [label];
  if (reverse && toks.length > 1 && !toks.every(t => TOKEN.test(t) && tokenOk(t))) {   // OTBS -> TBSO: a nickname turns as one
    const names = Object.keys(NICKNAMES).filter(k => /^[A-Z][A-Za-z]+$/.test(k)).sort((a, b) => b.length - a.length);
    const out = [];
    let i = 0, ok = true, nick = false;
    while (i < label.length && ok) {
      const n = names.find(k => label.startsWith(k, i) && !/[a-z]/.test(label[i + k.length] || ''));
      if (n) { const d = /^\d*/.exec(label.slice(i + n.length))[0]; out.push(n + d); i += n.length + d.length; nick = true; continue; }
      const m = /^[A-Z][a-z]*\d*/.exec(label.slice(i));
      if (!m || !tokenOk(m[0])) { ok = false; break; }
      out.push(m[0]); i += m[0].length;
    }
    if (ok && nick && out.length > 1) return [[...out].reverse(), out.length - 1];
  }
  if (reverse && toks.length > 1 && toks.join('') === label && toks.every(t => TOKEN.test(t) && tokenOk(t))) return [[...toks].reverse(), toks.length - 1];
  if (reverse && label.includes('(')) {                     // a bracket turns round as one piece
    const parts = label.match(/[A-Z][a-z]*\d*|\(|\)\d*/g) || [];
    const stack = [[]];
    let ok = parts.join('') === label;
    for (const t of parts) {
      if (!ok) break;
      if (t === '(') stack.push(['(']);
      else if (t.startsWith(')')) { if (stack.length < 2) { ok = false; break; } const g = stack.pop(); g.push(t); stack[stack.length - 1].push(g); }
      else if (TOKEN.test(t) && tokenOk(t)) stack[stack.length - 1].push(t);
      else ok = false;
    }
    if (ok && stack.length === 1 && stack[0].length > 1) {
      const turned = units => [...units].reverse().map(u => (Array.isArray(u) ? ['(', ...turned(u.slice(1, -1)), u[u.length - 1]] : u));
      const flat = units => units.flatMap(u => (Array.isArray(u) ? flat(u) : [u]));
      const out = flat(turned(stack[0]));
      const heavy = out.map((t, k) => [t, k]).filter(([t]) => TOKEN.test(t) && t.replace(/\d+$/, '') !== 'H').map(([, k]) => k);
      if (heavy.length) return [out, heavy[heavy.length - 1]];
    }
  }
  return [toks, 0];
}

/** {pieces [[text, x, baseline, size]], hw, hh, rect} of an atom's label, in px around its centre; null for a plain vertex. */
function labelOf(ctx, doc, a, BL, showCarbons = false) {
  const at = doc.atoms.get(a), nb = doc.neighbors(a);
  const cx = at.x * BL, cy = at.y * BL, fs = K.font * BL, small = fs * SUB, ch = fs * 0.72, margin = K.margin * BL;
  const width = (t, size) => { ctx.font = `${size}px ${FONT}`; return ctx.measureText(t).width; };
  const dxs = nb.map(([n]) => doc.atoms.get(n).x - at.x);
  const hasRight = dxs.some(d => d > 0.3), hasLeft = dxs.some(d => d < -0.3);
  const base = cy + ch / 2;
  const pieces = [];
  if (at.label) {
    const [toks, anchor] = tokens(at.label, hasRight && !hasLeft);
    const parts = toks.map(t => { const letters = t.replace(/\d+$/, ''); return { letters, digits: t.slice(letters.length), w: width(letters, fs), wd: t === letters ? 0 : width(t.slice(letters.length), small) }; });
    const aw = parts[anchor].w;
    let x = cx - aw / 2 - parts.slice(0, anchor).reduce((s, p) => s + p.w + p.wd, 0);
    for (const p of parts) {
      pieces.push([p.letters, x, base, fs]);
      if (p.wd) pieces.push([p.digits, x + p.w, base + 0.28 * fs, small]);
      x += p.w + p.wd;
    }
    if (at.charge || at.radical) pieces.push([signGlyph(at), x + 0.02 * BL, base - ch + 0.15 * fs, small]);
    return { pieces, hw: aw / 2 + margin, hh: ch / 2 + margin };
  }
  if (at.el === 'C' && nb.length && !at.charge && !at.radical && (at.hcount === null || at.hcount === undefined) && !showCarbons) return null;
  const ws = width(at.el, fs), x0 = cx - ws / 2;
  pieces.push([at.el, x0, base, fs]);
  let rightEnd = x0 + ws;
  const hs = doc.implicitH(a);
  if (hs) {
    const count = hs > 1 ? String(hs) : '', wh = width('H', fs), wc = count ? width(count, small) : 0;
    let side;
    if (!nb.length && ['O', 'S', 'Se'].includes(at.el)) side = 'left';
    else if (!hasRight) side = 'right';
    else if (!hasLeft) side = 'left';
    else side = doc.freeDirection(a)[1] > 0 ? 'down' : 'up';
    if (side === 'right') {
      pieces.push(['H', rightEnd, base, fs]);
      if (count) pieces.push([count, rightEnd + wh, base + 0.28 * fs, small]);
      rightEnd += wh + wc;
    } else if (side === 'left') {
      pieces.push(['H', x0 - wh - wc, base, fs]);
      if (count) pieces.push([count, x0 - wc, base + 0.28 * fs, small]);
    } else {
      const hy = side === 'down' ? base + ch + 1.4 * margin + 0.04 * BL : base - ch - 1.4 * margin - 0.04 * BL;
      pieces.push(['H', cx - wh / 2, hy, fs]);
      if (count) pieces.push([count, cx + wh / 2, hy + 0.28 * fs, small]);
    }
  }
  if (at.charge || at.radical) {
    const glyph = signGlyph(at);
    if (at.charge_pos) pieces.push([glyph, cx + at.charge_pos[0] * BL - width(glyph, small) / 2, cy + at.charge_pos[1] * BL + ch * SUB / 2, small]);
    else pieces.push([glyph, rightEnd + 0.02 * BL, base - ch + 0.1 * fs, small]);
  }
  return { pieces, hw: ws / 2 + margin, hh: ch / 2 + margin };
}

function exit(lab, ux, uy) {
  if (!lab) return 0;
  const a = Math.abs(ux) > 1e-6 ? lab.hw / Math.abs(ux) : Infinity, b = Math.abs(uy) > 1e-6 ? lab.hh / Math.abs(uy) : Infinity;
  return Math.min(a, b);
}

function line(ctx, x1, y1, x2, y2) { ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); }

function paintBond(ctx, doc, i, labels, BL, color) {
  const bd = doc.bonds.get(i), p = doc.atoms.get(bd.a), q = doc.atoms.get(bd.b);
  const x1 = p.x * BL, y1 = p.y * BL, x2 = q.x * BL, y2 = q.y * BL, L = Math.hypot(x2 - x1, y2 - y1);
  if (L < 1e-6) return;
  const ux = (x2 - x1) / L, uy = (y2 - y1) / L, nx = -uy, ny = ux;
  const la = labels.get(bd.a), lb = labels.get(bd.b);
  const ta = exit(la, ux, uy), tb = exit(lb, ux, uy);
  if (ta + tb >= L) return;
  const seg = (ox, oy, extraA = 0, extraB = 0) => line(ctx, x1 + ux * (ta + extraA) + ox, y1 + uy * (ta + extraA) + oy, x2 - ux * (tb + extraB) + ox, y2 - uy * (tb + extraB) + oy);
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = K.line * BL; ctx.lineCap = 'round';
  const sp = K.gap * BL;
  if (bd.order === 1 && (bd.stereo === 'wedge' || bd.stereo === 'hash')) {
    const wide = 1.5 * K.bold * BL / 2;
    if (bd.stereo === 'wedge') {
      const w2 = wide * (L - tb) / L, ex = x2 - ux * tb, ey = y2 - uy * tb;
      ctx.beginPath(); ctx.moveTo(x1 + ux * ta, y1 + uy * ta); ctx.lineTo(ex + nx * w2, ey + ny * w2); ctx.lineTo(ex - nx * w2, ey - ny * w2); ctx.closePath(); ctx.fill();
    } else {
      const step = K.hash * BL;
      for (let d = Math.max(ta, step * 0.6); d <= L - tb + 1e-6; d += step) {
        const w = Math.max(wide * d / L, K.line * BL * 0.6), cx = x1 + ux * d, cy = y1 + uy * d;
        line(ctx, cx + nx * w, cy + ny * w, cx - nx * w, cy - ny * w);
      }
    }
    return;
  }
  if (bd.order === 1) {
    if (bd.style === 'bold') { ctx.lineWidth = K.bold * BL; ctx.lineCap = 'butt'; seg(0, 0); return; }
    if (bd.style === 'dashed') { ctx.setLineDash([0.14 * BL, 0.11 * BL]); seg(0, 0); ctx.setLineDash([]); return; }
    if (bd.style === 'wavy') {
      const sx = x1 + ux * ta, sy = y1 + uy * ta, len = L - ta - tb, n = Math.max(2, Math.round(len / (0.22 * BL)));
      ctx.beginPath(); ctx.moveTo(sx, sy);
      for (let k = 1; k <= n; k++) {
        const d = len * k / n, m = len * (k - 0.5) / n, a = (k % 2 ? 1 : -1) * 0.2 * BL;
        ctx.quadraticCurveTo(sx + ux * m + nx * a, sy + uy * m + ny * a, sx + ux * d, sy + uy * d);
      }
      ctx.stroke();
      return;
    }
    seg(0, 0);
    return;
  }
  if (bd.order === 2) {
    const side = doc.doubleSide(i), dashed = bd.style === 'dashed';
    const second = (ox, oy, a, b) => { if (dashed) ctx.setLineDash([0.14 * BL, 0.11 * BL]); seg(ox, oy, a, b); ctx.setLineDash([]); };
    if (side === 0) { seg(nx * sp / 2, ny * sp / 2); second(-nx * sp / 2, -ny * sp / 2, 0, 0); return; }
    seg(0, 0);
    const trim = 0.75 * sp;
    const a = la || doc.neighbors(bd.a).length === 1 ? 0 : trim, b = lb || doc.neighbors(bd.b).length === 1 ? 0 : trim;
    second(nx * sp * side, ny * sp * side, a, b);
    return;
  }
  seg(0, 0); seg(nx * sp, ny * sp); seg(-nx * sp, -ny * sp);
}

function head(ctx, x, y, dx, dy, s, half = 0) {
  const L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L, nx = -uy, ny = ux, hl = 7 * s, hw = 2.6 * s, notch = 5.2 * s;
  ctx.beginPath(); ctx.moveTo(x, y);
  ctx.lineTo(x - ux * hl + nx * hw, y - uy * hl + ny * hw);
  if (half < 0) ctx.lineTo(x - ux * notch, y - uy * notch);
  else { ctx.lineTo(x - ux * notch, y - uy * notch); if (half === 0) ctx.lineTo(x - ux * hl - nx * hw, y - uy * hl - ny * hw); }
  ctx.closePath(); ctx.fill();
}

function paintArrow(ctx, ar, BL, color) {
  const s = BL / 14.4, x1 = ar.x1 * BL, y1 = ar.y1 * BL, x2 = ar.x2 * BL, y2 = ar.y2 * BL, dx = x2 - x1, dy = y2 - y1;
  const L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L, nx = -uy, ny = ux;
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = K.line * BL; ctx.lineCap = 'butt';
  if (ar.kind === 'curved' || ar.kind === 'curved_half') {
    const cx = (x1 + x2) / 2 + dy * ar.bend * 2, cy = (y1 + y2) / 2 - dx * ar.bend * 2;
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.quadraticCurveTo(cx, cy, x2, y2); ctx.stroke();
    head(ctx, x2, y2, x2 - cx, y2 - cy, s, ar.kind === 'curved_half' ? 1 : 0);
  } else if (ar.kind === 'equilibrium') {
    const g = 1.4 * s;
    line(ctx, x1 + nx * g, y1 + ny * g, x2 - ux * 3 * s + nx * g, y2 - uy * 3 * s + ny * g);
    line(ctx, x1 - nx * g, y1 - ny * g, x2 - nx * g, y2 - ny * g);
    head(ctx, x2 + nx * g, y2 + ny * g, dx, dy, s, 1);
    head(ctx, x1 - nx * g, y1 - ny * g, -dx, -dy, s, 1);
  } else if (ar.kind === 'retro') {
    const g = 2 * s, h = 5 * s, ex = x2 - ux * 2 * s, ey = y2 - uy * 2 * s;
    line(ctx, x1 + nx * g, y1 + ny * g, ex + nx * g, ey + ny * g);
    line(ctx, x1 - nx * g, y1 - ny * g, ex - nx * g, ey - ny * g);
    ctx.beginPath(); ctx.moveTo(x2 - ux * h + nx * 2 * h / 1.6, y2 - uy * h + ny * 2 * h / 1.6); ctx.lineTo(x2, y2); ctx.lineTo(x2 - ux * h - nx * 2 * h / 1.6, y2 - uy * h - ny * 2 * h / 1.6); ctx.stroke();
  } else {
    if (ar.kind === 'dashed') ctx.setLineDash([4 * s, 3 * s]);
    line(ctx, x1, y1, x2 - ux * 5 * s, y2 - uy * 5 * s);
    ctx.setLineDash([]);
    head(ctx, x2, y2, dx, dy, s);
    if (ar.kind === 'resonance') head(ctx, x1, y1, -dx, -dy, s);
  }
}

const TEXT_SUB = 0.72;                                     // a text's subscript size, as the desktop's render.SUB

/** [[text, isSubscript]] of a text line, as the desktop's render.runs: with formula on, digits right after a letter,
 *  ')' or ']' are subscripts (H2SO4, Pd(PPh3)4); "12 h", "1H", "100 °C" stay on the line. */
export function runs(line, formula) {
  if (!formula) return [[line, false]];
  return [...line.matchAll(/(?<=[A-Za-z)\]])\d+|[^\d]+|\d+/g)].map(m =>
    [m[0], m.index > 0 && /[A-Za-z)\]]/.test(line[m.index - 1]) && /^\d+$/.test(m[0])]);
}

function paintText(ctx, t, BL, color) {
  const size = (t.size || 10) * BL / 14.4, small = size * TEXT_SUB, lines = String(t.text).split('\n');
  const font = s => `${t.italic ? 'italic ' : ''}${t.bold ? 'bold ' : ''}${s}px ${FONT}`;
  ctx.fillStyle = color; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  ctx.font = font(size);
  const ascent = ctx.measureText('H').fontBoundingBoxAscent || size * 0.9;
  lines.forEach((l, k) => {
    const parts = runs(l, t.formula).map(([txt, sub]) => {
      ctx.font = font(sub ? small : size);
      return [txt, sub, ctx.measureText(txt).width];
    });
    const w = parts.reduce((s, p) => s + p[2], 0);
    let x = t.x * BL - (t.align === 'center' ? w / 2 : t.align === 'right' ? w : 0);
    const base = t.y * BL + ascent + k * size * 1.2;
    for (const [txt, sub, wd] of parts) {
      ctx.font = font(sub ? small : size);
      ctx.fillText(txt, x, base + (sub ? 0.28 * size : 0));
      x += wd;
    }
  });
}

function paintShape(ctx, s, BL, color) {
  const x = Math.min(s.x1, s.x2) * BL, y = Math.min(s.y1, s.y2) * BL, w = Math.abs(s.x2 - s.x1) * BL, h = Math.abs(s.y2 - s.y1) * BL;
  ctx.strokeStyle = s.color || color; ctx.lineWidth = K.line * BL;
  if (s.kind === 'line') line(ctx, s.x1 * BL, s.y1 * BL, s.x2 * BL, s.y2 * BL);
  else if (s.kind === 'rect') ctx.strokeRect(x, y, w, h);
  else if (s.kind === 'roundrect') { ctx.beginPath(); ctx.roundRect(x, y, w, h, 0.2 * BL); ctx.stroke(); }
  else if (s.kind === 'oval') { ctx.beginPath(); ctx.ellipse(x + w / 2, y + h / 2, w / 2, h / 2, 0, 0, 7); ctx.stroke(); }
  else if (s.kind === 'brackets') {
    const t = 0.15 * BL;
    ctx.beginPath(); ctx.moveTo(x + t, y); ctx.lineTo(x, y); ctx.lineTo(x, y + h); ctx.lineTo(x + t, y + h);
    ctx.moveTo(x + w - t, y); ctx.lineTo(x + w, y); ctx.lineTo(x + w, y + h); ctx.lineTo(x + w - t, y + h); ctx.stroke();
  }
}

export const DRAWN_SHAPES = ['line', 'rect', 'roundrect', 'oval', 'brackets'];

/** Paints the drawing; ctx is already scaled and moved so that one model unit is BL px. */
export function paint(ctx, doc, BL, opts = {}) {
  const color = opts.color || '#000000';
  const labels = new Map();
  for (const a of doc.atoms.keys()) { const l = labelOf(ctx, doc, a, BL); if (l) labels.set(a, l); }
  for (const s of doc.other.shapes || []) if (DRAWN_SHAPES.includes(s.kind)) paintShape(ctx, s, BL, color);
  for (const [i, b] of doc.bonds) {
    if (b.highlight) { const p = doc.atoms.get(b.a), q = doc.atoms.get(b.b); ctx.strokeStyle = b.highlight; ctx.lineWidth = 0.35 * BL; ctx.lineCap = 'round'; line(ctx, p.x * BL, p.y * BL, q.x * BL, q.y * BL); }
    paintBond(ctx, doc, i, labels, BL, b.color || color);
  }
  for (const [a, lab] of labels) {
    const at = doc.atoms.get(a);
    ctx.fillStyle = (opts.errors && doc.valenceProblem(a)) ? '#D32F2F' : at.color || color;
    ctx.textBaseline = 'alphabetic';
    for (const [text, x, y, size] of lab.pieces) { ctx.font = `${size}px ${FONT}`; ctx.fillText(text, x, y); }
  }
  if (opts.errors) {
    for (const a of doc.atoms.keys()) {
      if (labels.has(a) || !doc.valenceProblem(a)) continue;
      const at = doc.atoms.get(a);
      ctx.strokeStyle = '#D32F2F'; ctx.lineWidth = K.line * BL; ctx.strokeRect(at.x * BL - 0.18 * BL, at.y * BL - 0.18 * BL, 0.36 * BL, 0.36 * BL);
    }
  }
  for (const ar of doc.arrows.values()) paintArrow(ctx, ar, BL, ar.color || color);
  for (const t of doc.texts.values()) paintText(ctx, t, BL, t.color || color);
}

/** [x0, y0, x1, y1] of what is drawn in bond lengths, or null. */
export function bounds(doc) {
  const xs = [], ys = [];
  for (const a of doc.atoms.values()) { xs.push(a.x); ys.push(a.y); }
  for (const r of doc.arrows.values()) { xs.push(r.x1, r.x2); ys.push(r.y1, r.y2); }
  for (const t of doc.texts.values()) { xs.push(t.x, t.x + String(t.text).length * 0.4); ys.push(t.y, t.y + 0.8); }
  for (const s of doc.other.shapes || []) { xs.push(s.x1, s.x2); ys.push(s.y1, s.y2); }
  return xs.length ? [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)] : null;
}

/** The drawing on a new canvas, on white, at px per bond. */
export function toCanvas(doc, px = 60, margin = 1) {
  const b = bounds(doc);
  if (!b) return null;
  const c = document.createElement('canvas');
  c.width = Math.ceil((b[2] - b[0] + 2 * margin) * px); c.height = Math.ceil((b[3] - b[1] + 2 * margin) * px);
  const ctx = c.getContext('2d');
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, c.width, c.height);
  ctx.translate((margin - b[0]) * px, (margin - b[1]) * px);
  paint(ctx, doc, px);
  return c;
}
