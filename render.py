"""
render.py - draws a Doc with QPainter, in points, the same way on the screen and in every export (PNG, SVG, PDF, the
clipboard), so a figure looks exactly as it does in the editor.

    r = Renderer(doc)
    r.paint(painter)                       # painter already scaled so 1 unit = 1 pt
    r.bounds()                             # QRectF (pt) around everything drawn
    to_image(doc, dpi=600) / save_svg / save_pdf / save_png

The default style is the ACS 1996 document style that journals ask for: 14.4 pt bonds, 0.6 pt lines, 2 pt bold,
2.5 pt hash spacing, 18 % double-bond spacing, 1.6 pt margin around labels, Arial 10 pt. Text is drawn as outlines,
so exported figures look the same on a computer without the font.
"""
import math
import re

from PyQt6.QtCore import QPointF, QRectF, QSize, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QFontMetricsF, QImage, QPainter, QPainterPath, QPen, QPolygonF, QTransform

import model

DEFAULT_STYLE = dict(bond_length=14.4, line_width=0.6, bold_width=2.0, hash_spacing=2.5, bond_spacing=0.18, margin=1.6,
                     font='Arial', font_size=10.0, color='#000000', color_atoms=False, page='A4',
                     show_stereo=False, show_atom_numbers=False)
ACS_STYLE = dict(DEFAULT_STYLE)
PAGES = {'A4': (595.28, 841.89), 'Letter': (612.0, 792.0)}  # pt
ELEMENT_COLORS = {'N': '#2846C8', 'O': '#E01E1E', 'S': '#B89400', 'P': '#E07000', 'F': '#2A9A2A', 'Cl': '#1C9A1C',
                  'Br': '#A0291A', 'I': '#7A1FA0', 'B': '#C07070', 'Si': '#B08000', 'Se': '#B07000'}
OVERLAY = '#2D5FB8'                                         # stereo labels and atom numbers
REF = 100.0                                                 # text is shaped at 100 px and scaled down
SUB = 0.72                                                  # subscript / superscript size
_paths = {}


def style_of(doc):
    s = dict(DEFAULT_STYLE)
    s.update({k: v for k, v in (doc.style or {}).items() if k in DEFAULT_STYLE})
    return s


def _font(family, bold=False, italic=False):
    f = QFont(family)
    f.setPixelSize(int(REF))
    f.setBold(bold); f.setItalic(italic)
    f.setHintingPreference(QFont.HintingPreference.PreferNoHinting)
    return f


def text_path(text, family, size, x, y, bold=False, italic=False):
    """Outline of text with its baseline starting at (x, y), size in pt."""
    key = (text, family, bold, italic)
    if key not in _paths:
        p = QPainterPath()
        p.addText(0, 0, _font(family, bold, italic), text)
        _paths[key] = p
    k = size / REF
    return QTransform(k, 0, 0, k, x, y).map(_paths[key])


class Metrics:
    def __init__(self, family):
        self.family = family
        self.fms = {}

    def fm(self, bold=False, italic=False):
        if (bold, italic) not in self.fms:
            self.fms[(bold, italic)] = QFontMetricsF(_font(self.family, bold, italic))
        return self.fms[(bold, italic)]

    def width(self, text, size, bold=False, italic=False):
        return self.fm(bold, italic).horizontalAdvance(text) * size / REF

    def cap(self, size):
        return self.fm().capHeight() * size / REF

    def ascent(self, size):
        return self.fm().ascent() * size / REF

    def line(self, size):
        return self.fm().lineSpacing() * size / REF


def runs(line, formula):
    """[(text, is_subscript)] of a text line: with formula on, digits right after a letter or ')' are subscripts."""
    if not formula:
        return [(line, False)]
    out = []
    for m in re.finditer(r'(?<=[A-Za-z\)\]])\d+|[^\d]+|\d+', line):
        sub = m.start() > 0 and bool(re.match(r'[A-Za-z\)\]]', line[m.start() - 1])) and m.group().isdigit()
        out.append((m.group(), sub))
    return out


def exit_t(ox, oy, ux, uy, hw, hh):
    """How far along (ux, uy) a line starting at offset (ox, oy) from a box centre leaves the box (half sizes hw, hh)."""
    t = math.inf
    if abs(ux) > 1e-9:
        t = min(t, ((hw if ux > 0 else -hw) - ox) / ux)
    if abs(uy) > 1e-9:
        t = min(t, ((hh if uy > 0 else -hh) - oy) / uy)
    return max(t, 0.0) if t != math.inf else 0.0


def dashes(path, x1, y1, x2, y2, on=2.0, off=1.6):
    L = math.hypot(x2 - x1, y2 - y1)
    if L < 1e-6:
        return
    ux, uy = (x2 - x1) / L, (y2 - y1) / L
    d = 0.0
    while d < L:
        e = min(d + on, L)
        path.moveTo(x1 + ux * d, y1 + uy * d); path.lineTo(x1 + ux * e, y1 + uy * e)
        d = e + off


def wave(path, x1, y1, x2, y2, amp=1.4, period=3.2):
    L = math.hypot(x2 - x1, y2 - y1)
    if L < 1e-6:
        return
    ux, uy = (x2 - x1) / L, (y2 - y1) / L
    nx, ny = -uy, ux
    n = max(2, round(L / period * 2))
    step = L / n
    path.moveTo(x1, y1)
    for k in range(n):                                      # half waves as quadratic curves, alternating sides
        s = 1 if k % 2 == 0 else -1
        cx = x1 + ux * (k + 0.5) * step + nx * amp * 2 * s
        cy = y1 + uy * (k + 0.5) * step + ny * amp * 2 * s
        path.quadTo(cx, cy, x1 + ux * (k + 1) * step, y1 + uy * (k + 1) * step)


class Renderer:
    def __init__(self, doc):
        self.doc = doc
        self.s = style_of(doc)
        self.BL = self.s['bond_length']
        self.m = Metrics(self.s['font'])
        self.labels = {a: self._label(a) for a in doc.atoms}
        self._cip = None

    def pt(self, x, y):
        return QPointF(x * self.BL, y * self.BL)

    def color(self, item, a=None):
        if item.color:
            return QColor(item.color)
        if a is not None and self.s['color_atoms']:
            return QColor(ELEMENT_COLORS.get(self.doc.atoms[a].el, self.s['color']))
        return QColor(self.s['color'])

    # ------------------------------------------------------------------------------------------- atom labels
    def _label(self, a):
        """{'pieces': [(text, x, baseline, size)], 'clip': (hw, hh) or None, 'rect': QRectF} in pt; None: no label."""
        doc, s = self.doc, self.s
        at = doc.atoms[a]
        nb = doc.neighbors(a)
        cx, cy = at.x * self.BL, at.y * self.BL
        fs = s['font_size']
        small = fs * SUB
        m = self.m
        ch = m.cap(fs)
        dxs = [doc.atoms[n].x - at.x for n, _b in nb]
        from_right = any(d > 0.3 for d in dxs) and not any(d < -0.3 for d in dxs)
        if at.label:                                        # nickname: the attaching token sits on the atom
            toks, anchor = model.label_order(at.label, from_right)
            widths = []
            for t in toks:
                letters = t.rstrip('0123456789')
                widths.append((m.width(letters, fs), m.width(t[len(letters):], small) if t != letters else 0.0, letters))
            aw = widths[anchor][0]
            x = cx - aw / 2 - sum(w + wd for w, wd, _l in widths[:anchor])
            base = cy + ch / 2
            pieces = []
            x0 = x
            for t, (w, wd, letters) in zip(toks, widths):
                pieces.append((letters, x, base, fs))
                if wd:
                    pieces.append((t[len(letters):], x + w, base + 0.28 * fs, small))
                x += w + wd
            if at.charge:
                pieces.append((self._charge_text(at.charge), x + 0.3, base - ch + 0.15 * fs, small))
            rect = QRectF(x0, base - ch, x - x0, ch + 0.2 * small)
            return dict(pieces=pieces, clip=(aw / 2 + s['margin'], ch / 2 + s['margin']), rect=rect)
        if at.el == 'C' and nb:
            if not at.charge:
                return None
            fx, fy = doc.free_direction(a)                  # a carbocation / carbanion: the charge beside the vertex
            t = self._charge_text(at.charge)
            w = m.width(t, small)
            px, py = cx + fx * 0.32 * self.BL - w / 2, cy + fy * 0.32 * self.BL + m.cap(small) / 2
            return dict(pieces=[(t, px, py, small)], clip=None, rect=QRectF(px, py - m.cap(small), w, m.cap(small)))
        sym = at.el
        hs = doc.implicit_h(a)
        ws = m.width(sym, fs)
        base = cy + ch / 2
        x0 = cx - ws / 2
        pieces = [(sym, x0, base, fs)]
        right_end = x0 + ws
        if hs:
            hcount = str(hs) if hs > 1 else ''
            wh, wc = m.width('H', fs), m.width(hcount, small) if hcount else 0.0
            right = any(d > 0.3 for d in dxs)
            left = any(d < -0.3 for d in dxs)
            if not nb and sym in ('O', 'S', 'Se'):
                side = 'left'
            elif not right:
                side = 'right'
            elif not left:
                side = 'left'
            else:
                side = 'down' if doc.free_direction(a)[1] > 0 else 'up'
            if side == 'right':
                pieces.append(('H', right_end, base, fs))
                if hcount:
                    pieces.append((hcount, right_end + wh, base + 0.28 * fs, small))
                right_end += wh + wc
            elif side == 'left':
                hx = x0 - wh - wc
                pieces.append(('H', hx, base, fs))
                if hcount:
                    pieces.append((hcount, hx + wh, base + 0.28 * fs, small))
            else:
                hy = base + (ch + 1.4 * s['margin'] + 0.5) if side == 'down' else base - (ch + 1.4 * s['margin'] + 0.5)
                pieces.append(('H', cx - wh / 2, hy, fs))
                if hcount:
                    pieces.append((hcount, cx + wh / 2, hy + 0.28 * fs, small))
        if at.charge:                                       # top right, or top left when the bonds come from the right
            t = self._charge_text(at.charge)
            left_end = min(x for _t, x, _y, _s in pieces)
            if nb and doc.free_direction(a)[0] < -0.3 and not (hs and pieces[1][1] < x0):
                pieces.append((t, left_end - m.width(t, small) - 0.3, base - ch + 0.15 * fs, small))
            else:
                pieces.append((t, right_end + 0.3, base - ch + 0.15 * fs, small))
        rect = QRectF()
        for text, x, y, size in pieces:
            r = QRectF(x, y - m.cap(size), m.width(text, size), m.cap(size) + (0.2 * size if size == small else 0))
            rect = r if rect.isNull() else rect.united(r)
        return dict(pieces=pieces, clip=(ws / 2 + s['margin'], ch / 2 + s['margin']), rect=rect)

    @staticmethod
    def _charge_text(z):
        return ('' if abs(z) == 1 else str(abs(z))) + ('+' if z > 0 else '−')

    def label_rect(self, a):
        """QRectF (pt) of the atom's label, or None."""
        lab = self.labels.get(a)
        return lab['rect'] if lab else None

    def labelled(self, a):
        lab = self.labels.get(a)
        return bool(lab and lab['clip'])

    # ------------------------------------------------------------------------------------------- bonds
    def _end(self, a, ox, oy, ux, uy):
        """Distance along (ux, uy) from atom a (offset ox, oy) to where a line may start (0 = unlabelled atom)."""
        lab = self.labels.get(a)
        if not lab or not lab['clip']:
            return 0.0
        return exit_t(ox, oy, ux, uy, *lab['clip'])

    def _seg(self, a, b, ox=0.0, oy=0.0, trim_a=0.0, trim_b=0.0):
        """End points (pt) of a line between atoms a and b, offset by (ox, oy), stopped short of labels; or None."""
        p, q = self.doc.atoms[a], self.doc.atoms[b]
        x1, y1, x2, y2 = p.x * self.BL, p.y * self.BL, q.x * self.BL, q.y * self.BL
        L = math.hypot(x2 - x1, y2 - y1)
        if L < 1e-6:
            return None
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        t1 = max(self._end(a, ox, oy, ux, uy), trim_a)
        t2 = max(self._end(b, ox, oy, -ux, -uy), trim_b)
        if t1 + t2 >= L - 0.2:
            return None
        return (x1 + ox + ux * t1, y1 + oy + uy * t1, x2 + ox - ux * t2, y2 + oy - uy * t2)

    def _line(self, path, a, b, ox=0.0, oy=0.0, trim_a=0.0, trim_b=0.0, dashed=False):
        seg = self._seg(a, b, ox, oy, trim_a, trim_b)
        if seg is None:
            return
        if dashed:
            dashes(path, *seg)
        else:
            path.moveTo(seg[0], seg[1]); path.lineTo(seg[2], seg[3])

    def bond_paths(self, i):
        """(stroked QPainterPath, filled QPainterPath) of bond i, in pt."""
        doc = self.doc
        bd = doc.bonds[i]
        stroke, fill = QPainterPath(), QPainterPath()
        p, q = doc.atoms[bd.a], doc.atoms[bd.b]
        x1, y1, x2, y2 = p.x * self.BL, p.y * self.BL, q.x * self.BL, q.y * self.BL
        L = math.hypot(x2 - x1, y2 - y1)
        if L < 1e-6:
            return stroke, fill
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        nx, ny = -uy, ux
        sp = self.s['bond_spacing'] * self.BL
        if bd.order == 1 and bd.stereo in ('wedge', 'hash'):
            t1 = self._end(bd.a, 0, 0, ux, uy)
            t2 = self._end(bd.b, 0, 0, -ux, -uy)
            wide = 1.5 * self.s['bold_width'] / 2
            if bd.stereo == 'wedge':
                w2 = wide * (L - t2) / L
                ex, ey = x2 - ux * t2, y2 - uy * t2
                fill.addPolygon(QPolygonF([QPointF(x1 + ux * t1, y1 + uy * t1), QPointF(ex + nx * w2, ey + ny * w2),
                                           QPointF(ex - nx * w2, ey - ny * w2)]))
                fill.closeSubpath()
            else:
                step = self.s['hash_spacing']
                d = max(t1, step * 0.6)
                while d <= L - t2 + 1e-6:
                    w = max(wide * d / L, self.s['line_width'] * 0.6)
                    cx, cy = x1 + ux * d, y1 + uy * d
                    stroke.moveTo(cx + nx * w, cy + ny * w)
                    stroke.lineTo(cx - nx * w, cy - ny * w)
                    d += step
            return stroke, fill
        if bd.order == 1:
            seg = self._seg(bd.a, bd.b)
            if seg is None:
                return stroke, fill
            sx1, sy1, sx2, sy2 = seg
            if bd.style == 'bold':
                w = self.s['bold_width'] / 2
                fill.addPolygon(QPolygonF([QPointF(sx1 + nx * w, sy1 + ny * w), QPointF(sx2 + nx * w, sy2 + ny * w),
                                           QPointF(sx2 - nx * w, sy2 - ny * w), QPointF(sx1 - nx * w, sy1 - ny * w)]))
                fill.closeSubpath()
            elif bd.style == 'dashed':
                dashes(stroke, *seg)
            elif bd.style == 'wavy':
                wave(stroke, *seg)
            elif bd.style == 'hashed':
                w = self.s['bold_width'] / 2
                Ls = math.hypot(sx2 - sx1, sy2 - sy1)
                d = self.s['hash_spacing'] * 0.4
                while d <= Ls:
                    cx, cy = sx1 + ux * d, sy1 + uy * d
                    stroke.moveTo(cx + nx * w, cy + ny * w); stroke.lineTo(cx - nx * w, cy - ny * w)
                    d += self.s['hash_spacing']
            else:
                stroke.moveTo(sx1, sy1); stroke.lineTo(sx2, sy2)
        elif bd.order == 2:
            side = doc.double_side(i)
            dashed = bd.style == 'dashed'
            if side == 0:
                self._line(stroke, bd.a, bd.b, nx * sp / 2, ny * sp / 2)
                self._line(stroke, bd.a, bd.b, -nx * sp / 2, -ny * sp / 2, dashed=dashed)
            else:
                self._line(stroke, bd.a, bd.b)
                trim = 0.75 * sp
                ta = 0.0 if self.labelled(bd.a) or len(doc.neighbors(bd.a)) == 1 else trim
                tb = 0.0 if self.labelled(bd.b) or len(doc.neighbors(bd.b)) == 1 else trim
                self._line(stroke, bd.a, bd.b, nx * sp * side, ny * sp * side, ta, tb, dashed=dashed)
        else:
            self._line(stroke, bd.a, bd.b)
            self._line(stroke, bd.a, bd.b, nx * sp, ny * sp)
            self._line(stroke, bd.a, bd.b, -nx * sp, -ny * sp)
        return stroke, fill

    # ------------------------------------------------------------------------------------------- arrows
    HEAD_L, HEAD_W, NOTCH = 7.0, 2.6, 5.2                   # head length, half width, where its back meets the shaft

    def _head(self, fill, x, y, dx, dy, half=0):
        """A filled head with its tip at (x, y) pointing along (dx, dy) (unit); half: +1 / -1 one barb only."""
        nx, ny = -dy, dx
        hl, hw, notch = self.HEAD_L, self.HEAD_W, self.NOTCH
        left = QPointF(x - dx * hl + nx * hw, y - dy * hl + ny * hw)
        right = QPointF(x - dx * hl - nx * hw, y - dy * hl - ny * hw)
        back = QPointF(x - dx * notch, y - dy * notch)
        pts = [QPointF(x, y), left, back, right] if not half else \
            [QPointF(x, y), left if half > 0 else right, back]
        fill.addPolygon(QPolygonF(pts)); fill.closeSubpath()

    def arrow_paths(self, i):
        ar = self.doc.arrows[i]
        BL = self.BL
        x1, y1, x2, y2 = ar.x1 * BL, ar.y1 * BL, ar.x2 * BL, ar.y2 * BL
        stroke, fill = QPainterPath(), QPainterPath()
        L = math.hypot(x2 - x1, y2 - y1)
        if L < 1e-6:
            return stroke, fill
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        nx, ny = -uy, ux
        notch = self.NOTCH
        k = ar.kind
        if k in ('curved', 'curved_half') and ar.bend:
            cx, cy = model.curve_control(ar)
            cx, cy = cx * BL, cy * BL
            tx, ty = x2 - cx, y2 - cy                       # the arc's direction at the head
            tl = math.hypot(tx, ty) or 1.0
            tx, ty = tx / tl, ty / tl
            stroke.moveTo(x1, y1)
            stroke.quadTo(cx, cy, x2 - tx * notch, y2 - ty * notch)
            side = 1 if ar.bend > 0 else -1                  # a fish-hook's barb on the outside of the curve
            self._head(fill, x2, y2, tx, ty, half=(-side if k == 'curved_half' else 0))
            return stroke, fill
        if k == 'equilibrium':
            off = 1.8
            for sgn, (ax, ay, bx, by) in ((1, (x1, y1, x2, y2)), (-1, (x2, y2, x1, y1))):
                dx, dy = (bx - ax) / L, (by - ay) / L
                ox, oy = -nx * off * sgn, -ny * off * sgn
                sx, sy, ex, ey = ax + ox, ay + oy, bx + ox, by + oy
                stroke.moveTo(sx, sy); stroke.lineTo(ex - dx * notch * 0.5, ey - dy * notch * 0.5)
                self._head(fill, ex, ey, dx, dy, half=-1)
            return stroke, fill
        if k == 'retro':                                    # two shafts and an open chevron
            off, hl, hw = 1.8, 6.0, 5.0
            for s in (1, -1):
                stroke.moveTo(x1 + nx * off * s, y1 + ny * off * s)
                stroke.lineTo(x2 - ux * hl * (1 - off / hw) + nx * off * s, y2 - uy * hl * (1 - off / hw) + ny * off * s)
                stroke.moveTo(x2, y2)
                stroke.lineTo(x2 - ux * hl + nx * hw * s, y2 - uy * hl + ny * hw * s)
            return stroke, fill
        start = (x1 + ux * notch, y1 + uy * notch) if k == 'resonance' else (x1, y1)
        end = (x2 - ux * notch, y2 - uy * notch)
        if k == 'dashed':
            dashes(stroke, *start, *end, on=2.6, off=1.8)
        else:
            stroke.moveTo(*start); stroke.lineTo(*end)
        self._head(fill, x2, y2, ux, uy)
        if k == 'resonance':
            self._head(fill, x1, y1, -ux, -uy)
        return stroke, fill

    # ------------------------------------------------------------------------------------------- text
    def text_layout(self, i):
        """[(text, x, baseline, size)] pieces of text item i (subscripts as their own pieces) and its QRectF."""
        t = self.doc.texts[i]
        fs = t.size or self.s['font_size']
        small = fs * SUB
        m = self.m
        lh = m.line(fs)
        pieces = []
        rect = QRectF()
        for k, line in enumerate(t.text.split('\n')):
            parts = runs(line, t.formula)
            widths = [m.width(txt, small if sub else fs, t.bold, t.italic) for txt, sub in parts]
            w = sum(widths)
            x0 = t.x * self.BL - (w / 2 if t.align == 'center' else w if t.align == 'right' else 0)
            base = t.y * self.BL + m.ascent(fs) + k * lh
            x = x0
            for (txt, sub), wd in zip(parts, widths):
                pieces.append((txt, x, base + (0.28 * fs if sub else 0), small if sub else fs))
                x += wd
            r = QRectF(x0, t.y * self.BL + k * lh, max(w, fs * 0.6), lh)
            rect = r if rect.isNull() else rect.united(r)
        return pieces, rect

    def text_rect(self, i):
        """QRectF (pt) of text item i."""
        return self.text_layout(i)[1]

    # ------------------------------------------------------------------------------------------- shapes
    def shape_paths(self, i):
        sh = self.doc.shapes[i]
        BL = self.BL
        x1, x2 = sorted((sh.x1 * BL, sh.x2 * BL))
        y1, y2 = sorted((sh.y1 * BL, sh.y2 * BL))
        p = QPainterPath()
        w, h = x2 - x1, y2 - y1
        if sh.kind == 'line':
            p.moveTo(sh.x1 * BL, sh.y1 * BL); p.lineTo(sh.x2 * BL, sh.y2 * BL)
        elif sh.kind == 'rect':
            p.addRect(QRectF(x1, y1, w, h))
        elif sh.kind == 'roundrect':
            r = min(w, h) * 0.15
            p.addRoundedRect(QRectF(x1, y1, w, h), r, r)
        elif sh.kind == 'oval':
            p.addEllipse(QRectF(x1, y1, w, h))
        elif sh.kind == 'brackets':
            arm = min(max(h * 0.08, 2.5), 6.0)
            for x, s in ((x1, 1), (x2, -1)):
                p.moveTo(x + arm * s, y1); p.lineTo(x, y1); p.lineTo(x, y2); p.lineTo(x + arm * s, y2)
        elif sh.kind == 'parens':
            bow = min(max(h * 0.12, 3.0), 10.0)
            for x, s in ((x1, 1), (x2, -1)):
                p.moveTo(x + bow * s, y1)
                p.cubicTo(x - bow * 0.15 * s, y1 + h * 0.2, x - bow * 0.15 * s, y2 - h * 0.2, x + bow * s, y2)
        elif sh.kind == 'braces':
            bw = min(max(h * 0.08, 2.5), 7.0)
            ym = (y1 + y2) / 2
            for x, s in ((x1, 1), (x2, -1)):
                p.moveTo(x + bw * s, y1)
                p.cubicTo(x, y1, x + bw * s, ym, x - 0.0, ym)
                p.cubicTo(x + bw * s, ym, x, y2, x + bw * s, y2)
        return p

    # ------------------------------------------------------------------------------------------- painting
    def cip(self):
        if self._cip is None:
            try:
                import chem
                self._cip = chem.cip_labels(self.doc)
            except Exception:
                self._cip = ({}, {})
        return self._cip

    def paint(self, p, errors=False):
        """Draw everything (errors: a red box around atoms over their valence, for the screen)."""
        s = self.s
        doc = self.doc
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)

        def pen(color):
            return QPen(color, s['line_width'], Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
        for i, sh in doc.shapes.items():
            p.setPen(pen(self.color(sh))); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(self.shape_paths(i))
        for i, bd in doc.bonds.items():
            color = self.color(bd)
            stroke, fill = self.bond_paths(i)
            p.setPen(pen(color)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(stroke)
            if not fill.isEmpty():
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color))
                p.drawPath(fill)
        for a, lab in self.labels.items():
            if not lab:
                continue
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(self.color(doc.atoms[a], a)))
            for text, x, y, size in lab['pieces']:
                p.drawPath(text_path(text, s['font'], size, x, y))
        for i, ar in doc.arrows.items():
            color = self.color(ar)
            stroke, fill = self.arrow_paths(i)
            p.setPen(pen(color)); p.setBrush(Qt.BrushStyle.NoBrush); p.drawPath(stroke)
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color)); p.drawPath(fill)
        p.setPen(Qt.PenStyle.NoPen)
        for i, t in doc.texts.items():
            p.setBrush(QBrush(self.color(t)))
            for txt, x, y, size in self.text_layout(i)[0]:
                if txt.strip():
                    p.drawPath(text_path(txt, s['font'], size, x, y, t.bold, t.italic))
        if s['show_stereo'] or s['show_atom_numbers']:
            self.paint_overlays(p)
        if errors:
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.setPen(QPen(QColor('#E0201B'), 0.8))
            for a, at in doc.atoms.items():
                if doc.valence_error(a):
                    r = 0.28 * self.BL
                    p.drawRect(QRectF(at.x * self.BL - r, at.y * self.BL - r, 2 * r, 2 * r))

    def paint_overlays(self, p):
        """(R) / (S) and (E) / (Z) labels and atom numbers, small and blue, beside their atoms and bonds."""
        s, doc, BL = self.s, self.doc, self.BL
        size = s['font_size'] * 0.7
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor(OVERLAY)))
        notes = {}
        if s['show_stereo']:
            atoms, bonds = self.cip()
            for a, code in atoms.items():
                notes.setdefault(('atom', a), []).append(f'({code})')
            for b, code in bonds.items():
                notes[('bond', b)] = [f'({code})']
        if s['show_atom_numbers']:
            for k, a in enumerate(doc.atoms, 1):
                notes.setdefault(('atom', a), []).insert(0, str(k))
        for (kind, i), parts in notes.items():
            text = ' '.join(parts)
            w = self.m.width(text, size)
            if kind == 'atom':
                at = doc.atoms[i]
                fx, fy = doc.free_direction(i)
                x, y = at.x * BL + fx * 0.45 * BL, at.y * BL + fy * 0.45 * BL
            else:
                bd = doc.bonds[i]
                pa, pb = doc.atoms[bd.a], doc.atoms[bd.b]
                L = math.hypot(pb.x - pa.x, pb.y - pa.y) or 1.0
                side = -(doc.double_side(i) or 1)
                x = (pa.x + pb.x) / 2 * BL - (pb.y - pa.y) / L * 0.45 * BL * side
                y = (pa.y + pb.y) / 2 * BL + (pb.x - pa.x) / L * 0.45 * BL * side
            p.drawPath(text_path(text, s['font'], size, x - w / 2, y + self.m.cap(size) / 2))

    def bounds(self):
        """QRectF (pt) around everything drawn (null when the drawing is empty)."""
        r = QRectF()

        def add(rect):
            nonlocal r
            r = rect if r.isNull() else r.united(rect)
        pad = self.s['bold_width']
        for a, at in self.doc.atoms.items():
            x, y = at.x * self.BL, at.y * self.BL
            add(QRectF(x - pad, y - pad, 2 * pad, 2 * pad))
            lab = self.labels.get(a)
            if lab:
                add(lab['rect'])
        for i in self.doc.arrows:
            stroke, fill = self.arrow_paths(i)
            add(stroke.boundingRect().united(fill.boundingRect()).adjusted(-1, -1, 1, 1))
        for i in self.doc.texts:
            add(self.text_rect(i))
        for i in self.doc.shapes:
            add(self.shape_paths(i).boundingRect().adjusted(-1, -1, 1, 1))
        if not r.isNull() and (self.s['show_stereo'] or self.s['show_atom_numbers']):
            r = r.adjusted(-6, -6, 6, 6)
        return r

    def item_rect(self, kind, i):
        """QRectF (pt) of one item (an atom's label or point, a bond, an arrow, a text, a shape)."""
        BL = self.BL
        if kind == 'atom':
            at = self.doc.atoms[i]
            lab = self.label_rect(i)
            pt = QRectF(at.x * BL - 1, at.y * BL - 1, 2, 2)
            return lab.united(pt) if lab is not None else pt
        if kind == 'bond':
            bd = self.doc.bonds[i]
            return self.item_rect('atom', bd.a).united(self.item_rect('atom', bd.b))
        if kind == 'arrow':
            stroke, fill = self.arrow_paths(i)
            return stroke.boundingRect().united(fill.boundingRect())
        if kind == 'text':
            return self.text_rect(i)
        return self.shape_paths(i).boundingRect()


# ------------------------------------------------------------------------------------------------ export
MARGIN_PT = 4.0


def _frame(doc):
    r = Renderer(doc)
    b = r.bounds()
    if b.isNull():
        b = QRectF(0, 0, 10, 10)
    return r, b.adjusted(-MARGIN_PT, -MARGIN_PT, MARGIN_PT, MARGIN_PT)


def to_image(doc, dpi=600, transparent=False):
    r, b = _frame(doc)
    k = dpi / 72.0
    img = QImage(QSize(max(1, math.ceil(b.width() * k)), max(1, math.ceil(b.height() * k))),
                 QImage.Format.Format_ARGB32_Premultiplied)
    img.setDotsPerMeterX(round(dpi / 0.0254)); img.setDotsPerMeterY(round(dpi / 0.0254))
    img.fill(QColor(0, 0, 0, 0) if transparent else QColor('white'))
    p = QPainter(img)
    p.scale(k, k); p.translate(-b.left(), -b.top())
    r.paint(p)
    p.end()
    return img


def save_png(doc, path, dpi=600, transparent=False):
    if not to_image(doc, dpi, transparent).save(path, 'PNG'):
        raise OSError(f'Could not write {path}')


def svg_bytes(doc):
    from PyQt6.QtCore import QBuffer, QByteArray, QIODevice
    from PyQt6.QtSvg import QSvgGenerator
    r, b = _frame(doc)
    data = QByteArray()
    buf = QBuffer(data); buf.open(QIODevice.OpenModeFlag.WriteOnly)
    gen = QSvgGenerator()
    gen.setOutputDevice(buf)
    gen.setResolution(72)                                    # 1 unit = 1 pt
    gen.setSize(QSize(math.ceil(b.width()), math.ceil(b.height())))
    gen.setViewBox(QRectF(0, 0, math.ceil(b.width()), math.ceil(b.height())))
    gen.setTitle('Structure')
    p = QPainter(gen)
    p.translate(-b.left(), -b.top())
    r.paint(p)
    p.end()
    buf.close()
    return bytes(data)


def save_svg(doc, path):
    with open(path, 'wb') as f:
        f.write(svg_bytes(doc))


def save_pdf(doc, path):
    from PyQt6.QtCore import QMarginsF
    from PyQt6.QtGui import QPageLayout, QPageSize, QPdfWriter
    r, b = _frame(doc)
    w = QPdfWriter(path)
    w.setResolution(72)
    w.setPageSize(QPageSize(b.size(), QPageSize.Unit.Point, 'Structure'))
    w.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Unit.Point)
    p = QPainter(w)
    p.translate(-b.left(), -b.top())
    r.paint(p)
    p.end()
