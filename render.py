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

from PyQt6.QtCore import QPointF, QRectF, QSize, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QFontMetricsF, QImage, QPainter, QPainterPath, QPen, QPolygonF, QTransform

DEFAULT_STYLE = dict(bond_length=14.4, line_width=0.6, bold_width=2.0, hash_spacing=2.5, bond_spacing=0.18, margin=1.6,
                     font='Arial', font_size=10.0, color='#000000', color_atoms=False)
ACS_STYLE = dict(DEFAULT_STYLE)
ELEMENT_COLORS = {'N': '#2846C8', 'O': '#E01E1E', 'S': '#B89400', 'P': '#E07000', 'F': '#2A9A2A', 'Cl': '#1C9A1C',
                  'Br': '#A0291A', 'I': '#7A1FA0', 'B': '#C07070', 'Si': '#B08000', 'Se': '#B07000'}
REF = 100.0                                                 # text is shaped at 100 px and scaled down
_paths = {}


def style_of(doc):
    s = dict(DEFAULT_STYLE)
    s.update({k: v for k, v in (doc.style or {}).items() if k in DEFAULT_STYLE})
    return s


def _font(family):
    f = QFont(family)
    f.setPixelSize(int(REF))
    f.setHintingPreference(QFont.HintingPreference.PreferNoHinting)
    return f


def text_path(text, family, size, x, y):
    """Outline of text with its baseline starting at (x, y), size in pt."""
    key = (text, family)
    if key not in _paths:
        p = QPainterPath()
        p.addText(0, 0, _font(family), text)
        _paths[key] = p
    k = size / REF
    return QTransform(k, 0, 0, k, x, y).map(_paths[key])


class Metrics:
    def __init__(self, family):
        self.fm = QFontMetricsF(_font(family))

    def width(self, text, size):
        return self.fm.horizontalAdvance(text) * size / REF

    def cap(self, size):
        return self.fm.capHeight() * size / REF

    def ascent(self, size):
        return self.fm.ascent() * size / REF

    def line(self, size):
        return self.fm.lineSpacing() * size / REF


def exit_t(ox, oy, ux, uy, hw, hh):
    """How far along (ux, uy) a line starting at offset (ox, oy) from a box centre leaves the box (half sizes hw, hh)."""
    t = math.inf
    if abs(ux) > 1e-9:
        t = min(t, ((hw if ux > 0 else -hw) - ox) / ux)
    if abs(uy) > 1e-9:
        t = min(t, ((hh if uy > 0 else -hh) - oy) / uy)
    return max(t, 0.0) if t != math.inf else 0.0


class Renderer:
    def __init__(self, doc):
        self.doc = doc
        self.s = style_of(doc)
        self.BL = self.s['bond_length']
        self.m = Metrics(self.s['font'])
        self.labels = {a: self._label(a) for a in doc.atoms}

    def pt(self, x, y):
        return QPointF(x * self.BL, y * self.BL)

    def atom_color(self, a):
        if self.s['color_atoms']:
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
        small = fs * 0.72
        m = self.m
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
        ch = m.cap(fs)
        ws = m.width(sym, fs)
        base = cy + ch / 2
        x0 = cx - ws / 2
        pieces = [(sym, x0, base, fs)]
        right_end = x0 + ws
        dxs = [doc.atoms[n].x - at.x for n, _b in nb]
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
        if at.charge:
            t = self._charge_text(at.charge)
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

    # ------------------------------------------------------------------------------------------- bonds
    def _end(self, a, ox, oy, ux, uy):
        """Distance along (ux, uy) from atom a (offset ox, oy) to where a line may start (0 = unlabelled atom)."""
        lab = self.labels.get(a)
        if not lab or not lab['clip']:
            return 0.0
        return exit_t(ox, oy, ux, uy, *lab['clip'])

    def _line(self, path, a, b, ox=0.0, oy=0.0, trim_a=0.0, trim_b=0.0):
        p, q = self.doc.atoms[a], self.doc.atoms[b]
        x1, y1, x2, y2 = p.x * self.BL, p.y * self.BL, q.x * self.BL, q.y * self.BL
        L = math.hypot(x2 - x1, y2 - y1)
        if L < 1e-6:
            return
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        t1 = max(self._end(a, ox, oy, ux, uy), trim_a)
        t2 = max(self._end(b, ox, oy, -ux, -uy), trim_b)
        if t1 + t2 >= L - 0.2:
            return
        path.moveTo(x1 + ox + ux * t1, y1 + oy + uy * t1)
        path.lineTo(x2 + ox - ux * t2, y2 + oy - uy * t2)

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
            self._line(stroke, bd.a, bd.b)
        elif bd.order == 2:
            side = doc.double_side(i)
            if side == 0:
                self._line(stroke, bd.a, bd.b, nx * sp / 2, ny * sp / 2)
                self._line(stroke, bd.a, bd.b, -nx * sp / 2, -ny * sp / 2)
            else:
                self._line(stroke, bd.a, bd.b)
                trim = 0.75 * sp
                ta = 0.0 if bd.a in self.labels and self.labels[bd.a] and self.labels[bd.a]['clip'] else trim
                tb = 0.0 if bd.b in self.labels and self.labels[bd.b] and self.labels[bd.b]['clip'] else trim
                if len(doc.neighbors(bd.a)) == 1:
                    ta = 0.0
                if len(doc.neighbors(bd.b)) == 1:
                    tb = 0.0
                self._line(stroke, bd.a, bd.b, nx * sp * side, ny * sp * side, ta, tb)
        else:
            self._line(stroke, bd.a, bd.b)
            self._line(stroke, bd.a, bd.b, nx * sp, ny * sp)
            self._line(stroke, bd.a, bd.b, -nx * sp, -ny * sp)
        return stroke, fill

    # ------------------------------------------------------------------------------------------- arrows, text
    def arrow_paths(self, i):
        ar = self.doc.arrows[i]
        x1, y1, x2, y2 = ar.x1 * self.BL, ar.y1 * self.BL, ar.x2 * self.BL, ar.y2 * self.BL
        stroke, fill = QPainterPath(), QPainterPath()
        L = math.hypot(x2 - x1, y2 - y1)
        if L < 1e-6:
            return stroke, fill
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        nx, ny = -uy, ux
        hl, hw, notch = 7.0, 2.6, 5.2                       # head length, half width, where its back meets the line
        if ar.kind == 'equilibrium':
            off = 1.8
            for sgn, (ax, ay, bx, by) in ((1, (x1, y1, x2, y2)), (-1, (x2, y2, x1, y1))):
                dx, dy = (bx - ax) / L, (by - ay) / L
                ox, oy = -nx * off * sgn, -ny * off * sgn
                sx, sy, ex, ey = ax + ox, ay + oy, bx + ox, by + oy
                stroke.moveTo(sx, sy); stroke.lineTo(ex - dx * notch * 0.5, ey - dy * notch * 0.5)
                bxx, byy = -nx * sgn, -ny * sgn                  # the barb points away from the other arrow
                fill.addPolygon(QPolygonF([QPointF(ex, ey), QPointF(ex - dx * hl + bxx * hw, ey - dy * hl + byy * hw),
                                           QPointF(ex - dx * notch, ey - dy * notch)]))
                fill.closeSubpath()
            return stroke, fill
        stroke.moveTo(x1, y1); stroke.lineTo(x2 - ux * notch, y2 - uy * notch)
        fill.addPolygon(QPolygonF([QPointF(x2, y2), QPointF(x2 - ux * hl + nx * hw, y2 - uy * hl + ny * hw),
                                   QPointF(x2 - ux * notch, y2 - uy * notch),
                                   QPointF(x2 - ux * hl - nx * hw, y2 - uy * hl - ny * hw)]))
        fill.closeSubpath()
        return stroke, fill

    def text_lines(self, i):
        t = self.doc.texts[i]
        fs = self.s['font_size']
        x, y = t.x * self.BL, t.y * self.BL
        lh = self.m.line(fs)
        return [(line, x, y + self.m.ascent(fs) + k * lh) for k, line in enumerate(t.text.split('\n'))]

    def text_rect(self, i):
        """QRectF (pt) of text item i."""
        t = self.doc.texts[i]
        fs = self.s['font_size']
        lines = t.text.split('\n')
        w = max([self.m.width(line, fs) for line in lines] + [fs * 0.6])
        return QRectF(t.x * self.BL, t.y * self.BL, w, self.m.line(fs) * len(lines))

    # ------------------------------------------------------------------------------------------- painting
    def paint(self, p, errors=False):
        """Draw everything (errors: a red box around atoms over their valence, for the screen)."""
        s = self.s
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        color = QColor(s['color'])
        pen = QPen(color, s['line_width'], Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
        for i in self.doc.bonds:
            stroke, fill = self.bond_paths(i)
            p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(stroke)
            if not fill.isEmpty():
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color))
                p.drawPath(fill)
        for a, lab in self.labels.items():
            if not lab:
                continue
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(self.atom_color(a)))
            for text, x, y, size in lab['pieces']:
                p.drawPath(text_path(text, s['font'], size, x, y))
        for i in self.doc.arrows:
            stroke, fill = self.arrow_paths(i)
            p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush); p.drawPath(stroke)
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color)); p.drawPath(fill)
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color))
        for i in self.doc.texts:
            for line, x, y in self.text_lines(i):
                if line:
                    p.drawPath(text_path(line, s['font'], s['font_size'], x, y))
        if errors:
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.setPen(QPen(QColor('#E0201B'), 0.8))
            for a, at in self.doc.atoms.items():
                if self.doc.valence_error(a):
                    r = 0.28 * self.BL
                    p.drawRect(QRectF(at.x * self.BL - r, at.y * self.BL - r, 2 * r, 2 * r))

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
            add(stroke.boundingRect().united(fill.boundingRect()))
        for i in self.doc.texts:
            add(self.text_rect(i))
        return r


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
