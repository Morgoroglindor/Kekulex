"""
icons.py - the tool and toolbar icons, drawn as line glyphs in the theme's text colour (no image files to ship).

    icons.icon('benzene')   QIcon at 1x and 2x; call again after a theme change

Each is drawn on a 24 x 24 grid.
"""
import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QIcon, QPainter, QPainterPath, QPen, QPixmap, QPolygonF

import theme


def _poly(n, cx=12, cy=12.5, r=8.0, start=-90):
    return [QPointF(cx + r * math.cos(math.radians(start + 360 * k / n)), cy + r * math.sin(math.radians(start + 360 * k / n)))
            for k in range(n)]


def _ring(n):
    def draw(p):
        p.drawPolygon(QPolygonF(_poly(n, r=8.5 if n > 4 else 9)))
    return draw


def _benzene(p):
    pts = _poly(6, r=8.5)
    p.drawPolygon(QPolygonF(pts))
    inner = _poly(6, r=6.0)
    for k in (0, 2, 4):
        a, b = inner[k], inner[(k + 1) % 6]
        p.drawLine(a + (b - a) * 0.12, b - (b - a) * 0.12)


def _bond(n):
    def draw(p):
        a, b = QPointF(4, 19), QPointF(20, 5)
        d = b - a
        L = math.hypot(d.x(), d.y())
        nx, ny = -d.y() / L, d.x() / L
        offs = {1: [0], 2: [-1.8, 1.8], 3: [-3, 0, 3]}[n]
        for o in offs:
            p.drawLine(a + QPointF(nx * o, ny * o), b + QPointF(nx * o, ny * o))
    return draw


def _wedge(p):
    p.setBrush(p.pen().color())
    p.drawPolygon(QPolygonF([QPointF(4, 19), QPointF(18.2, 3.4), QPointF(21.6, 6.6)]))


def _hash(p):
    a, b = QPointF(4, 19), QPointF(20, 5)
    for k in range(1, 7):
        t = k / 6
        c = a + (b - a) * t
        w = 3.2 * t
        p.drawLine(c + QPointF(-w * 0.66, -w * 0.75), c + QPointF(w * 0.66, w * 0.75))


def _chain(p):
    path = QPainterPath(QPointF(2.5, 16))
    for k, (x, y) in enumerate(((7.5, 8), (12, 16), (16.5, 8), (21.5, 16))):
        path.lineTo(x, y)
    p.drawPath(path)


def _select(p):                                            # an arrow pointer and a dashed box
    pen = QPen(p.pen()); pen.setStyle(Qt.PenStyle.DashLine); pen.setWidthF(1.0)
    old = p.pen(); p.setPen(pen)
    p.drawRect(QRectF(3, 3, 13, 11))
    p.setPen(old)
    p.setBrush(old.color())
    p.drawPolygon(QPolygonF([QPointF(11, 9), QPointF(11, 21.5), QPointF(14, 18.5), QPointF(16.5, 22.5),
                             QPointF(18.5, 21.5), QPointF(16, 17.5), QPointF(20, 17)]))


def _erase(p):
    p.save()
    p.translate(12, 12); p.rotate(-40)
    p.drawRoundedRect(QRectF(-8, -4.5, 16, 9), 1.5, 1.5)
    p.drawLine(QPointF(-2, -4.5), QPointF(-2, 4.5))
    p.restore()
    p.drawLine(QPointF(9, 21), QPointF(21, 21))


def _charge(sign):
    def draw(p):
        p.drawEllipse(QPointF(12, 12), 7.5, 7.5)
        p.drawLine(QPointF(8, 12), QPointF(16, 12))
        if sign > 0:
            p.drawLine(QPointF(12, 8), QPointF(12, 16))
    return draw


def _arrow(p):
    p.drawLine(QPointF(3, 12), QPointF(18, 12))
    p.setBrush(p.pen().color())
    p.drawPolygon(QPolygonF([QPointF(21.5, 12), QPointF(16, 8.5), QPointF(17.5, 12), QPointF(16, 15.5)]))


def _equilibrium(p):
    p.drawLine(QPointF(3, 9.5), QPointF(21, 9.5)); p.drawLine(QPointF(21, 9.5), QPointF(16.5, 6))
    p.drawLine(QPointF(3, 14.5), QPointF(21, 14.5)); p.drawLine(QPointF(3, 14.5), QPointF(7.5, 18))


def _text(p):
    f = QFont('Arial'); f.setPixelSize(17); f.setBold(False)
    p.setFont(f)
    p.drawText(QRectF(0, 0, 24, 24), Qt.AlignmentFlag.AlignCenter, 'T')


def _new(p):
    p.drawPolygon(QPolygonF([QPointF(5, 3), QPointF(14, 3), QPointF(19, 8), QPointF(19, 21), QPointF(5, 21)]))
    p.drawPolyline(QPolygonF([QPointF(14, 3), QPointF(14, 8), QPointF(19, 8)]))


def _open(p):
    p.drawPolygon(QPolygonF([QPointF(3, 6), QPointF(3, 19.5), QPointF(21, 19.5), QPointF(21, 8.5), QPointF(11.5, 8.5),
                             QPointF(9.5, 6)]))
    p.drawLine(QPointF(3, 11.5), QPointF(21, 11.5))


def _save(p):
    p.drawRoundedRect(QRectF(4, 4, 16, 16), 1.5, 1.5)
    p.drawRect(QRectF(8, 4, 8, 5))
    p.drawRect(QRectF(7, 13, 10, 7))


def _undo(p):
    path = QPainterPath(QPointF(7, 10))
    path.cubicTo(11, 5, 20, 7, 19, 15)
    p.drawPath(path)
    p.drawPolyline(QPolygonF([QPointF(6.5, 4.5), QPointF(6.5, 10.5), QPointF(12.5, 10.5)]))


def _redo(p):
    p.save(); p.translate(24, 0); p.scale(-1, 1); _undo(p); p.restore()


def _copy(p):
    p.drawRoundedRect(QRectF(4, 4, 11, 13), 1.2, 1.2)
    p.drawRoundedRect(QRectF(9, 8, 11, 13), 1.2, 1.2)


def _clean(p):                                             # a hexagon with sparkles
    p.drawPolygon(QPolygonF(_poly(6, cx=10, cy=13, r=7)))
    for cx, cy, r in ((19, 5, 3), (21, 12, 1.8)):
        p.drawLine(QPointF(cx - r, cy), QPointF(cx + r, cy)); p.drawLine(QPointF(cx, cy - r), QPointF(cx, cy + r))


def _fit(p):
    for (x, y, dx, dy) in ((3, 3, 1, 1), (21, 3, -1, 1), (3, 21, 1, -1), (21, 21, -1, -1)):
        p.drawPolyline(QPolygonF([QPointF(x + 6 * dx, y), QPointF(x, y), QPointF(x, y + 6 * dy)]))
    p.drawPolygon(QPolygonF(_poly(6, r=5, cy=12)))


def _image(p):
    p.drawRoundedRect(QRectF(3, 4.5, 18, 15), 1.5, 1.5)
    p.drawPolyline(QPolygonF([QPointF(4, 18), QPointF(9, 12), QPointF(13, 16), QPointF(16, 13), QPointF(20, 17.5)]))
    p.drawEllipse(QPointF(16, 9), 1.6, 1.6)


def _bold(p):
    p.save(); pen = QPen(p.pen()); pen.setWidthF(4.0); pen.setCapStyle(Qt.PenCapStyle.FlatCap); p.setPen(pen)
    p.drawLine(QPointF(5, 18), QPointF(19, 6)); p.restore()


def _dashed(p):
    pen = QPen(p.pen()); pen.setDashPattern([2.0, 1.6]); p.save(); p.setPen(pen)
    p.drawLine(QPointF(4, 19), QPointF(20, 5)); p.restore()


def _wavy(p):
    path = QPainterPath(QPointF(4, 19))
    for k in range(6):                                      # half waves along the diagonal
        x0, y0 = 4 + 16 * k / 6, 19 - 14 * k / 6
        x1, y1 = 4 + 16 * (k + 1) / 6, 19 - 14 * (k + 1) / 6
        s = 2.6 if k % 2 == 0 else -2.6
        path.quadTo((x0 + x1) / 2 + s * 0.66, (y0 + y1) / 2 + s * 0.75, x1, y1)
    p.drawPath(path)


def _hashed(p):
    a, b = QPointF(4, 19), QPointF(20, 5)
    for k in range(1, 7):
        c = a + (b - a) * (k / 6.5)
        p.drawLine(c + QPointF(-1.8, -2.1), c + QPointF(1.8, 2.1))


def _cyclopentadiene(p):
    pts = _poly(5, r=8.5)
    p.drawPolygon(QPolygonF(pts))
    inner = _poly(5, r=6.0)
    for k in (0, 2):
        a, b = inner[k], inner[(k + 1) % 5]
        p.drawLine(a + (b - a) * 0.15, b - (b - a) * 0.15)


def _head(p, x, y, dx):
    p.save(); p.setBrush(p.pen().color())
    p.drawPolygon(QPolygonF([QPointF(x, y), QPointF(x - 5.5 * dx, y - 3.5), QPointF(x - 4 * dx, y), QPointF(x - 5.5 * dx, y + 3.5)]))
    p.restore()


def _retro(p):
    p.drawLine(QPointF(3, 10), QPointF(17, 10)); p.drawLine(QPointF(3, 14), QPointF(17, 14))
    p.drawPolyline(QPolygonF([QPointF(15, 6), QPointF(21, 12), QPointF(15, 18)]))


def _resonance(p):
    p.drawLine(QPointF(6, 12), QPointF(18, 12))
    _head(p, 21.5, 12, 1); _head(p, 2.5, 12, -1)


def _arrow_dashed(p):
    pen = QPen(p.pen()); pen.setDashPattern([2.2, 1.6]); p.save(); p.setPen(pen)
    p.drawLine(QPointF(3, 12), QPointF(17.5, 12)); p.restore()
    _head(p, 21.5, 12, 1)


def _curved(half):
    def draw(p):
        path = QPainterPath(QPointF(3, 18))
        path.quadTo(12, 0, 19, 14)
        p.drawPath(path)
        p.save(); p.setBrush(p.pen().color())
        pts = [QPointF(20.5, 19), QPointF(15.5, 13.5), QPointF(19, 14.5)] if half else \
            [QPointF(20.5, 19), QPointF(15.5, 13.5), QPointF(19, 14.5), QPointF(22.5, 12.5)]
        p.drawPolygon(QPolygonF(pts)); p.restore()
    return draw


def _shape(kind):
    def draw(p):
        r = QRectF(4, 5, 16, 14)
        if kind == 'rect':
            p.drawRect(r)
        elif kind == 'roundrect':
            p.drawRoundedRect(r, 4, 4)
        elif kind == 'oval':
            p.drawEllipse(r)
        elif kind == 'line':
            p.drawLine(QPointF(4, 19), QPointF(20, 5))
        elif kind == 'brackets':
            p.drawPolyline(QPolygonF([QPointF(8, 4), QPointF(5, 4), QPointF(5, 20), QPointF(8, 20)]))
            p.drawPolyline(QPolygonF([QPointF(16, 4), QPointF(19, 4), QPointF(19, 20), QPointF(16, 20)]))
        elif kind == 'parens':
            for x, s in ((8, 1), (16, -1)):
                path = QPainterPath(QPointF(x, 4)); path.quadTo(x - 5 * s, 12, x, 20); p.drawPath(path)
        else:
            for x, s in ((8, 1), (16, -1)):
                path = QPainterPath(QPointF(x, 4)); path.quadTo(x - 3 * s, 4, x - 3 * s, 8); path.lineTo(x - 3 * s, 10)
                path.quadTo(x - 3 * s, 12, x - 5 * s, 12); path.quadTo(x - 3 * s, 12, x - 3 * s, 14)
                path.lineTo(x - 3 * s, 16); path.quadTo(x - 3 * s, 20, x, 20); p.drawPath(path)
    return draw


def _align(how):
    def draw(p):
        for k, (w, y) in enumerate(((14, 6), (9, 11), (12, 16))):
            x = 4 if how == 'left' else 20 - w if how == 'right' else 12 - w / 2
            p.drawLine(QPointF(x, y), QPointF(x + w, y))
    return draw


def _lasso(p):                                            # a loose loop with its tail
    path = QPainterPath(QPointF(8, 16))
    path.cubicTo(1, 13, 3, 4, 12, 4)
    path.cubicTo(21, 4, 22, 12, 15, 14)
    path.cubicTo(11, 15, 8, 15, 8, 16)
    path.cubicTo(8, 19, 10, 20, 9, 22)
    p.drawPath(path)


def _templates(p):                                         # a stack of hexagons
    p.drawPolygon(QPolygonF(_poly(6, cx=9, cy=9, r=5.5)))
    p.drawPolygon(QPolygonF(_poly(6, cx=15, cy=15, r=5.5)))


def _bio(kind):
    """A small copy of the bio object itself (in its own colours), fitted into the icon."""
    def draw(p):
        import biodraw
        import model
        import render
        doc = model.Doc()
        i = doc.add_shape(*biodraw.default_shape(kind, 0, 0), kind)
        if kind == 'plasmid':
            doc.shapes[i].params = {'bp': 100, 'regions': [{'start': 5, 'end': 40}, {'start': 55, 'end': 85}]}
        if kind in model.BIO_LINEAR:                        # a short piece, so its subunits stay visible
            sh = doc.shapes[i]
            sh.x1, sh.x2 = -1.6, 1.6
        r = render.Renderer(doc)
        b = biodraw.outline(r, doc.shapes[i]).boundingRect()
        k = 20 / max(b.width(), b.height(), 1e-6)
        p.save()
        p.translate(12, 12); p.scale(k, k); p.translate(-b.center().x(), -b.center().y())
        biodraw.paint_bio(r, p, i)
        p.restore()
    return draw


def _element(name):
    """A drawing element or orbital tool: the shape itself, drawn by the renderer, fitted into the icon."""
    def draw(p):
        import canvas
        import model
        import render
        kind, style, params = canvas.SHAPE_TOOLS[name]
        doc = model.Doc()
        doc.style = {'line_width': 1.0, 'bold_width': 3.0}
        if kind == 'circle':
            box = (0, 0, 0.62, 0)
        elif kind == 'line':
            box = (-0.7, 0.55, 0.7, -0.55)
        elif kind == 'arc':
            box = (-0.62, 0.25, 0.62, 0.25)
        elif kind == 'orbital':
            box = (0, 0, 0, -0.62) if params['type'] not in ('s', 'sigma') else (0, 0, 0.55, 0)
        else:
            box = (-0.68, -0.48, 0.68, 0.48)
        i = doc.add_shape(*box, kind)
        doc.shapes[i].style, doc.shapes[i].params = style, dict(params)
        r = render.Renderer(doc)
        b = r.shape_paths(i).boundingRect()
        k = 19 / max(b.width(), b.height(), 1e-6)
        p.save()
        p.translate(12, 12); p.scale(k, k); p.translate(-b.center().x(), -b.center().y())
        r.paint(p)
        p.restore()
        if kind == 'arc':
            f = QFont(); f.setPixelSize(7); f.setBold(True); p.setFont(f)
            p.drawText(QRectF(0, 15, 24, 9), Qt.AlignmentFlag.AlignCenter, str(params['angle']))
    return draw


def _massfrag(p):                                         # a bond cut by a dashed line
    p.drawLine(QPointF(3, 15), QPointF(11, 10)); p.drawLine(QPointF(11, 10), QPointF(20, 15))
    p.save(); pen = p.pen(); pen.setStyle(Qt.PenStyle.DashLine); pen.setWidthF(1.0); p.setPen(pen)
    p.drawLine(QPointF(16, 4), QPointF(14, 21)); p.restore()


def _retrosyn(p):                                         # the open (retrosynthetic) arrow
    p.drawLine(QPointF(3, 9), QPointF(18, 9)); p.drawLine(QPointF(3, 15), QPointF(18, 15))
    p.drawPolyline(QPolygonF([QPointF(15, 5), QPointF(21, 12), QPointF(15, 19)]))


def _synthesis(p):
    p.drawLine(QPointF(3, 12), QPointF(20, 12))
    p.drawPolyline(QPolygonF([QPointF(15, 8), QPointF(20, 12), QPointF(15, 16)]))


def _tlc(p):                                              # a plate with spots
    p.drawRect(QRectF(5, 2, 14, 20))
    p.save(); pen = p.pen(); pen.setWidthF(0.8); p.setPen(pen)
    p.drawLine(QPointF(5, 18.5), QPointF(19, 18.5))
    p.setBrush(p.pen().color())
    for x, y in ((9, 7), (15, 11), (9, 13)):
        p.drawEllipse(QPointF(x, y), 1.6, 1.2)
    p.restore()


def _send3d(p):                                           # a ball and stick model
    p.drawLine(QPointF(7, 16), QPointF(12, 8)); p.drawLine(QPointF(12, 8), QPointF(18, 15))
    p.save()
    for x, y, r, c in ((7, 16, 3.4, '#7FA7E0'), (18, 15, 3.0, '#E08A8A'), (12, 8, 4.2, '#B8B8B8')):
        p.setBrush(QColor(c)); p.drawEllipse(QPointF(x, y), r, r)
    p.restore()


def _sequence(p):                                         # three monomers on a chain
    p.drawLine(QPointF(3, 12), QPointF(21, 12))
    p.save()
    for x, c in ((6, '#E9B872'), (12, '#7FA7E0'), (18, '#8CCB9B')):
        p.setBrush(QColor(c)); p.drawEllipse(QPointF(x, 12), 3.6, 3.6)
    p.restore()


DRAW = {'select': _select, 'bond1': _bond(1), 'bond2': _bond(2), 'bond3': _bond(3), 'wedge': _wedge, 'hash': _hash,
        'chain': _chain, 'ring3': _ring(3), 'ring4': _ring(4), 'ring5': _ring(5), 'ring6': _ring(6), 'ring7': _ring(7),
        'ring8': _ring(8), 'benzene': _benzene, 'erase': _erase, 'charge+': _charge(1), 'charge-': _charge(-1),
        'arrow': _arrow, 'equilibrium': _equilibrium, 'text': _text, 'new': _new, 'open': _open, 'save': _save,
        'undo': _undo, 'redo': _redo, 'copy': _copy, 'clean': _clean, 'fit': _fit, 'image': _image,
        'bold': _bold, 'dashed': _dashed, 'wavy': _wavy, 'hashed': _hashed, 'cyclopentadiene': _cyclopentadiene,
        'retro': _retro, 'resonance': _resonance, 'arrow_dashed': _arrow_dashed, 'curved': _curved(False),
        'curved_half': _curved(True), 'rect': _shape('rect'), 'roundrect': _shape('roundrect'), 'oval': _shape('oval'),
        'line': _shape('line'), 'brackets': _shape('brackets'), 'parens': _shape('parens'), 'braces': _shape('braces'),
        'align_left': _align('left'), 'align_center': _align('center'), 'align_right': _align('right'),
        'templates': _templates, 'lasso': _lasso, 'sequence': _sequence, 'massfrag': _massfrag,
        'retrosyn': _retrosyn, 'synthesis': _synthesis, 'tlc': _tlc, 'send3d': _send3d}
DRAW.update({k: _bio(k) for k in ('membrane', 'membrane_arc', 'micelle', 'dna', 'helix', 'receptor', 'channel', 'enzyme',
                                  'gprotein', 'antibody', 'ribosome', 'trna', 'mrna', 'mito', 'nucleus', 'cell',
                                  'bacterium', 'virus', 'cloud', 'plasmid')})


def _pixmap(name, scale, color):
    pm = QPixmap(round(24 * scale), round(24 * scale))
    pm.setDevicePixelRatio(scale)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(QPen(color, 1.5, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
    p.setBrush(Qt.BrushStyle.NoBrush)
    (DRAW.get(name) or _element(name))(p)                 # drawing elements / orbitals: drawn by the renderer
    p.end()
    return pm


def icon(name):
    ic = QIcon()
    col = QColor(theme.T['text'])
    off = QColor(theme.T['disabled'])
    for s in (1.0, 2.0):
        ic.addPixmap(_pixmap(name, s, col), QIcon.Mode.Normal)
        ic.addPixmap(_pixmap(name, s, off), QIcon.Mode.Disabled)
    return ic
