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


DRAW = {'select': _select, 'bond1': _bond(1), 'bond2': _bond(2), 'bond3': _bond(3), 'wedge': _wedge, 'hash': _hash,
        'chain': _chain, 'ring3': _ring(3), 'ring4': _ring(4), 'ring5': _ring(5), 'ring6': _ring(6), 'ring7': _ring(7),
        'ring8': _ring(8), 'benzene': _benzene, 'erase': _erase, 'charge+': _charge(1), 'charge-': _charge(-1),
        'arrow': _arrow, 'equilibrium': _equilibrium, 'text': _text, 'new': _new, 'open': _open, 'save': _save,
        'undo': _undo, 'redo': _redo, 'copy': _copy, 'clean': _clean, 'fit': _fit, 'image': _image}


def _pixmap(name, scale, color):
    pm = QPixmap(round(24 * scale), round(24 * scale))
    pm.setDevicePixelRatio(scale)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(QPen(color, 1.5, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
    p.setBrush(Qt.BrushStyle.NoBrush)
    DRAW[name](p)
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
