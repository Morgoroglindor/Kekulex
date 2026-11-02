"""
elements.py - ChemDraw-like drawing elements and orbitals for render.Renderer.

Drawing elements: circle (dragged out from its centre), oval, rounded rectangle and rectangle, each plain, dashed,
shaded (lit like a sphere), filled or with a shadow; lines plain, dashed, bold or wavy; arcs of 90, 120, 180 or 270
degrees (plain or dashed) whose ends are the drag's ends - left to right bulges up (convex), right to left down.

Orbitals (Shape kind 'orbital', params['type']): s, sigma, single lobe, p, hybrid (sp3: a big and a small lobe),
dxy (four lobes) and dz2 (two lobes and a ring). x1,y1 is the node (usually an atom), x2,y2 the tip of the main lobe.
Shape.style gives the phase look: '' both lobes outlined, 'shaded' the + lobe shaded, 'filled' the + lobe solid.
"""
import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QPainterPath, QPen, QRadialGradient, QTransform

ORBITALS = {'s': 's orbital', 'sigma': 'σ orbital', 'lobe': 'Single lobe', 'p': 'p orbital',
            'hybrid': 'Hybrid (sp³) orbital', 'dxy': 'dxy orbital', 'dz2': 'dz² orbital'}
ARC_ANGLES = (270, 180, 120, 90)
CLOSED = ('circle', 'oval', 'rect', 'roundrect')


def _color(r, sh):
    return QColor(sh.color) if sh.color else QColor(r.s.get('color', '#000000') or '#000000')


def _box(sh, BL):
    if sh.kind == 'circle':
        R = math.hypot(sh.x2 - sh.x1, sh.y2 - sh.y1) * BL
        return QRectF(sh.x1 * BL - R, sh.y1 * BL - R, 2 * R, 2 * R)
    x1, x2 = sorted((sh.x1 * BL, sh.x2 * BL))
    y1, y2 = sorted((sh.y1 * BL, sh.y2 * BL))
    return QRectF(x1, y1, max(x2 - x1, 0.5), max(y2 - y1, 0.5))


def _closed_path(sh, b):
    p = QPainterPath()
    if sh.kind in ('circle', 'oval'):
        p.addEllipse(b)
    elif sh.kind == 'roundrect':
        rr = min(b.width(), b.height()) * 0.15
        p.addRoundedRect(b, rr, rr)
    else:
        p.addRect(b)
    return p


def arc_points(sh, BL, n=48):
    """Points (pt) along an arc from x1,y1 to x2,y2 sweeping params['angle'] degrees, bulging to the left of the
    direction of the drag (up for a left-to-right drag)."""
    x1, y1, x2, y2 = sh.x1 * BL, sh.y1 * BL, sh.x2 * BL, sh.y2 * BL
    theta = math.radians(float(sh.params.get('angle', 180)))
    c = math.hypot(x2 - x1, y2 - y1) or 1e-6
    R = c / (2 * math.sin(theta / 2))
    ux, uy = (x2 - x1) / c, (y2 - y1) / c
    nx, ny = uy, -ux                                        # left of the drag direction (y down)
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    d = R * math.cos(theta / 2)
    cx, cy = mx - nx * d, my - ny * d                       # centre on the far side for arcs under 180 degrees
    a1 = math.atan2(y1 - cy, x1 - cx)
    bulge = (mx + nx * (R - d), my + ny * (R - d))
    best = None
    for sign in (1, -1):
        am = a1 + sign * theta / 2
        dist = math.hypot(cx + R * math.cos(am) - bulge[0], cy + R * math.sin(am) - bulge[1])
        if best is None or dist < best[0]:
            best = (dist, sign)
    sign = best[1]
    return [(cx + R * math.cos(a1 + sign * theta * k / n), cy + R * math.sin(a1 + sign * theta * k / n))
            for k in range(n + 1)]


def _wavy(x1, y1, x2, y2, amp, period):
    L = math.hypot(x2 - x1, y2 - y1) or 1e-6
    ux, uy = (x2 - x1) / L, (y2 - y1) / L
    p = QPainterPath(QPointF(x1, y1))
    steps = max(8, int(L / 1.0))
    for k in range(1, steps + 1):
        t = L * k / steps
        off = amp * math.sin(2 * math.pi * t / period)
        p.lineTo(x1 + ux * t - uy * off, y1 + uy * t + ux * off)
    return p


# ------------------------------------------------------------------------------------------------ orbitals
def _lobe(L, W):
    """A lobe along +x from the node (0, 0) to its round end at (L, 0)."""
    p = QPainterPath(QPointF(0, 0))
    p.cubicTo(0.08 * L, -1.05 * W, 0.92 * L, -1.25 * W, L, 0)
    p.cubicTo(0.92 * L, 1.25 * W, 0.08 * L, 1.05 * W, 0, 0)
    return p


def orbital_parts(sh, BL):
    """[(QPainterPath in pt, phase +1 / -1)] of an orbital, back to front."""
    x1, y1, x2, y2 = sh.x1 * BL, sh.y1 * BL, sh.x2 * BL, sh.y2 * BL
    L = max(math.hypot(x2 - x1, y2 - y1), 2.0)
    t = QTransform()
    t.translate(x1, y1)
    t.rotate(math.degrees(math.atan2(y2 - y1, x2 - x1)))
    kind = sh.params.get('type', 'p')
    W = 0.42 * L
    parts = []
    if kind == 's':
        c = QPainterPath(); c.addEllipse(QPointF(0, 0), L, L)
        parts.append((c, 1))
    elif kind == 'sigma':
        c = QPainterPath(); c.addEllipse(QPointF(0, 0), L, 0.55 * L)
        parts.append((c, 1))
    elif kind == 'lobe':
        parts.append((_lobe(L, W), 1))
    elif kind == 'p':
        parts += [(_lobe(L, W), 1), (QTransform().rotate(180).map(_lobe(L, W)), -1)]
    elif kind == 'hybrid':
        parts += [(QTransform().rotate(180).map(_lobe(0.38 * L, 0.42 * 0.38 * L * 1.15)), -1), (_lobe(L, W * 1.15), 1)]
    elif kind == 'dxy':
        for k, ang in enumerate((45, 135, 225, 315)):
            parts.append((QTransform().rotate(ang).map(_lobe(L, W * 0.85)), 1 if k % 2 == 0 else -1))
    elif kind == 'dz2':
        ring = QPainterPath(); ring.addEllipse(QPointF(0, 0), 0.24 * L, 0.7 * L)
        parts += [(ring, -1), (_lobe(L, W * 0.8), 1), (QTransform().rotate(180).map(_lobe(L, W * 0.8)), 1)]
    return [(t.map(path), phase) for path, phase in parts]


# ------------------------------------------------------------------------------------------------ public
def outline(r, sh):
    """Outline path (pt) for selection marks, bounds and hit tests."""
    BL = r.BL
    if sh.kind == 'arc':
        pts = arc_points(sh, BL)
        p = QPainterPath(QPointF(*pts[0]))
        for x, y in pts[1:]:
            p.lineTo(x, y)
        return p
    if sh.kind == 'orbital':
        p = QPainterPath()
        for path, _phase in orbital_parts(sh, BL):
            p.addPath(path)
        return p
    if sh.kind == 'line':
        p = QPainterPath(QPointF(sh.x1 * BL, sh.y1 * BL)); p.lineTo(sh.x2 * BL, sh.y2 * BL)
        return p
    b = _box(sh, BL)
    p = _closed_path(sh, b)
    if sh.style == 'shadow':
        d = _shadow_offset(r)
        p.addPath(_closed_path(sh, b.translated(d, d)))
    return p


def _shadow_offset(r):
    return max(1.5, r.s['line_width'] * 3.5)


def _shade_brush(rect, color, sphere=True):
    """A gradient lit from the top left (spheres) or from the middle (boxes): white -> the colour."""
    if sphere:
        g = QRadialGradient(QPointF(rect.left() + rect.width() * 0.32, rect.top() + rect.height() * 0.3),
                            max(rect.width(), rect.height()) * 0.85)
    else:
        g = QRadialGradient(rect.center(), max(rect.width(), rect.height()) * 0.62)
    light = QColor(color).lighter(100)
    g.setColorAt(0.0, QColor('#FFFFFF'))
    mid = QColor(color); mid.setAlpha(255)
    g.setColorAt(0.75 if sphere else 1.0, mid if color.lightness() > 40 else QColor('#3A3A3A'))
    g.setColorAt(1.0, light.darker(170) if color.lightness() > 40 else QColor('#000000'))
    return QBrush(g)


def paint(r, p, i):
    sh = r.doc.shapes[i]
    BL = r.BL
    color = _color(r, sh)
    lw = r.s['line_width']
    pen = QPen(color, lw, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin)
    if sh.style == 'dashed':
        pen.setStyle(Qt.PenStyle.DashLine)
    p.save()
    if sh.kind == 'orbital':
        _paint_orbital(r, p, sh, color, lw)
    elif sh.kind in ('line', 'arc'):
        if sh.style == 'bold':
            pen.setWidthF(r.s['bold_width'])
            pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
        if sh.kind == 'line' and sh.style == 'wavy':
            p.drawPath(_wavy(sh.x1 * BL, sh.y1 * BL, sh.x2 * BL, sh.y2 * BL, BL * 0.07, BL * 0.25))
        else:
            p.drawPath(outline(r, sh))
    else:
        b = _box(sh, BL)
        path = _closed_path(sh, b)
        sphere = sh.kind in ('circle', 'oval')
        if sh.style == 'shadow':
            d = _shadow_offset(r)
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor('#A6A6A6')))
            p.drawPath(_closed_path(sh, b.translated(d, d)))
            p.setBrush(QBrush(QColor('#FFFFFF')))
        elif sh.style == 'filled':
            p.setBrush(QBrush(color))
        elif sh.style == 'shaded':
            p.setBrush(_shade_brush(b, color, sphere))
        else:
            p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(pen)
        p.drawPath(path)
    p.restore()


def _paint_orbital(r, p, sh, color, lw):
    for path, phase in orbital_parts(sh, r.BL):
        b = path.boundingRect()
        if sh.style == 'filled':
            p.setBrush(QBrush(color if phase > 0 else QColor('#FFFFFF')))
        elif sh.style == 'shaded':
            p.setBrush(_shade_brush(b, color, True) if phase > 0 else QBrush(QColor('#FFFFFF')))
        else:
            p.setBrush(QBrush(QColor('#FFFFFF')))
        p.setPen(QPen(color, lw, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
        p.drawPath(path)


def contains(sh, x, y, tol, BL=14.4):
    """Hit test in model units: filled / shaded / shadowed shapes and orbitals count inside, the others on the
    outline (they may frame a structure)."""
    from model import seg_dist
    if sh.kind == 'orbital':
        from PyQt6.QtCore import QPointF as P
        for path, _ph in orbital_parts(sh, BL):
            if path.contains(P(x * BL, y * BL)):
                return True
        return math.hypot(x - sh.x1, y - sh.y1) <= tol
    if sh.kind == 'arc':
        pts = [(px / BL, py / BL) for px, py in arc_points(sh, BL, 24)]
        return min(seg_dist(x, y, *pts[k], *pts[k + 1]) for k in range(len(pts) - 1)) <= tol
    if sh.kind == 'line':
        return seg_dist(x, y, sh.x1, sh.y1, sh.x2, sh.y2) <= tol
    if sh.kind == 'circle':
        R = math.hypot(sh.x2 - sh.x1, sh.y2 - sh.y1)
        d = math.hypot(x - sh.x1, y - sh.y1)
        return d <= R + tol if sh.style in ('filled', 'shaded', 'shadow') else abs(d - R) <= tol
    return False
