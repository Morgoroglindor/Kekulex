"""
tlc.py - thin-layer chromatography plates (model.TLC), as ChemDraw's TLC tool draws them: a plate with an origin line,
a solvent front and lanes; each lane holds spots placed by their Rf (0 at the origin, 1 at the front).

    i = new_plate(doc, x, y, w, h)        # lanes from the width (about one per bond length), one spot each at Rf 0.5
    layout(renderer, i)                   # geometry in pt: plate rect, origin / front y, lane x, spot rects
    spot_at(renderer, i, x_pt, y_pt)      # (lane, spot index) under a point, or None
    rf_at(renderer, i, y_pt)              # the Rf of a height on the plate (dragging a spot)
"""
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QPen, QRadialGradient

import model

SPOT_STYLES = (('filled', 'Filled'), ('shaded', 'Shaded'), ('hollow', 'Hollow'))


def new_spot(rf=0.5, **kw):
    return {'rf': round(float(rf), 3), 'w': 1.0, 'h': 1.0, 'style': 'filled', 'color': '', 'show_rf': False, **kw}


def new_plate(doc, x, y, w=3.0, h=4.0):
    x, w = (x, w) if w >= 0 else (x + w, -w)
    y, h = (y, h) if h >= 0 else (y + h, -h)
    w, h = max(w, 0.8), max(h, 1.2)
    n = max(1, int(round(w / 1.0)))
    i = doc._id()
    doc.tlcs[i] = model.TLC(x, y, w, h, lanes=[[new_spot(0.5)] for _ in range(n)])
    return i


def layout(r, i):
    t = r.doc.tlcs[i]
    BL = r.BL
    rect = QRectF(t.x * BL, t.y * BL, max(t.w * BL, 4), max(t.h * BL, 4))
    oy = rect.bottom() - t.origin * rect.height()
    fy = rect.top() + t.front * rect.height()
    n = max(1, len(t.lanes))
    lane_w = rect.width() / n
    lanes = [rect.left() + lane_w * (k + 0.5) for k in range(len(t.lanes))]
    base = min(lane_w * 0.17, BL * 0.22)
    spots = {}
    for k, lane in enumerate(t.lanes):
        for j, s in enumerate(lane):
            cy = oy - s['rf'] * (oy - fy)
            rx, ry = base * s.get('w', 1.0), base * 0.72 * s.get('h', 1.0)
            spots[(k, j)] = QRectF(lanes[k] - rx, cy - ry, 2 * rx, 2 * ry)
    full = QRectF(rect)
    for k, j in spots:
        if t.lanes[k][j].get('show_rf'):
            sr = spots[(k, j)]
            full = full.united(QRectF(sr.right(), sr.top(), r.m.width('0.00', r.s['font_size'] * 0.75) + 4, sr.height()))
    return {'rect': full, 'plate': rect, 'origin': oy, 'front': fy, 'lanes': lanes, 'lane_w': lane_w, 'spots': spots}


def spot_at(r, i, x, y):
    L = layout(r, i)
    pt = QPointF(x, y)
    for key, sr in L['spots'].items():
        if sr.adjusted(-1.5, -1.5, 1.5, 1.5).contains(pt):
            return key
    return None


def lane_at(r, i, x):
    L = layout(r, i)
    k = int((x - L['plate'].left()) / L['lane_w']) if L['lane_w'] else 0
    return max(0, min(len(L['lanes']) - 1, k)) if L['lanes'] else 0


def rf_at(r, i, y):
    L = layout(r, i)
    span = (L['origin'] - L['front']) or 1.0
    return max(0.0, min(1.0, (L['origin'] - y) / span))


def paint(r, p, i):
    import render
    t = r.doc.tlcs[i]
    L = layout(r, i)
    plate = L['plate']
    ink = QColor(t.color) if t.color else QColor('#000000')
    lw = r.s['line_width']
    p.save()
    if t.show_border:
        p.setPen(QPen(ink, lw)); p.setBrush(QBrush(QColor('#FFFFFF')))
        p.drawRect(plate)
    tick = min(L['lane_w'] * 0.12, r.BL * 0.15)
    for y, show, dashed in ((L['origin'], t.show_origin, False), (L['front'], t.show_front, True)):
        if not show:
            continue
        pen = QPen(ink, lw * 0.8)
        if dashed:
            pen.setStyle(Qt.PenStyle.DashLine)
        p.setPen(pen)
        p.drawLine(QPointF(plate.left(), y), QPointF(plate.right(), y))
        if t.show_ticks:                                    # side ticks, outside the plate
            p.setPen(QPen(ink, lw * 0.8))
            p.drawLine(QPointF(plate.left() - tick, y), QPointF(plate.left(), y))
            p.drawLine(QPointF(plate.right(), y), QPointF(plate.right() + tick, y))
    if t.show_origin and t.show_ticks:                      # where each lane was spotted
        p.setPen(QPen(ink, lw * 0.8))
        for x in L['lanes']:
            p.drawLine(QPointF(x, L['origin'] - tick * 0.6), QPointF(x, L['origin'] + tick * 0.6))
    for (k, j), sr in L['spots'].items():
        s = t.lanes[k][j]
        col = QColor(s.get('color') or '#000000')
        style = s.get('style', 'filled')
        p.setPen(QPen(col, lw))
        if style == 'filled':
            p.setBrush(QBrush(col))
        elif style == 'shaded':
            g = QRadialGradient(sr.center(), max(sr.width(), sr.height()) * 0.6)
            g.setColorAt(0.0, col)
            fade = QColor(col); fade.setAlpha(40)
            g.setColorAt(1.0, fade)
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(g))
        else:
            p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawEllipse(sr)
        if s.get('show_rf'):
            size = r.s['font_size'] * 0.75
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(ink))
            p.drawPath(render.text_path(f"{s['rf']:.2f}", r.s['font'], size, sr.right() + 2,
                                        sr.center().y() + r.m.cap(size) / 2))
    p.restore()
