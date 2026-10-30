"""
biodraw.py - how render.Renderer draws biopolymer sequences and BioDraw-like objects.

Sequences (model.Sequence): a row of monomers wrapped every per_line residues - peptides as circles (pills in the
three-letter style), DNA / RNA as rounded squares - coloured by residue (biopolymer.colour), joined by short bonds,
with the termini (H- / -OH, 5' / 3'), optional residue numbers, the head-to-tail link of a cyclic peptide and the
disulfide / lactam bridges as arcs.

Bio objects (model.Shape with a kind in model.BIO_KINDS): membranes, micelle, DNA double helix, helix protein,
receptor, ion channel, enzyme, G protein, antibody, ribosome, tRNA, mRNA, mitochondrion, nucleus, cell, bacterium,
virus, cloud and the plasmid map. The "linear" kinds (model.BIO_LINEAR) run along the line x1,y1 -> x2,y2 and grow
more subunits as they get longer (like ChemDraw's); the others fill their box. The artwork is this program's own:
simple shapes from paths, each kind with a default fill (Shape.fill overrides it, Shape.color the outline).
"""
import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QPainterPath, QPen, QTransform

import biopolymer
import model

INK = '#3A4049'
NAMES = {'membrane': 'Membrane (bilayer)', 'membrane_arc': 'Membrane arc', 'micelle': 'Micelle', 'dna': 'DNA double helix',
         'helix': 'Helix protein (transmembrane)', 'receptor': 'Receptor', 'channel': 'Ion channel', 'enzyme': 'Enzyme',
         'gprotein': 'G protein (α β γ)', 'antibody': 'Antibody (IgG)', 'ribosome': 'Ribosome', 'trna': 'tRNA',
         'mrna': 'mRNA', 'mito': 'Mitochondrion', 'nucleus': 'Nucleus', 'cell': 'Cell', 'bacterium': 'Bacterium',
         'virus': 'Virus', 'cloud': 'Cloud (pathway label)', 'plasmid': 'Plasmid map'}
FILL = {'membrane': '#F2A65A', 'membrane_arc': '#F2A65A', 'micelle': '#F2A65A', 'dna': '#4A7BC8', 'helix': '#7FB685',
        'receptor': '#6FA8DC', 'channel': '#E6B35A', 'enzyme': '#9BC59D', 'gprotein': '#B39DDB', 'antibody': '#5B8DD9',
        'ribosome': '#D9A066', 'trna': '#E8A0A0', 'mrna': '#C77DB5', 'mito': '#F4B183', 'nucleus': '#C9B6E4',
        'cell': '#FBE3D6', 'bacterium': '#9FD3C7', 'virus': '#E27D7D', 'cloud': '#FFFFFF', 'plasmid': '#5A6270'}
# default size in bond lengths: (length, thickness) for linear kinds, (width, height) for the others
SIZE = {'membrane': (6, 1.5), 'dna': (6, 1.3), 'helix': (5, 1.8), 'mrna': (5, 0.6), 'membrane_arc': (6, 2.5),
        'micelle': (3, 3), 'receptor': (1.4, 2.4), 'channel': (1.6, 2.2), 'enzyme': (2.2, 1.6), 'gprotein': (2.4, 1.4),
        'antibody': (2.2, 2.4), 'ribosome': (2.2, 2.0), 'trna': (1.8, 2.4), 'mito': (3.6, 1.8), 'nucleus': (2.4, 2.0),
        'cell': (5, 3.6), 'bacterium': (3.6, 1.4), 'virus': (2.2, 2.2), 'cloud': (3.2, 1.8), 'plasmid': (8, 8)}
PLASMID_COLOURS = ('#E07B53', '#4A7BC8', '#7FB685', '#B39DDB', '#E6B35A', '#5BB5C8')


def default_shape(kind, x, y):
    """(x1, y1, x2, y2) of an object of the default size placed with a click at x, y (bond lengths)."""
    a, b = SIZE.get(kind, (2, 2))
    if kind in model.BIO_LINEAR:
        return x - a / 2, y, x + a / 2, y
    return x - a / 2, y - b / 2, x + a / 2, y + b / 2


def thickness(sh):
    """Thickness of a linear object, bond lengths."""
    return float(sh.params.get('height', SIZE.get(sh.kind, (0, 1.2))[1]))


# ------------------------------------------------------------------------------------------------ geometry
def _linear_frame(sh, BL):
    """(QTransform from local (along, across) to pt, length pt, thickness pt) of a linear object."""
    x1, y1, x2, y2 = sh.x1 * BL, sh.y1 * BL, sh.x2 * BL, sh.y2 * BL
    t = QTransform()
    t.translate(x1, y1)
    t.rotate(math.degrees(math.atan2(y2 - y1, x2 - x1)))
    return t, max(math.hypot(x2 - x1, y2 - y1), 1.0), thickness(sh) * BL


def _box(sh, BL):
    x1, x2 = sorted((sh.x1 * BL, sh.x2 * BL))
    y1, y2 = sorted((sh.y1 * BL, sh.y2 * BL))
    return QRectF(x1, y1, max(x2 - x1, 1.0), max(y2 - y1, 1.0))


def outline(r, sh):
    """The object's outline path in pt (selection marks, hit tests, bounds)."""
    p = QPainterPath()
    if sh.kind in model.BIO_LINEAR:
        t, L, H = _linear_frame(sh, r.BL)
        p.addRect(QRectF(0, -H / 2, L, H))
        return t.map(p)
    b = _box(sh, r.BL)
    if sh.kind == 'plasmid':                                # labels stick out of the ring
        return _plasmid_bounds(r, sh, b)
    p.addRect(b)
    return p


def contains(sh, x, y, tol):
    """Whether model point x, y (bond lengths) is on the object (its inside counts)."""
    if sh.kind in model.BIO_LINEAR:
        return model.seg_dist(x, y, sh.x1, sh.y1, sh.x2, sh.y2) <= thickness(sh) / 2 + tol
    x1, x2 = sorted((sh.x1, sh.x2)); y1, y2 = sorted((sh.y1, sh.y2))
    return x1 - tol <= x <= x2 + tol and y1 - tol <= y <= y2 + tol


# ------------------------------------------------------------------------------------------------ painting
class _Ops:
    """Collects what to draw: filled / stroked paths and text, in the object's local coordinates."""

    def __init__(self, r, ink):
        self.r, self.ink, self.items = r, ink, []

    def path(self, path, fill=None, stroke=None, width=None):
        self.items.append(('path', path, fill, self.ink if stroke is None else stroke, width))

    def text(self, text, x, y, size, color=None, align='center', bold=False):
        self.items.append(('text', text, x, y, size, color or self.ink, align, bold))

    def paint(self, p, transform=None):
        import render
        r = self.r
        p.save()
        if transform is not None:
            p.setTransform(transform, True)
        lw = r.s['line_width']
        for it in self.items:
            if it[0] == 'path':
                _k, path, fill, stroke, width = it
                p.setBrush(QBrush(QColor(fill)) if fill else Qt.BrushStyle.NoBrush)
                if stroke == '':
                    p.setPen(Qt.PenStyle.NoPen)
                else:
                    p.setPen(QPen(QColor(stroke), lw if width is None else width, Qt.PenStyle.SolidLine,
                                  Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
                p.drawPath(path)
            else:
                _k, text, x, y, size, color, align, bold = it
                w = r.m.width(text, size, bold)
                x0 = x - w / 2 if align == 'center' else x - w if align == 'right' else x
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor(color)))
                p.drawPath(render.text_path(text, r.s['font'], size, x0, y + r.m.cap(size) / 2, bold=bold))
        p.restore()


def _shade(color, factor):
    c = QColor(color)
    return (c.darker(factor) if factor > 100 else c.lighter(int(10000 / factor))).name()


def paint_bio(r, p, i):
    sh = r.doc.shapes[i]
    fill = sh.fill or FILL.get(sh.kind, '#CCCCCC')
    ink = sh.color or INK
    ops = _Ops(r, ink)
    if sh.kind in model.BIO_LINEAR:
        t, L, H = _linear_frame(sh, r.BL)
        {'membrane': _membrane, 'dna': _dna, 'helix': _helix, 'mrna': _mrna}[sh.kind](ops, L, H, fill, sh)
        ops.paint(p, t)
        return
    b = _box(sh, r.BL)
    DRAW[sh.kind](ops, b, fill, sh)
    ops.paint(p)


def _lipid(ops, x, y, nx, ny, H, rh, fill, tails=True):
    """One lipid: head at (x, y) + normal * (H/2 - rh), two tails towards the middle."""
    hx, hy = x + nx * (H / 2 - rh), y + ny * (H / 2 - rh)
    head = QPainterPath(); head.addEllipse(QPointF(hx, hy), rh, rh)
    tail_end = 0.05 * H
    tx, ty = -ny, nx                                        # along the membrane
    path = QPainterPath()
    for s in (-0.45, 0.45):
        x0, y0 = hx - nx * rh + tx * s * rh, hy - ny * rh + ty * s * rh
        x1, y1 = x + nx * tail_end + tx * s * rh, y + ny * tail_end + ty * s * rh
        mx, my = (x0 + x1) / 2 + tx * s * rh * 0.5, (y0 + y1) / 2 + ty * s * rh * 0.5
        path.moveTo(x0, y0); path.quadTo(mx, my, x1, y1)
    if tails:
        ops.path(path, None, '#8A939E', max(0.5, rh * 0.22))
    ops.path(head, fill, _shade(fill, 135), max(0.4, rh * 0.15))


def _membrane(ops, L, H, fill, sh):
    rh = H * 0.11
    n = max(1, round(L / (rh * 2.3)))
    step = L / n
    for k in range(n):
        x = (k + 0.5) * step
        _lipid(ops, x, 0, 0, -1, H, rh, fill)
        _lipid(ops, x, 0, 0, 1, H, rh, fill)


def _dna(ops, L, H, fill, sh):
    period = H * 2.2
    amp = H * 0.4
    off = 0.36 * 2 * math.pi
    steps = max(24, int(L / 1.5))
    strand2 = sh.params.get('fill2', '#E07B53')
    rung_colours = ('#8CCB7E', '#F0B860', '#7FA7E0', '#E58C8A')
    k = 0
    u = period / 10
    while u < L:                                            # base pairs behind the strands
        a = 2 * math.pi * u / period
        y1, y2 = amp * math.sin(a), amp * math.sin(a + off)
        rung = QPainterPath(); rung.moveTo(u, y1); rung.lineTo(u, y2)
        ops.path(rung, None, rung_colours[k % 4], H * 0.07)
        k += 1
        u += period / 10
    for phase, colour in ((off, strand2), (0.0, fill)):
        path = QPainterPath()
        for j in range(steps + 1):
            u = L * j / steps
            y = amp * math.sin(2 * math.pi * u / period + phase)
            path.moveTo(u, y) if j == 0 else path.lineTo(u, y)
        ops.path(path, None, colour, H * 0.13)


def _helix(ops, L, H, fill, sh):
    s = H * 0.42
    n = max(1, int(L / s))
    s = L / n
    cw, ch = s * 0.62, H * 0.78
    loops = QPainterPath()
    for k in range(n - 1):                                  # loops between the cylinders, top and bottom in turn
        x0, x1 = (k + 0.5) * s, (k + 1.5) * s
        y = -ch / 2 if k % 2 == 0 else ch / 2
        bulge = -H * 0.18 if k % 2 == 0 else H * 0.18
        loops.moveTo(x0, y); loops.cubicTo(x0, y + bulge, x1, y + bulge, x1, y)
    ops.path(loops, None, _shade(fill, 150), H * 0.04)
    for k in range(n):
        x = (k + 0.5) * s
        cyl = QPainterPath(); cyl.addRoundedRect(QRectF(x - cw / 2, -ch / 2, cw, ch), cw / 2, cw / 2)
        ops.path(cyl, fill, _shade(fill, 140), H * 0.025)
        shine = QPainterPath(); shine.addRoundedRect(QRectF(x - cw * 0.28, -ch / 2 + cw * 0.3, cw * 0.18, ch - cw * 0.6),
                                                     cw * 0.09, cw * 0.09)
        ops.path(shine, _shade(fill, 70), '', None)


def _mrna(ops, L, H, fill, sh):
    amp, period = H * 0.28, H * 1.6
    cap = H * 0.32
    path = QPainterPath()
    steps = max(20, int(L / 1.5))
    start = cap * 2
    for j in range(steps + 1):
        u = start + (L - start) * j / steps
        y = amp * math.sin(2 * math.pi * (u - start) / period)
        path.moveTo(u, y) if j == 0 else path.lineTo(u, y)
    ops.path(path, None, fill, H * 0.16)
    c = QPainterPath(); c.addEllipse(QPointF(cap, 0), cap, cap)
    ops.path(c, _shade(fill, 130), _shade(fill, 170), H * 0.05)


def _micelle(ops, b, fill, sh):
    cx, cy = b.center().x(), b.center().y()
    R = min(b.width(), b.height()) / 2
    rh = R * 0.1
    H = R * 1.1                                             # one layer: heads out, tails reaching the middle
    n = max(6, int(2 * math.pi * (R - rh) / (rh * 2.3)))
    for k in range(n):
        a = 2 * math.pi * k / n
        nx, ny = math.cos(a), math.sin(a)
        x, y = cx + nx * (R - H / 2), cy + ny * (R - H / 2)
        _lipid(ops, x, y, nx, ny, H, rh, fill)


def _membrane_arc(ops, b, fill, sh):
    H = min(b.height() * 0.45, b.width() * 0.25)
    rh = H * 0.11
    cx, cy = b.center().x(), b.bottom() - H / 2
    rx, ry = b.width() / 2 - H / 2, b.height() - H
    length = math.pi * (rx + ry) / 2
    n = max(3, int(length / (rh * 2.3)))
    for k in range(n):
        a = math.pi + math.pi * (k + 0.5) / n               # the top half of the ellipse
        x, y = cx + rx * math.cos(a), cy + ry * math.sin(a)
        nx, ny = math.cos(a) / max(rx, 1e-6), math.sin(a) / max(ry, 1e-6)
        d = math.hypot(nx, ny) or 1.0
        _lipid(ops, x, y, nx / d, ny / d, H, rh, fill)
        _lipid(ops, x, y, -nx / d, -ny / d, H, rh, fill)


def _rounded(x, y, w, h, r=None):
    p = QPainterPath(); r = min(w, h) / 2 if r is None else r
    p.addRoundedRect(QRectF(x, y, w, h), r, r)
    return p


def _ellipse(x, y, w, h):
    p = QPainterPath(); p.addEllipse(QRectF(x, y, w, h))
    return p


def _receptor(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    cup = _rounded(x, y, w * 0.42, h * 0.5).united(_rounded(x + w * 0.58, y, w * 0.42, h * 0.5))
    cup = cup.united(_rounded(x, y + h * 0.3, w, h * 0.24, w * 0.12))
    stems = _rounded(x + w * 0.26, y + h * 0.45, w * 0.2, h * 0.55).united(_rounded(x + w * 0.54, y + h * 0.45, w * 0.2,
                                                                                    h * 0.55))
    ops.path(stems, _shade(fill, 120), _shade(fill, 160))
    ops.path(cup, fill, _shade(fill, 160))


def _channel(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    for k, (px, pw) in enumerate(((0.0, 0.2), (0.16, 0.22), (0.62, 0.22), (0.8, 0.2))):
        ops.path(_rounded(x + w * px, y + h * (0.04 if k in (1, 2) else 0.1), w * pw, h * (0.92 if k in (1, 2) else 0.8),
                          w * 0.09), fill if k in (1, 2) else _shade(fill, 120), _shade(fill, 160))
    for k in range(3):                                      # ions passing through
        ion = QPainterPath(); ion.addEllipse(QPointF(x + w * 0.5, y + h * (0.25 + 0.25 * k)), w * 0.05, w * 0.05)
        ops.path(ion, '#5BB5C8', '#2E7F8F', None)


def _enzyme(ops, b, fill, sh):
    body = _ellipse(b.x(), b.y(), b.width(), b.height())
    cx, cy = b.center().x(), b.center().y()
    notch = QPainterPath()
    notch.moveTo(cx + b.width() * 0.12, cy)
    notch.lineTo(b.right() + 2, cy - b.height() * 0.3); notch.lineTo(b.right() + 2, cy + b.height() * 0.3)
    notch.closeSubpath()
    ops.path(body.subtracted(notch), fill, _shade(fill, 160))


def _gprotein(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    parts = (('α', _ellipse(x, y + h * 0.18, w * 0.55, h * 0.82), fill),
             ('β', _ellipse(x + w * 0.46, y + h * 0.12, w * 0.4, h * 0.7), _shade(fill, 125)),
             ('γ', _ellipse(x + w * 0.8, y, w * 0.2, h * 0.5), _shade(fill, 150)))
    size = ops.r.s['font_size']
    for name, path, colour in parts:
        ops.path(path, colour, _shade(colour, 150))
        c = path.boundingRect().center()
        ops.text(name, c.x(), c.y(), size, '#FFFFFF' if name != 'α' else INK)


def _antibody(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    cx, cy = x + w / 2, y + h * 0.48
    arm_l, arm_w = min(w, h) * 0.55, min(w, h) * 0.11
    light = _shade(fill, 60)
    for k in (-1, 1):                                       # the stem (Fc): the two heavy chains
        ops.path(_rounded(cx + (0.05 * arm_w if k > 0 else -1.05 * arm_w), cy - arm_w * 0.3, arm_w,
                          y + h - cy + arm_w * 0.3), fill, _shade(fill, 160))
    for side in (-1, 1):                                    # the arms (Fab): heavy chain inside, light chain outside
        t = QTransform(); t.translate(cx, cy); t.rotate(side * 40)
        heavy = _rounded(-arm_w / 2, -arm_l, arm_w, arm_l)
        light_chain = _rounded(-arm_w / 2 + side * arm_w * 1.05, -arm_l, arm_w, arm_l * 0.62)
        ops.path(t.map(heavy), fill, _shade(fill, 160))
        ops.path(t.map(light_chain), light, _shade(fill, 160))


def _ribosome(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    ops.path(_ellipse(x + w * 0.1, y + h * 0.5, w * 0.8, h * 0.5), _shade(fill, 80), _shade(fill, 160))
    big = _ellipse(x, y, w, h * 0.66).united(_ellipse(x + w * 0.15, y - h * 0.02, w * 0.4, h * 0.3))
    ops.path(big, fill, _shade(fill, 160))


def _trna(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    cx, cy = x + w / 2, y + h * 0.5
    stroke = _shade(fill, 170)
    lw = max(1.0, w * 0.06)
    stem = QPainterPath(); stem.moveTo(cx, y); stem.lineTo(cx, cy)
    ops.path(stem, None, stroke, lw)
    loops = ((x + w * 0.17, cy, w * 0.17), (x + w * 0.83, cy, w * 0.17), (cx, y + h * 0.84, w * 0.18))
    for lx, ly, lr in loops:
        arm = QPainterPath(); arm.moveTo(cx, cy); arm.lineTo(lx, ly)
        ops.path(arm, None, stroke, lw)
    for lx, ly, lr in loops:
        c = QPainterPath(); c.addEllipse(QPointF(lx, ly), lr, lr)
        ops.path(c, fill, stroke, lw * 0.6)
    ops.text('CCA', cx + w * 0.08, y + h * 0.05, ops.r.s['font_size'] * 0.6, stroke, 'left')


def _mito(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    ops.path(_rounded(x, y, w, h), fill, _shade(fill, 160))
    inner = _rounded(x + h * 0.1, y + h * 0.1, w - h * 0.2, h * 0.8)
    ops.path(inner, _shade(fill, 85), _shade(fill, 140))
    folds = QPainterPath()                                  # cristae: the inner membrane folding in
    n = max(3, int(w / (h * 0.45)))
    for k in range(n):
        fx = x + h * 0.35 + (w - h * 0.7) * (k + 0.5) / n
        top = k % 2 == 0
        y0 = y + h * (0.12 if top else 0.88)
        y1 = y + h * (0.62 if top else 0.38)
        folds.moveTo(fx - h * 0.06, y0); folds.lineTo(fx - h * 0.06, y1)
        folds.arcTo(QRectF(fx - h * 0.06, y1 - h * 0.06, h * 0.12, h * 0.12), 180 if top else 0, 180)
        folds.lineTo(fx + h * 0.06, y0)
    ops.path(folds, None, _shade(fill, 150), max(0.8, h * 0.03))


def _nucleus(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    ops.path(_ellipse(x, y, w, h), fill, _shade(fill, 160), max(1.0, w * 0.025))
    ops.path(_ellipse(x + w * 0.04, y + h * 0.04, w * 0.92, h * 0.92), None, _shade(fill, 130), max(0.5, w * 0.01))
    ops.path(_ellipse(x + w * 0.52, y + h * 0.3, w * 0.24, h * 0.24), _shade(fill, 150), '', None)


def _cell(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    ops.path(_ellipse(x, y, w, h), fill, _shade(fill, 170), max(1.0, w * 0.01))
    _nucleus(ops, QRectF(x + w * 0.33, y + h * 0.3, w * 0.3, h * 0.36), FILL['nucleus'], sh)
    for fx, fy in ((0.18, 0.4), (0.75, 0.3), (0.72, 0.7), (0.25, 0.68), (0.5, 0.82)):
        ops.path(_ellipse(x + w * fx, y + h * fy, w * 0.07, h * 0.05), '#F4B183', '#C98B5E', None)


def _bacterium(ops, b, fill, sh):
    x, y, w, h = b.x(), b.y(), b.width(), b.height()
    body_w = w * 0.78
    flag = QPainterPath()
    for k in (-1, 0, 1):
        y0 = y + h / 2 + k * h * 0.18
        flag.moveTo(x + body_w - h * 0.1, y0)
        steps = 16
        for j in range(1, steps + 1):
            u = j / steps
            flag.lineTo(x + body_w - h * 0.1 + (w - body_w + h * 0.1) * u,
                        y0 + k * h * 0.15 * u + h * 0.1 * math.sin(u * 3 * math.pi))
    ops.path(flag, None, _shade(fill, 170), max(0.6, h * 0.025))
    ops.path(_rounded(x, y, body_w, h), fill, _shade(fill, 160))
    ops.path(_rounded(x + h * 0.12, y + h * 0.12, body_w - h * 0.24, h * 0.76), None, _shade(fill, 125), None)


def _virus(ops, b, fill, sh):
    cx, cy = b.center().x(), b.center().y()
    R = min(b.width(), b.height()) * 0.36
    spikes = QPainterPath()
    heads = QPainterPath()
    for k in range(12):
        a = 2 * math.pi * k / 12 + math.pi / 12
        spikes.moveTo(cx + R * 0.9 * math.cos(a), cy + R * 0.9 * math.sin(a))
        spikes.lineTo(cx + R * 1.28 * math.cos(a), cy + R * 1.28 * math.sin(a))
        heads.addEllipse(QPointF(cx + R * 1.3 * math.cos(a), cy + R * 1.3 * math.sin(a)), R * 0.1, R * 0.1)
    ops.path(spikes, None, _shade(fill, 160), max(0.8, R * 0.06))
    ops.path(heads, _shade(fill, 120), _shade(fill, 160))
    hexagon = QPainterPath()
    for k in range(7):
        a = math.pi / 6 + 2 * math.pi * k / 6
        pt = (cx + R * math.cos(a), cy + R * math.sin(a))
        hexagon.moveTo(*pt) if k == 0 else hexagon.lineTo(*pt)
    ops.path(hexagon, fill, _shade(fill, 160))
    inner = QPainterPath()
    for k in range(6):
        a = math.pi / 6 + 2 * math.pi * k / 6
        inner.moveTo(cx, cy); inner.lineTo(cx + R * math.cos(a), cy + R * math.sin(a))
    ops.path(inner, None, _shade(fill, 140), max(0.4, R * 0.03))


def _cloud(ops, b, fill, sh):
    cx, cy = b.center().x(), b.center().y()
    rx, ry = b.width() / 2, b.height() / 2
    path = _ellipse(cx - rx * 0.78, cy - ry * 0.7, rx * 1.56, ry * 1.4)
    n = 9
    for k in range(n):
        a = 2 * math.pi * k / n
        r = min(rx, ry) * (0.42 if k % 2 else 0.5)
        px, py = cx + (rx - r) * math.cos(a), cy + (ry - r) * math.sin(a)
        path = path.united(_ellipse(px - r, py - r, 2 * r, 2 * r))
    ops.path(path, fill, None)


def _plasmid_geometry(b):
    cx, cy = b.center().x(), b.center().y()
    return cx, cy, min(b.width(), b.height()) * 0.36


def _bp_angle(pos, bp):
    return -math.pi / 2 + 2 * math.pi * (pos - 1) / max(bp, 1)


def _plasmid(ops, b, fill, sh):
    cx, cy, R = _plasmid_geometry(b)
    bp = int(sh.params.get('bp', 5000))
    lw = max(1.2, R * 0.035)
    ring = QPainterPath(); ring.addEllipse(QPointF(cx, cy), R, R)
    ops.path(ring, None, fill, lw)
    size = ops.r.s['font_size']
    for k, reg in enumerate(sh.params.get('regions', [])):  # regions: coloured arcs with an arrowhead at the end
        a0, a1 = _bp_angle(reg['start'], bp), _bp_angle(reg['end'], bp)
        if a1 <= a0:
            a1 += 2 * math.pi
        colour = reg.get('color') or PLASMID_COLOURS[k % len(PLASMID_COLOURS)]
        band = R * 0.11
        head = min((a1 - a0) * 0.4, band * 1.4 / R)
        arc = QPainterPath()
        steps = max(4, int((a1 - a0 - head) / 0.05))
        for j in range(steps + 1):
            a = a0 + (a1 - head - a0) * j / steps
            pt = (cx + R * math.cos(a), cy + R * math.sin(a))
            arc.moveTo(*pt) if j == 0 else arc.lineTo(*pt)
        ops.path(arc, None, colour, band)
        tip = QPainterPath()
        ab = a1 - head
        tip.moveTo(cx + (R + band * 0.9) * math.cos(ab), cy + (R + band * 0.9) * math.sin(ab))
        tip.lineTo(cx + R * math.cos(a1), cy + R * math.sin(a1))
        tip.lineTo(cx + (R - band * 0.9) * math.cos(ab), cy + (R - band * 0.9) * math.sin(ab))
        tip.closeSubpath()
        ops.path(tip, colour, '', None)
        if reg.get('label'):
            am = (a0 + a1) / 2
            rr = R - band * 1.6
            ops.text(reg['label'], cx + rr * math.cos(am), cy + rr * math.sin(am), size * 0.75, _shade(colour, 140))
    for mk in sh.params.get('markers', []):                 # markers: a tick and a label outside the ring
        a = _bp_angle(mk['pos'], bp)
        tick = QPainterPath()
        tick.moveTo(cx + R * math.cos(a), cy + R * math.sin(a))
        tick.lineTo(cx + R * 1.12 * math.cos(a), cy + R * 1.12 * math.sin(a))
        ops.path(tick, None, INK, max(0.6, lw * 0.5))
        text = f"{mk.get('label', '')} ({mk['pos']})" if mk.get('label') else str(mk['pos'])
        rr = R * 1.17
        tx, ty = cx + rr * math.cos(a), cy + rr * math.sin(a)
        ops.text(text, tx, ty, size * 0.75, INK, 'left' if math.cos(a) > 0.2 else 'right' if math.cos(a) < -0.2
                 else 'center')
    name = sh.params.get('name', '')
    if name:
        ops.text(name, cx, cy - size * 0.7, size, INK, bold=True)
    ops.text(f'{bp:,} bp'.replace(',', ' '), cx, cy + (size * 0.6 if name else 0), size * 0.85, INK)


def _plasmid_bounds(r, sh, b):
    cx, cy, R = _plasmid_geometry(b)
    p = QPainterPath()
    rect = QRectF(cx - R * 1.15, cy - R * 1.15, R * 2.3, R * 2.3)
    size = r.s['font_size'] * 0.75
    bp = int(sh.params.get('bp', 5000))
    for mk in sh.params.get('markers', []):
        a = _bp_angle(mk['pos'], bp)
        text = f"{mk.get('label', '')} ({mk['pos']})" if mk.get('label') else str(mk['pos'])
        w = r.m.width(text, size)
        tx, ty = cx + R * 1.17 * math.cos(a), cy + R * 1.17 * math.sin(a)
        x0 = tx if math.cos(a) > 0.2 else tx - w if math.cos(a) < -0.2 else tx - w / 2
        rect = rect.united(QRectF(x0, ty - size, w, size * 2))
    p.addRect(rect.united(b))
    return p


DRAW = {'membrane_arc': _membrane_arc, 'micelle': _micelle, 'receptor': _receptor, 'channel': _channel,
        'enzyme': _enzyme, 'gprotein': _gprotein, 'antibody': _antibody, 'ribosome': _ribosome, 'trna': _trna,
        'mito': _mito, 'nucleus': _nucleus, 'cell': _cell, 'bacterium': _bacterium, 'virus': _virus, 'cloud': _cloud,
        'plasmid': _plasmid}


# ------------------------------------------------------------------------------------------------ sequences
def seq_layout(r, i):
    """Geometry of sequence i: {'cells': [QRectF of each residue], 'rect': QRectF around everything, 'd': monomer
    size, 'left' / 'right': (text, x, y) of the termini, 'top': y above which the bridges arc} - all in pt."""
    q = r.doc.seqs[i]
    fs = r.s['font_size']
    d = fs * 1.9
    gap = fs * 0.55
    three = q.polymer == 'PEPTIDE' and q.three_letter
    cell_w = d
    if three:
        cell_w = max(d, max([r.m.width(biopolymer.label(q, k), fs * 0.8) for k in range(len(q.residues))] or [0])
                     + fs * 0.9)
    if q.polymer == 'PEPTIDE':
        left = '' if q.cyclic else ('Ac' if q.n_term == 'Ac' else 'H')
        right = '' if q.cyclic else ('NH₂' if q.c_term == 'NH2' else 'OH')
    else:
        left = '5′' + ('-p' if q.p5 else '')
        right = ('p-' if q.p3 else '') + '3′'
    term = fs * 0.85
    lead = (r.m.width(left, term) + gap) if left else 0.0
    id_h = fs * 0.95 if q.show_ids else 0.0
    bridge_h = d * (0.5 + 0.22 * len(q.bridges)) if q.bridges else 0.0
    per = max(1, int(q.per_line or 10))
    row_h = id_h + d + gap * 1.8
    x0, y0 = q.x * r.BL, q.y * r.BL + bridge_h
    cells = []
    for k in range(len(q.residues)):
        row, col = divmod(k, per)
        cells.append(QRectF(x0 + lead + col * (cell_w + gap), y0 + row * row_h + id_h, cell_w, d))
    rect = QRectF(q.x * r.BL, q.y * r.BL, max(lead, 1.0), max(d + bridge_h, 1.0))
    for c in cells:
        rect = rect.united(c)
    lt = rt = None
    if cells:
        if left:
            lt = (left, cells[0].left() - gap * 0.6, cells[0].center().y())
        if right:
            rt = (right, cells[-1].right() + gap * 0.6, cells[-1].center().y())
            rect = rect.united(QRectF(rt[1], rt[2] - term, r.m.width(right, term) + 1, 2 * term))
        if q.cyclic and len(cells) > 1:
            rect = rect.united(QRectF(cells[-1].left(), cells[-1].bottom(), 1, d * 0.9))
    if id_h:
        rect = rect.adjusted(0, 0, 0, 0)
    return {'cells': cells, 'rect': rect, 'd': d, 'gap': gap, 'left': lt, 'right': rt, 'term': term,
            'per': per, 'top': y0}


def residue_at(r, i, x_pt, y_pt):
    """Index of the residue of sequence i under the point (pt), or None."""
    for k, c in enumerate(seq_layout(r, i)['cells']):
        if c.adjusted(-1, -1, 1, 1).contains(QPointF(x_pt, y_pt)):
            return k
    return None


def paint_seq(r, p, i):
    q = r.doc.seqs[i]
    L = seq_layout(r, i)
    cells, d, gap = L['cells'], L['d'], L['gap']
    ink = q.color or INK
    ops = _Ops(r, ink)
    fs = r.s['font_size']
    lw = r.s['line_width'] * 1.4
    links = QPainterPath()
    for k in range(len(cells) - 1):
        a, b = cells[k], cells[k + 1]
        if abs(a.center().y() - b.center().y()) < 1e-6:
            links.moveTo(a.right(), a.center().y()); links.lineTo(b.left(), b.center().y())
        else:                                               # wrap: round the end of the line to the next one
            yb = a.bottom() + gap * 0.9
            links.moveTo(a.right(), a.center().y())
            links.lineTo(a.right() + gap * 0.5, a.center().y())
            links.lineTo(a.right() + gap * 0.5, yb)
            links.lineTo(b.left() - gap * 0.5, yb)
            links.lineTo(b.left() - gap * 0.5, b.center().y())
            links.lineTo(b.left(), b.center().y())
    for t in (L['left'], L['right']):
        if t:
            c = cells[0] if t is L['left'] else cells[-1]
            x_edge = c.left() if t is L['left'] else c.right()
            links.moveTo(x_edge, c.center().y()); links.lineTo(t[1], c.center().y())
    if q.cyclic and len(cells) > 1:                         # head to tail, under the sequence
        a, b = cells[-1], cells[0]
        low = max(c.bottom() for c in cells) + d * 0.7
        links.moveTo(a.center().x(), a.bottom())
        links.cubicTo(a.center().x(), low, b.center().x(), low, b.center().x(), b.bottom())
    ops.path(links, None, _shade(ink, 60) if not q.color else ink, lw)
    for n, (bi, bj) in enumerate(q.bridges):                # disulfide / lactam bridges, arcing above
        if not (0 <= bi < len(cells) and 0 <= bj < len(cells)):
            continue
        a, b = cells[bi], cells[bj]
        ss = q.residues[bi] in 'CU' and q.residues[bj] in 'CU'
        colour = '#C9A400' if ss else '#7A5CC4'
        arc = QPainterPath()
        if abs(a.center().y() - b.center().y()) < 1e-6:     # same line: an arc above it
            top = a.top() - d * (0.35 + 0.22 * n)
            arc.moveTo(a.center().x(), a.top())
            arc.cubicTo(a.center().x(), top, b.center().x(), top, b.center().x(), b.top())
            lx, ly = (a.center().x() + b.center().x()) / 2, top + d * 0.05
        else:                                               # lines apart: from under the upper to over the lower
            up, lo = (a, b) if a.top() < b.top() else (b, a)
            mid = (up.bottom() + lo.top()) / 2
            arc.moveTo(up.center().x(), up.bottom())
            arc.cubicTo(up.center().x(), mid, lo.center().x(), mid, lo.center().x(), lo.top())
            lx, ly = (up.center().x() + lo.center().x()) / 2 + fs * 0.9, mid
        ops.path(arc, None, colour, lw * 1.2)
        ops.text('S–S' if ss else 'lactam', lx, ly, fs * 0.55, colour)
    for k, c in enumerate(cells):
        res = q.residues[k]
        colour = biopolymer.colour(q.polymer, res)
        if q.polymer == 'PEPTIDE':
            shape = _rounded(c.x(), c.y(), c.width(), c.height())
        else:
            shape = _rounded(c.x(), c.y(), c.width(), c.height(), c.height() * 0.22)
        ops.path(shape, colour, _shade(colour, 150), r.s['line_width'])
        text = biopolymer.label(q, k)
        size = fs * (0.8 if len(text) > 1 else 0.95)
        ops.text(text, c.center().x(), c.center().y(), size, '#1F2328', bold=len(text) == 1)
        if q.show_ids:
            ops.text(str(k + 1), c.center().x(), c.top() - fs * 0.5, fs * 0.6, '#6A7380')
    for t, align in ((L['left'], 'right'), (L['right'], 'left')):
        if t:
            ops.text(t[0], t[1], t[2], L['term'], ink, align)
    ops.paint(p)
