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
                     font='Arial', font_size=10.0, label_bold=False, caption_font='Arial', caption_size=10.0,
                     color='#000000', color_atoms=False, page='A4', orientation='portrait',
                     show_stereo=False, show_atom_numbers=False, show_terminal_carbons=False)
ACS_STYLE = dict(DEFAULT_STYLE)
# Drawing settings of the journal / document styles ChemDraw ships (numbers read from its style sheets): bond length,
# line, bold and hash widths, margin (pt), double-bond spacing (fraction of the bond), label and caption fonts.
STYLE_SHEETS = {
    'ACS Document 1996': dict(bond_length=14.4, line_width=0.6, bold_width=2.0, hash_spacing=2.5, margin=1.6,
                              bond_spacing=0.18, font='Arial', font_size=10, caption_font='Arial', caption_size=10),
    'New Document': dict(bond_length=30.0, line_width=1.0, bold_width=4.0, hash_spacing=2.7, margin=2.0,
                         bond_spacing=0.12, font='Arial', font_size=10, caption_font='Times New Roman', caption_size=12),
    'New Slide': dict(bond_length=30.0, line_width=1.6, bold_width=4.0, hash_spacing=2.7, margin=2.0,
                      bond_spacing=0.15, font='Arial', font_size=12, caption_font='Times New Roman', caption_size=16),
    'RSC Document': dict(bond_length=12.2, line_width=0.45, bold_width=1.6, hash_spacing=1.75, margin=1.25,
                         bond_spacing=0.20, font='Arial', font_size=7, caption_font='Arial', caption_size=7),
    'Wiley Document': dict(bond_length=17.0, line_width=0.75, bold_width=2.6, hash_spacing=2.6, margin=2.0,
                           bond_spacing=0.18, font='Arial', font_size=12, caption_font='Arial', caption_size=12),
    'Nature Research Document': dict(bond_length=10.8, line_width=0.595, bold_width=1.56, hash_spacing=1.7,
                                     margin=1.19, bond_spacing=0.18, font='Arial', font_size=6, caption_font='Arial',
                                     caption_size=6),
    'SYNTHESIS, SYNLETT Document': dict(bond_length=17.0, line_width=0.8, bold_width=2.0, hash_spacing=2.5, margin=1.3,
                                        bond_spacing=0.18, font='Arial', font_size=10, caption_font='Arial',
                                        caption_size=10),
    'Science of Synthesis': dict(bond_length=17.0, line_width=0.8, bold_width=2.0, hash_spacing=2.5, margin=2.25,
                                 bond_spacing=0.18, font='Arial', font_size=10, caption_font='Arial', caption_size=10),
    'Adv. Synth. Catal. Document': dict(bond_length=17.0, line_width=1.0, bold_width=2.0, hash_spacing=2.5,
                                        margin=1.6, bond_spacing=0.18, font='Arial', font_size=10,
                                        caption_font='Times New Roman', caption_size=12),
    'Helv. Chim. Acta Document': dict(bond_length=17.0, line_width=0.6, bold_width=2.0, hash_spacing=2.0, margin=2.0,
                                      bond_spacing=0.15, font='Arial', font_size=10, caption_font='Arial',
                                      caption_size=12),
    'J. Chin. Chem. Soc. Document': dict(bond_length=18.0, line_width=1.0, bold_width=2.5, hash_spacing=2.5,
                                         margin=2.0, bond_spacing=0.20, font='Arial', font_size=12,
                                         caption_font='Arial', caption_size=12),
    'J. Het. Chem.': dict(bond_length=10.8, line_width=0.605, bold_width=2.06, hash_spacing=1.3, margin=1.08,
                          bond_spacing=0.16, font='Times New Roman', font_size=7, caption_font='Times New Roman',
                          caption_size=7),
    'J. Mol. Mod.': dict(bond_length=14.4, line_width=0.6, bold_width=2.0, hash_spacing=2.5, margin=1.6,
                         bond_spacing=0.18, font='Times New Roman', font_size=10, caption_font='Times New Roman',
                         caption_size=10),
    'Phytomedicine Document': dict(bond_length=20.0, line_width=1.0, bold_width=1.33, hash_spacing=3.0, margin=1.25,
                                   bond_spacing=0.08, font='Arial', font_size=12, label_bold=True,
                                   caption_font='Arial', caption_size=12),
}
for _sheet in STYLE_SHEETS.values():                       # what a sheet does not set goes back to plain defaults
    _sheet.setdefault('label_bold', False)
PAGES = {'A0': (2383.94, 3370.39), 'A1': (1683.78, 2383.94), 'A2': (1190.55, 1683.78), 'A3': (841.89, 1190.55),
         'A4': (595.28, 841.89), 'A5': (419.53, 595.28), 'B5': (498.9, 708.66),
         'Letter': (612.0, 792.0), 'Legal': (612.0, 1008.0), 'Tabloid': (792.0, 1224.0)}    # pt, portrait


def page_size(style):
    """(width, height) in pt of the page a style asks for, turned for landscape."""
    w, h = PAGES.get(style.get('page', 'A4'), PAGES['A4'])
    return (h, w) if style.get('orientation') == 'landscape' else (w, h)
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
    def __init__(self, family, bold=False):
        self.family = family
        self.bold = bold                                    # the default weight (bold atom labels in some styles)
        self.fms = {}

    def fm(self, bold=False, italic=False):
        if (bold, italic) not in self.fms:
            self.fms[(bold, italic)] = QFontMetricsF(_font(self.family, bold, italic))
        return self.fms[(bold, italic)]

    def width(self, text, size, bold=None, italic=False):
        bold = self.bold if bold is None else bold
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
        self.m = Metrics(self.s['font'], bold=bool(self.s['label_bold']))     # atom labels
        self.mc = Metrics(self.s['caption_font'])                          # captions
        self._first = True
        self.labels = {a: self._label(a) for a in doc.atoms}
        self._first = False
        for a, lab in list(self.labels.items()):            # again, now that the other labels are known: the H, the
            if lab and len(lab['pieces']) > 1:              # sign and a nickname's reading way go where they are clear
                self.labels[a] = self._label(a)
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
            def build(toks, anchor):
                widths = []
                for t in toks:
                    letters = t.rstrip('0123456789')
                    widths.append((m.width(letters, fs), m.width(t[len(letters):], small) if t != letters else 0.0,
                                   letters))
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
                if at.charge or at.radical:
                    pieces.append((self._sign_text(at), x + 0.3, base - ch + 0.15 * fs, small))
                rect = QRectF(x0, base - ch, x - x0, ch + 0.2 * small)
                return dict(pieces=pieces, clip=(aw / 2 + s['margin'], ch / 2 + s['margin']), rect=rect)
            best = build(*model.label_order(at.label, from_right))
            if not self._first and self._collisions(a, best['pieces']):
                other = model.label_order(at.label, not from_right)    # read the other way if that is clear
                if other[0] != model.label_order(at.label, from_right)[0]:
                    alt = build(*other)
                    if self._collisions(a, alt['pieces']) < self._collisions(a, best['pieces']):
                        best = alt
            return best
        if (at.el == 'C' and nb and not at.charge and not at.radical and at.hcount is None and not at.iso
                and not (s['show_terminal_carbons'] and len(nb) == 1)):
            return None                                     # a plain vertex (a charged carbon or radical shows C, its
                                                            # H and its charge / dot: CH2•, CH3+)
        sym = at.el
        mass = ''                                           # a labelled isotope: its mass number up before the symbol
        if at.iso and model.isotope_mass(at.el, at.iso):
            if at.el == 'H' and at.iso in (2, 3):
                sym = 'DT'[at.iso - 2]
            else:
                mass = str(at.iso)
        hs = doc.implicit_h(a)
        ws = m.width(sym, fs)
        base = cy + ch / 2
        x0 = cx - ws / 2
        pieces = [(sym, x0, base, fs)]
        wm = m.width(mass, small) + 0.2 if mass else 0.0
        if mass:
            pieces.append((mass, x0 - wm, base - ch + 0.32 * small, small))
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
            side = self._h_side(a, side, right, left)

            def h_pieces(side):
                if side == 'right':
                    out = [('H', right_end, base, fs)]
                    if hcount:
                        out.append((hcount, right_end + wh, base + 0.28 * fs, small))
                elif side == 'left':
                    hx = x0 - wm - wh - wc                  # H3¹³C: the mass number stays on its element
                    out = [('H', hx, base, fs)]
                    if hcount:
                        out.append((hcount, hx + wh, base + 0.28 * fs, small))
                else:
                    hy = base + (ch + 1.4 * s['margin'] + 0.5) if side == 'down' else base - (ch + 1.4 * s['margin'] + 0.5)
                    out = [('H', cx - wh / 2, hy, fs)]
                    if hcount:
                        out.append((hcount, cx + wh / 2, hy + 0.28 * fs, small))
                return out
            chosen = h_pieces(side)
            if not self._first and self._collisions(a, chosen):   # the usual side is taken: the clearest open one
                ups = [doc.atoms[n].y - at.y for n, _b in nb]
                open_ = [x for x, ok in (('right', not right), ('left', not left),
                                         ('down', not any(d > 0.5 for d in ups)), ('up', not any(d < -0.5 for d in ups)))
                         if ok and x != side]
                for alt in open_:
                    cand = h_pieces(alt)
                    if self._collisions(a, cand) < self._collisions(a, chosen):
                        chosen, side = cand, alt
            pieces += chosen
            if side == 'right':
                right_end += wh + wc
        charge_rect = None
        if at.charge or at.radical:
            piece = self._place_charge(a, at.charge, pieces, x0, ws, base, ch, fs, small, right_end)
            pieces.append(piece)
            charge_rect = QRectF(piece[1], piece[2] - m.cap(small), m.width(piece[0], small), m.cap(small))
        rect = QRectF()
        for text, x, y, size in pieces:
            r = QRectF(x, y - m.cap(size), m.width(text, size), m.cap(size) + (0.2 * size if size == small else 0))
            rect = r if rect.isNull() else rect.united(r)
        avoid = []                                          # ¹⁵N: a bond coming at the 15 stops before it
        if mass:
            mp = pieces[1]
            avoid.append(QRectF(mp[1] - cx, mp[2] - m.cap(small) - cy, wm, m.cap(small)).adjusted(
                -s['margin'], -s['margin'], s['margin'], s['margin']))
        return dict(pieces=pieces, clip=(ws / 2 + s['margin'], ch / 2 + s['margin']), rect=rect, avoid=avoid,
                    charge_rect=charge_rect)

    def _collisions(self, a, pieces):
        """How badly label pieces of atom a sit on other atoms' letters or on bonds (its own bonds near it excepted)."""
        doc, m, BL = self.doc, self.m, self.BL
        at = doc.atoms[a]
        rects = [QRectF(x, y - m.cap(size), m.width(t, size), m.cap(size)).adjusted(0.3, 0.3, -0.3, -0.3)
                 for t, x, y, size in pieces if t.strip()]
        n = 0
        for o, lab in self.labels.items():
            if o == a or not lab:
                continue
            ot = doc.atoms[o]
            if abs(ot.x - at.x) > 4 or abs(ot.y - at.y) > 3:
                continue
            for t, x, y, size in lab['pieces']:
                if t.strip():
                    r2 = QRectF(x, y - m.cap(size), m.width(t, size), m.cap(size))
                    n += sum(1 for r in rects if r.intersects(r2))
        for b, bd in doc.bonds.items():
            p, q = doc.atoms[bd.a], doc.atoms[bd.b]
            if min(abs(p.x - at.x) + abs(p.y - at.y), abs(q.x - at.x) + abs(q.y - at.y)) > 5:
                continue
            own = a in (bd.a, bd.b)
            for k in range(1, 20):
                px, py = p.x + (q.x - p.x) * k / 20, p.y + (q.y - p.y) * k / 20
                if own and (px - at.x) ** 2 + (py - at.y) ** 2 < 0.45 ** 2:
                    continue                                # the bond leaves from under its own letters
                pt = QPointF(px * BL, py * BL)
                if any(r.contains(pt) for r in rects):
                    n += 1
                    break
        return n

    def _h_side(self, a, side, right, left):
        """The side for an atom's H: the usual one, unless another atom sits there (OH beside HO in a diphosphate) -
        then the other side, below or above, whichever has room and no bond of its own."""
        doc = self.doc
        at = doc.atoms[a]
        mine = {n for n, _b in doc.neighbors(a)}
        others = [(o.x - at.x, o.y - at.y) for i, o in doc.atoms.items() if i != a and i not in mine
                  and abs(o.x - at.x) < 1.6 and abs(o.y - at.y) < 1.4]
        if not others:
            return side
        ups = [doc.atoms[n].y - at.y for n in mine]
        open_ = {'right': not right, 'left': not left,
                 'down': not any(d > 0.5 for d in ups), 'up': not any(d < -0.5 for d in ups)}

        def crowded(s):
            for dx, dy in others:
                if s == 'right' and 0.2 < dx < 1.5 and abs(dy) < 0.6:
                    return True
                if s == 'left' and -1.5 < dx < -0.2 and abs(dy) < 0.6:
                    return True
                if s == 'down' and 0.2 < dy < 1.2 and abs(dx) < 0.6:
                    return True
                if s == 'up' and -1.2 < dy < -0.2 and abs(dx) < 0.6:
                    return True
            return False
        if not crowded(side):
            return side
        for s in ('right', 'left', 'down', 'up'):
            if s != side and open_[s] and not crowded(s):
                return s
        return side

    def charge_rect(self, a):
        """QRectF (pt) of the atom's charge sign, or None (what the mouse grabs to move it)."""
        lab = self.labels.get(a)
        return lab.get('charge_rect') if lab else None

    def _place_charge(self, a, z, pieces, x0, ws, base, ch, fs, small, right_end):
        """The charge as a label piece, on the atom's own symbol: the first of top right, bottom right, top left,
        bottom left, right, left (ChemDraw's order) that is clear of the bonds and of the atom's H."""
        doc, m, BL = self.doc, self.m, self.BL
        at = doc.atoms[a]
        cx, cy = at.x * BL, at.y * BL
        t = self._sign_text(at)
        w, hgt = m.width(t, small), m.cap(small)
        if at.charge_pos:                                   # dragged by hand: centred where it was dropped
            dx, dy = at.charge_pos
            return (t, cx + dx * BL - w / 2, cy + dy * BL + hgt / 2, small)
        gap = 0.3
        row_right = right_end if any(abs(p[2] - base) < 0.1 for p in pieces[1:]) else x0 + ws
        top, bottom, mid = base - ch + 0.15 * fs, base + 0.35 * fs, base - ch / 2 + hgt / 2
        left = min(p[1] for p in pieces if abs(p[2] - base) < 0.1 or p[2] > base) - w - gap
        candidates = [(row_right + gap, top), (x0 + ws + gap, bottom), (left, top), (x0 - w - gap, bottom),
                      (row_right + gap, mid), (left, mid)]
        others = [QRectF(x, y - m.cap(size), m.width(txt, size), m.cap(size)) for txt, x, y, size in pieces]
        hw, hh = ws / 2 + self.s['margin'], ch / 2 + self.s['margin']
        segments = []                                       # the drawn lines of each bond, from the label's edge
        sp = self.s['bond_spacing'] * BL
        for n, b in doc.neighbors(a):
            o = doc.atoms[n]
            dx, dy = o.x * BL - cx, o.y * BL - cy
            L = (dx * dx + dy * dy) ** 0.5 or 1.0
            start = max(0.0, exit_t(0, 0, dx / L, dy / L, hw, hh) - 1.0)
            order = doc.bonds[b].order
            offsets = [0.0] if order == 1 else [-sp, 0.0, sp] if order == 3 else [-sp, -sp / 2, 0.0, sp / 2, sp]
            nx, ny = -dy / L, dx / L
            for off in offsets:                             # second and third lines too (either side, to be safe)
                segments.append([(cx + nx * off + dx * (start + (L - start) * k / 40) / L,
                                  cy + ny * off + dy * (start + (L - start) * k / 40) / L) for k in range(41)])

        mine = {n for n, _b in doc.neighbors(a)}
        for b, bd in doc.bonds.items():                     # other bonds nearby (a ring closing past the atom ...)
            if a in (bd.a, bd.b):
                continue
            p, q = doc.atoms[bd.a], doc.atoms[bd.b]
            if min(math.hypot(p.x - at.x, p.y - at.y), math.hypot(q.x - at.x, q.y - at.y)) > 2.2:
                continue
            segments.append([(p.x * BL + (q.x - p.x) * BL * k / 40, p.y * BL + (q.y - p.y) * BL * k / 40)
                             for k in range(41)])
        for o, lab in getattr(self, 'labels', {}).items():  # other atoms' letters (known on the second pass)
            if o != a and lab and abs(doc.atoms[o].x - at.x) < 2.5 and abs(doc.atoms[o].y - at.y) < 2.5:
                others += [QRectF(px, py - m.cap(size), m.width(txt, size), m.cap(size))
                           for txt, px, py, size in lab['pieces'] if txt.strip()]

        def clear(x, y):
            r = QRectF(x, y - hgt, w, hgt)
            if any(r.adjusted(0.1, 0.1, -0.1, -0.1).intersects(o) for o in others):   # not on the label's letters
                return False
            near = r.adjusted(-1.4, -1.4, 1.4, 1.4)                                 # a visible gap to any bond
            return not any(near.contains(QPointF(px, py)) for seg in segments for px, py in seg)
        for x, y in candidates:
            if clear(x, y):
                return (t, x, y, small)
        # crowded (four bonds, a ring and an H ...): the middle of the widest gap between the bonds and the H
        angles = [math.atan2(doc.atoms[n].y - at.y, doc.atoms[n].x - at.x) for n, _b in doc.neighbors(a)]
        for txt, px_, py_, _size in pieces[1:]:
            if txt == 'H':
                angles.append(math.atan2(py_ - ch / 2 - cy, px_ + m.width('H', fs) / 2 - cx))
        angles.sort()
        gaps = [((angles[(k + 1) % len(angles)] - angles[k]) % (2 * math.pi) or 2 * math.pi, angles[k])
                for k in range(len(angles))]
        size, start = max(gaps)
        ang = start + size / 2
        dist = max(hw, hh) + max(w, hgt) * 0.4
        mx, my = cx + math.cos(ang) * dist, cy + math.sin(ang) * dist
        return (t, mx - w / 2, my + hgt / 2, small)

    @classmethod
    def _sign_text(cls, at):
        """What goes beside the atom: its radical dots and its charge (•, ••, +, •+, 2−)."""
        return '•' * at.radical + (cls._charge_text(at.charge) if at.charge else '')

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
    def _end(self, a, ox, oy, ux, uy, avoid=True):
        """Distance along (ux, uy) from atom a (offset ox, oy) to where a line may start (0 = unlabelled atom)."""
        lab = self.labels.get(a)
        if not lab or not lab['clip']:
            return 0.0
        t = exit_t(ox, oy, ux, uy, *lab['clip'])
        for r in lab.get('avoid', ()) if avoid else ():     # a box the line must also clear (an isotope's mass)
            lo, hi = 0.0, math.inf
            for o, u, a0, a1 in ((ox, ux, r.left(), r.right()), (oy, uy, r.top(), r.bottom())):
                if abs(u) < 1e-9:
                    if not a0 <= o <= a1:
                        lo, hi = 1.0, 0.0
                    continue
                t0, t1 = sorted(((a0 - o) / u, (a1 - o) / u))
                lo, hi = max(lo, t0), min(hi, t1)
            if lo <= hi and hi > t:
                t = hi
        return t

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
        if t1 + t2 >= L - 2.5:                              # crowded by an isotope's mass number: rather a bond that
            t1 = max(self._end(a, ox, oy, ux, uy, False), trim_a)     # meets it than no bond at all
            t2 = max(self._end(b, ox, oy, -ux, -uy, False), trim_b)
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
        fs = t.size or self.s['caption_size']
        small = fs * SUB
        m = self.mc
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
        if sh.kind in model.BIO_KINDS:
            import biodraw
            return biodraw.outline(self, sh)
        if sh.kind in model.ELEMENT_KINDS or sh.style:
            import elements
            return elements.outline(self, sh)
        x1, x2 = sorted((sh.x1 * BL, sh.x2 * BL))
        y1, y2 = sorted((sh.y1 * BL, sh.y2 * BL))
        p = QPainterPath()
        w, h = x2 - x1, y2 - y1
        if sh.kind == 'line':
            p.moveTo(sh.x1 * BL, sh.y1 * BL); p.lineTo(sh.x2 * BL, sh.y2 * BL)
        elif sh.kind in ('rect', 'image'):
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

    def repeat_path(self, i):
        """The repeat count of polymer brackets (n, 3 ...), small at the foot of the right bracket; None without."""
        sh = self.doc.shapes[i]
        text = str(sh.params.get('repeat', '')).strip() if sh.kind in model.BRACKET_KINDS else ''
        if not text:
            return None
        size = self.s['font_size'] * 0.85
        x = max(sh.x1, sh.x2) * self.BL + 1.2
        y = max(sh.y1, sh.y2) * self.BL + self.m.cap(size) * 0.45
        return text_path(text, self.s['font'], size, x, y, italic=not text.isdigit())

    # ------------------------------------------------------------------------------------------- stoichiometry grid
    GRID_COLORS = {'border': '#C5CED8', 'head': '#EEF2F7', 'alt': '#F6F8FA', 'rule': '#E3E8EE', 'label': '#5B6673',
                   'value': '#1F2933', 'typed': '#1659C7', 'unit': '#7B8794', 'limit_fill': '#FFF4D6',
                   'limit': '#B7791F', 'reactants': ('#E6EEF9', '#1F4E99'), 'products': ('#E5F4EA', '#1E7A45')}

    def grid_layout(self, i):
        """What a stoichiometry grid draws: ([('text', text, x, baseline, size, bold, italic, colour) |
        ('fill', QRectF, colour, radius) | ('frame', QRectF, colour, radius) | ('line', x1, y1, x2, y2, colour)],
        QRectF, {(side, column index, row key | 'header'): QRectF of the cell, ('title', side): QRectF}).
        A card per side: a coloured title, a header row with each compound's formula (the limiting reactant
        marked), then one row per quantity - typed values in blue, calculated ones in ink, units in grey."""
        import stoich
        from PyQt6.QtCore import QLocale
        g = self.doc.grids[i]
        C = self.GRID_COLORS
        own = g.color or ''
        table = stoich.compute(self.doc, g)
        dec = QLocale.system().decimalPoint() or ','
        fs = self.s['font_size']
        small = fs * 0.72
        m = self.m
        row_h = fs * 1.55
        head_h = fs * 2.3
        pad = fs * 0.55
        ops, cells = [], {}
        x = g.x * self.BL
        top = g.y * self.BL
        bounds = QRectF()

        def text_w(t, size, bold=False, italic=False, formula=False):
            return sum(m.width(txt, size * (SUB if sub else 1), bold, italic) for txt, sub in runs(t, formula))

        def put_text(t, x0, base, size, color, bold=False, italic=False, formula=False):
            for txt, sub in runs(t, formula):
                sz = size * (SUB if sub else 1)
                ops.append(('text', txt, x0, base + (0.28 * size if sub else 0), sz, bold, italic, color))
                x0 += m.width(txt, sz, bold, italic)
            return x0
        for side, rows in (('reactants', stoich.REACTANT_ROWS), ('products', stoich.PRODUCT_ROWS)):
            cols = table[side]
            if not cols:
                continue
            rows = [(k, label) for k, label in rows if k not in g.hidden and k not in ('formula', 'limiting')]
            inputs = stoich.REACTANT_INPUTS if side == 'reactants' else stoich.PRODUCT_INPUTS
            texts = []                                      # per column, per row: (number, unit, typed)
            for col in cols:
                col_texts = {}
                for key, _label in rows:
                    cell = col.cells.get(key)
                    num, unit = stoich.display(key, cell, dec) if cell is not None else ('', '')
                    col_texts[key] = (num, unit, bool(cell and cell.typed))
                texts.append(col_texts)
            label_w = max(m.width(label, small) for _k, label in rows) + 2 * pad if rows else fs * 4
            col_ws = []
            for col, col_texts in zip(cols, texts):
                w = max([m.width(n, fs, t) + (m.width(u, small) + fs * 0.2 if u else 0)
                         for n, u, t in col_texts.values()] + [text_w(col.formula, fs, True, formula=True),
                                                               m.width('limiting', small), fs * 3.6])
                col_ws.append(w + 2 * pad)
            w_total = label_w + sum(col_ws)
            h_total = head_h + row_h * len(rows)
            fill, ink = C[side]
            title = 'Reactants' if side == 'reactants' else 'Products'
            tw = m.width(title, small, True) + 2 * pad
            chip = QRectF(x, top - fs * 1.55, tw, fs * 1.2)
            ops.append(('fill', chip, fill, fs * 0.6))
            put_text(title, x + pad, chip.bottom() - fs * 0.33, small, own or ink, bold=True)
            cells[('title', side)] = chip
            card = QRectF(x, top, w_total, h_total)
            ops.append(('fill', card, '#FFFFFF', fs * 0.45))
            ops.append(('fill', QRectF(x, top, w_total, head_h), C['head'], fs * 0.45))
            ops.append(('fill', QRectF(x, top + head_h - fs * 0.5, w_total, fs * 0.5), C['head'], 0))
            for r in range(len(rows)):                      # zebra rows
                if r % 2:
                    ops.append(('fill', QRectF(x, top + head_h + r * row_h, w_total, row_h), C['alt'],
                                fs * 0.45 if r == len(rows) - 1 else 0))
            cx = x + label_w
            for k, (col, cw) in enumerate(zip(cols, col_ws)):
                limiting = side == 'reactants' and col.cells.get('limiting') and col.cells['limiting'].value
                if limiting:
                    ops.append(('fill', QRectF(cx, top, cw, head_h), C['limit_fill'], 0))
                fw = text_w(col.formula, fs, True, formula=True)
                put_text(col.formula, cx + (cw - fw) / 2, top + head_h * (0.5 if limiting else 0.62), fs,
                         own or C['value'], bold=True, formula=True)
                if limiting:
                    lw = m.width('limiting', small, False, True)
                    put_text('limiting', cx + (cw - lw) / 2, top + head_h * 0.86, small, C['limit'], italic=True)
                cells[(side, k, 'header')] = QRectF(cx, top, cw, head_h)
                for r, (key, _label) in enumerate(rows):
                    cells[(side, k, key)] = QRectF(cx, top + head_h + r * row_h, cw, row_h)
                cx += cw
            ops.append(('line', x, top + head_h, x + w_total, top + head_h, C['border']))
            for r in range(1, len(rows)):
                y = top + head_h + r * row_h
                ops.append(('line', x, y, x + w_total, y, C['rule']))
            cx = x + label_w
            for cw in col_ws:
                ops.append(('line', cx, top, cx, top + h_total, C['rule']))
                cx += cw
            ops.append(('frame', card, C['border'], fs * 0.45))
            for r, (key, label) in enumerate(rows):
                base = top + head_h + r * row_h + row_h * 0.66
                put_text(label, x + pad, base, small, C['label'])
                cx = x + label_w
                for k, cw in enumerate(col_ws):
                    num, unit, typed = texts[k][key]
                    if num:
                        color = own or (C['typed'] if typed else C['value'])
                        w = m.width(num, fs, typed) + (m.width(unit, small) + fs * 0.2 if unit else 0)
                        end = put_text(num, cx + cw - pad - w, base, fs, color, bold=typed)
                        if unit:
                            put_text(unit, end + fs * 0.2, base, small, C['unit'])
                    elif key in inputs:                     # an empty cell that takes a value
                        put_text('·', cx + cw - pad - m.width('·', fs), base, fs, C['rule'])
                    cx += cw
            r_all = card.united(chip)
            bounds = r_all if bounds.isNull() else bounds.united(r_all)
            x += w_total + fs * 1.2
        if bounds.isNull():                                 # every compound of the grid was deleted
            bounds = QRectF(g.x * self.BL, g.y * self.BL, fs * 10, fs * 2)
            ops.append(('text', 'Stoichiometry: no reaction', g.x * self.BL + pad, g.y * self.BL + fs * 1.3, fs,
                        False, True, C['label']))
        return ops, bounds, cells

    def grid_rect(self, i):
        return self.grid_layout(i)[1]

    def paint_grid(self, p, i):
        ops, _bounds, _cells = self.grid_layout(i)
        width = max(0.5, self.s['line_width'] * 0.6)
        for op in ops:
            if op[0] == 'text':
                _t, txt, x, y, size, bold, italic, color = op
                if txt.strip():
                    p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor(color)))
                    p.drawPath(text_path(txt, self.s['font'], size, x, y, bold, italic))
            elif op[0] == 'line':
                p.setPen(QPen(QColor(op[5]), width)); p.drawLine(QPointF(op[1], op[2]), QPointF(op[3], op[4]))
            elif op[0] == 'fill':
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor(op[2])))
                p.drawRoundedRect(op[1], op[3], op[3])
            else:
                p.setPen(QPen(QColor(op[2]), width)); p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(op[1], op[3], op[3])

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
        self.paint_fills_and_highlights(p)
        for i, sh in doc.shapes.items():                    # pictures first: everything else goes over them
            if sh.kind == 'image':
                import pictures
                pictures.paint(self, p, i)
        for i, sh in doc.shapes.items():
            if sh.kind == 'image':
                continue
            if sh.kind in model.BIO_KINDS:
                import biodraw
                biodraw.paint_bio(self, p, i)
                continue
            if sh.kind in model.ELEMENT_KINDS or sh.style:
                import elements
                elements.paint(self, p, i)
                continue
            p.setPen(pen(self.color(sh))); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(self.shape_paths(i))
            rep = self.repeat_path(i)
            if rep is not None:
                p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(self.color(sh)))
                p.drawPath(rep)
        if doc.seqs:
            import biodraw
            for i in doc.seqs:
                biodraw.paint_seq(self, p, i)
        if doc.tlcs:
            import tlc
            for i in doc.tlcs:
                tlc.paint(self, p, i)
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
                p.drawPath(text_path(text, s['font'], size, x, y, bold=bool(s['label_bold'])))
        for i, ar in doc.arrows.items():
            color = self.color(ar)
            stroke, fill = self.arrow_paths(i)
            p.setPen(pen(color)); p.setBrush(Qt.BrushStyle.NoBrush); p.drawPath(stroke)
            p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color)); p.drawPath(fill)
        for i in doc.grids:
            self.paint_grid(p, i)
        p.setPen(Qt.PenStyle.NoPen)
        for i, t in doc.texts.items():
            p.setBrush(QBrush(self.color(t)))
            for txt, x, y, size in self.text_layout(i)[0]:
                if txt.strip():
                    p.drawPath(text_path(txt, s['caption_font'], size, x, y, t.bold, t.italic))
        self.paint_lone_pairs(p)
        if s['show_stereo'] or s['show_atom_numbers'] or any(at.map for at in doc.atoms.values()):
            self.paint_overlays(p)
        if errors:
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.setPen(QPen(QColor('#E0201B'), 0.8))
            for a, at in doc.atoms.items():
                if doc.valence_error(a):
                    r = 0.28 * self.BL
                    p.drawRect(QRectF(at.x * self.BL - r, at.y * self.BL - r, 2 * r, 2 * r))

    def paint_fills_and_highlights(self, p):
        """Behind the structures: rings filled with their colour, then highlight bands along bonds and discs on atoms
        (round-ended, so the bands of one structure join up)."""
        doc, BL = self.doc, self.BL
        p.setPen(Qt.PenStyle.NoPen)
        for f in doc.fills.values():
            if len(f.atoms) < 3 or any(a not in doc.atoms for a in f.atoms):
                continue
            p.setBrush(QBrush(QColor(f.color)))
            p.drawPolygon(QPolygonF([QPointF(doc.atoms[a].x * BL, doc.atoms[a].y * BL) for a in f.atoms]))
        width = 0.42 * BL
        for b in doc.bonds.values():
            if b.highlight:
                p.setPen(QPen(QColor(b.highlight), width, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
                pa, pb = doc.atoms[b.a], doc.atoms[b.b]
                p.drawLine(QPointF(pa.x * BL, pa.y * BL), QPointF(pb.x * BL, pb.y * BL))
        p.setPen(Qt.PenStyle.NoPen)
        for a in doc.atoms.values():
            if a.highlight:
                p.setBrush(QBrush(QColor(a.highlight)))
                p.drawEllipse(QPointF(a.x * BL, a.y * BL), width / 2 * 1.15, width / 2 * 1.15)

    def lone_pair_dots(self, a):
        """[(x, y)] in pt of the dots of atom a's lone pairs, two per pair, in the widest gaps between its bonds
        and its H."""
        doc, BL = self.doc, self.BL
        at = doc.atoms[a]
        n = doc.lone_pair_count(a) if at.lone_pairs < 0 else min(at.lone_pairs, 4)
        if n <= 0:
            return []
        cx, cy = at.x * BL, at.y * BL
        taken = [math.atan2(doc.atoms[o].y - at.y, doc.atoms[o].x - at.x) for o, _b in doc.neighbors(a)]
        lab = self.labels.get(a)
        if lab and lab['clip']:                             # just round the ink of the letters themselves (Cl: the l)
            sym, sx, sy, size = lab['pieces'][0]
            ink = text_path(sym, self.s['font'], size, sx, sy, bool(self.s['label_bold'])).boundingRect()
        else:
            ink = QRectF(cx - 0.12 * BL, cy - 0.12 * BL, 0.24 * BL, 0.24 * BL)
        gap = max(self.s['line_width'] * 1.6, 1.4)
        dot = max(self.s['line_width'] * 0.75, 0.55)
        letters = []                                        # what the dots must keep off: other letters, own H ...
        for o, ol in self.labels.items():
            if ol and abs(doc.atoms[o].x - at.x) < 2.5 and abs(doc.atoms[o].y - at.y) < 2.5:
                for k, (txt, x, y, size) in enumerate(ol['pieces']):
                    if txt.strip() and not (o == a and k == 0):
                        letters.append(QRectF(x, y - self.m.cap(size), self.m.width(txt, size), self.m.cap(size)))

        def spot(ang):
            ux, uy = math.cos(ang), math.sin(ang)
            kx = ((ink.right() if ux > 0 else ink.left()) - cx) / ux if abs(ux) > 1e-6 else 1e9
            ky = ((ink.bottom() if uy > 0 else ink.top()) - cy) / uy if abs(uy) > 1e-6 else 1e9
            k = max(0.0, min(kx, ky))
            px, py = cx + ux * (k + 1.8), cy + uy * (k + 1.8)
            return [(px - uy * gap, py + ux * gap), (px + uy * gap, py - ux * gap)]

        def clear(pts):
            return all(self.is_free(QRectF(x - dot - 0.3, y - dot - 0.3, 2 * dot + 0.6, 2 * dot + 0.6), letters)
                       for x, y in pts)
        if lab:
            for txt, x, y, size in lab['pieces'][1:]:          # its H, charge and dots take room too
                if txt.strip():
                    tx, ty = x + self.m.width(txt, size) / 2 - cx, y - self.m.cap(size) / 2 - cy
                    if abs(tx) > 0.5 or abs(ty) > 0.5:
                        taken.append(math.atan2(ty, tx))
        chosen = []
        for _k in range(n):                                 # the free direction farthest from bonds, H and others:
            best, best_d = 0.0, -1.0                        # up / down / left / right first, then the diagonals
            for step, ang in enumerate((-math.pi / 2, math.pi / 2, 0.0, math.pi, -math.pi / 4, -3 * math.pi / 4,
                                        math.pi / 4, 3 * math.pi / 4)):
                d = min((abs((ang - t + math.pi) % (2 * math.pi) - math.pi) for t in taken + chosen), default=math.pi)
                d += 0.15 if step < 4 else 0.0
                if not clear(spot(ang)):                    # on a bond or on letters: only if nothing else is left
                    d -= 10
                if d > best_d + 1e-9:
                    best, best_d = ang, d
            chosen.append(best)
        out = []
        for ang in chosen:
            out += spot(ang)
        return out

    def paint_lone_pairs(self, p):
        r = max(self.s['line_width'] * 0.75, 0.55)
        p.setPen(Qt.PenStyle.NoPen)
        for a, at in self.doc.atoms.items():
            if at.lone_pairs:
                p.setBrush(QBrush(self.color(at, a)))
                for x, y in self.lone_pair_dots(a):
                    p.drawEllipse(QPointF(x, y), r, r)

    def paint_overlays(self, p):
        """(R) / (S) and (E) / (Z) labels, atom numbers and atom-map numbers, small and blue."""
        s = self.s
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(QColor(OVERLAY)))
        for text, x, y, size in self.overlay_layout():
            p.drawPath(text_path(text, s['font'], size, x, y))

    def overlay_layout(self):
        """[(text, x, baseline, size)] of the overlays (cached): beside their atoms and bonds."""
        if getattr(self, '_overlays', None) is None:
            self._overlays = self._overlay_layout()
        return self._overlays

    def _overlay_layout(self):
        s, doc, BL = self.s, self.doc, self.BL
        size = s['font_size'] * 0.7
        out = []
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
        for a, at in doc.atoms.items():                     # reaction atom map: :3
            if at.map:
                notes.setdefault(('atom', a), []).append(f':{at.map}')
        if not notes:
            return out
        taken = self.obstacles()                            # letters, dots, captions: rects; lines: bond_regions()
        cap = self.m.cap(size)
        for (kind, i), parts in sorted(notes.items(), key=lambda kv: kv[0][0] == 'bond'):
            text = ' '.join(parts)
            w = self.m.width(text, size)
            ink = text_path(text, s['font'], size, 0.0, 0.0).boundingRect()   # brackets reach below the line

            def box_at(x, y):
                return ink.translated(x - w / 2, y + cap / 2)
            spots = []
            if kind == 'atom':                              # around the atom, nearest and most open first
                at = doc.atoms[i]
                fx, fy = doc.free_direction(i)
                lab = self.labels.get(i)
                clip = max(lab['clip']) if lab and lab['clip'] else 0.0
                base = math.atan2(fy, fx)
                for reach in (0.42, 0.6, 0.8, 1.0, 1.25):
                    for k in range(24):
                        ang = base + math.radians(15 * ((k + 1) // 2) * (1 if k % 2 else -1))
                        r = reach * BL + clip + max(w, cap) * 0.25 * abs(math.cos(ang))
                        spots.append((at.x * BL + math.cos(ang) * r, at.y * BL + math.sin(ang) * r))
            else:                                           # beside the bond's middle, the open side first
                bd = doc.bonds[i]
                pa, pb = doc.atoms[bd.a], doc.atoms[bd.b]
                L = math.hypot(pb.x - pa.x, pb.y - pa.y) or 1.0
                ux, uy = (pb.x - pa.x) / L, (pb.y - pa.y) / L
                side = -(doc.double_side(i) or 1)
                mx, my = (pa.x + pb.x) / 2 * BL, (pa.y + pb.y) / 2 * BL
                for dist in (0.42, 0.62, 0.85):
                    for sgn in (side, -side):
                        for t in (0.0, 0.22, -0.22):
                            spots.append((mx - uy * dist * BL * sgn + ux * t * BL, my + ux * dist * BL * sgn + uy * t * BL))
            place, least, found = spots[0], None, False
            for gx, gy in ((1.0, 0.8), (0.4, 0.4)):         # a visible gap if there is room, else just clear
                for x, y in spots:
                    box = box_at(x, y).adjusted(-gx, -gy, gx, gy)
                    if self.is_free(box, taken):
                        place, found = (x, y), True
                        break
                    hits = self.hits(box, taken)            # nowhere clear: the least crowded spot
                    if least is None or hits < least:
                        place, least = (x, y), hits
                if found:
                    break
            x, y = place
            taken.append(box_at(x, y))
            out.append((text, x - w / 2, y + cap / 2, size))
        return out

    def bond_regions(self):
        """[(('bond' | 'arrow', id), QPainterPath of its drawn lines widened by their width, its bounding rect)]."""
        if getattr(self, '_regions', None) is None:
            from PyQt6.QtGui import QPainterPathStroker
            st = QPainterPathStroker()
            st.setWidth(max(self.s['line_width'], 0.5))
            out = []
            for kind, table, paths in (('bond', self.doc.bonds, self.bond_paths), ('arrow', self.doc.arrows, self.arrow_paths)):
                for i in table:
                    stroke, fill = paths(i)
                    region = st.createStroke(stroke)
                    region.addPath(fill)
                    out.append(((kind, i), region, region.boundingRect()))
            self._regions = out
        return self._regions

    def obstacles(self):
        """Rects (pt) of what is written or dotted: label letters, lone pairs, captions, repeat counts."""
        s = self.s
        out = []
        for a, lab in self.labels.items():
            if lab:
                for text, x, y, size in lab['pieces']:
                    if text.strip():
                        out.append(text_path(text, s['font'], size, x, y, bool(s['label_bold'])).boundingRect())
        r = max(s['line_width'] * 0.75, 0.55)
        for a, at in self.doc.atoms.items():
            if at.lone_pairs:
                for x, y in self.lone_pair_dots(a):
                    out.append(QRectF(x - r, y - r, 2 * r, 2 * r))
        for i in self.doc.texts:
            out.append(self.text_rect(i))
        for i in self.doc.shapes:
            p = self.repeat_path(i)
            if p is not None:
                out.append(p.boundingRect())
        return out

    def hits(self, box, rects):
        """How many of rects, bonds and arrows box touches."""
        n = sum(1 for r in rects if box.intersects(r))
        return n + sum(1 for _k, region, brect in self.bond_regions() if brect.intersects(box) and region.intersects(box))

    def is_free(self, box, rects, lines=True):
        """True when box touches none of rects nor (lines) any bond or arrow."""
        if any(box.intersects(r) for r in rects):
            return False
        if lines:
            for _key, region, brect in self.bond_regions():
                if brect.intersects(box) and region.intersects(box):
                    return False
        return True

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
            rep = self.repeat_path(i)
            if rep is not None:
                add(rep.boundingRect().adjusted(-1, -1, 1, 1))
        for i in self.doc.grids:
            add(self.grid_rect(i).adjusted(-1, -1, 1, 1))
        for i in self.doc.seqs:
            add(self.item_rect('seq', i).adjusted(-1, -1, 1, 1))
        for i in self.doc.tlcs:
            add(self.item_rect('tlc', i).adjusted(-1, -1, 1, 1))
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
        if kind == 'grid':
            return self.grid_rect(i)
        if kind == 'seq':
            import biodraw
            return biodraw.seq_layout(self, i)['rect']
        if kind == 'tlc':
            import tlc
            return tlc.layout(self, i)['rect']
        r = self.shape_paths(i).boundingRect()
        rep = self.repeat_path(i)
        return r.united(rep.boundingRect()) if rep is not None else r


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
