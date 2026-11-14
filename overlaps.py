"""Finds where things drawn sit on top of each other and make the drawing hard to read: letters of different labels,
numbers (atom numbers, R/S, map numbers, charges, H counts) on bonds or on other letters, lone-pair dots on bonds,
captions over structures, atoms laid on top of each other. Used by the tests to keep the drawing readable.

    issues = overlaps.find(doc)        # [Issue(kind, a, b, text)]
"""
from dataclasses import dataclass

from PyQt6.QtCore import QRectF

import render

SHRINK = 0.25                                               # pt: touching is fine, overlapping is not


@dataclass
class Issue:
    kind: str                                               # 'glyph-glyph', 'glyph-bond', 'dot-bond', 'dot-glyph',
    a: tuple                                                # 'caption', 'clash' ...
    b: tuple
    text: str


def _ink(text, family, size, x, y, bold=False, italic=False):
    return render.text_path(text, family, size, x, y, bold, italic).boundingRect()


def glyphs(rend):
    """[(owner, piece text, QRectF ink box)] of everything written: label pieces, overlays, captions, repeat counts."""
    s, doc = rend.s, rend.doc
    out = []
    for a, lab in rend.labels.items():
        if not lab:
            continue
        for k, (text, x, y, size) in enumerate(lab['pieces']):
            if text.strip():
                out.append((('atom', a, k), text, _ink(text, s['font'], size, x, y, bool(s['label_bold']))))
    for k, (text, x, y, size) in enumerate(rend.overlay_layout()):
        out.append((('overlay', k), text, _ink(text, s['font'], size, x, y)))
    for i, t in doc.texts.items():
        for k, (txt, x, y, size) in enumerate(rend.text_layout(i)[0]):
            if txt.strip():
                out.append((('text', i, k), txt, _ink(txt, s['caption_font'], size, x, y, t.bold, t.italic)))
    for i in doc.shapes:
        p = rend.repeat_path(i)
        if p is not None:
            out.append((('repeat', i), 'n', p.boundingRect()))
    return out


def dots(rend):
    """[(owner, QRectF)] of each lone pair (both dots)."""
    out = []
    for a, at in rend.doc.atoms.items():
        if at.lone_pairs:
            pts = rend.lone_pair_dots(a)
            r = max(rend.s['line_width'] * 0.75, 0.55)
            for k, (x, y) in enumerate(pts):              # each dot on its own
                out.append((('lp', a, k // 2), QRectF(x - r, y - r, 2 * r, 2 * r)))
    return out


def bond_regions(rend):
    """[(('bond', id), QPainterPath of the drawn lines, widened by their width)] and the same for arrows."""
    return [(key, region) for key, region, _r in rend.bond_regions()]


def _small(r, d=SHRINK):
    return r.adjusted(d, d, -d, -d)


def find(doc, rend=None):
    rend = rend or render.Renderer(doc)
    out = []
    gl = glyphs(rend)
    regions = bond_regions(rend)
    owner_atom = lambda o: o[1] if o[0] == 'atom' else None

    # letters / numbers on each other (pieces of one label may touch: H and its count, symbol and charge)
    for k in range(len(gl)):
        o1, t1, r1 = gl[k]
        for j in range(k + 1, len(gl)):
            o2, t2, r2 = gl[j]
            if o1[:2] == o2[:2] and o1[0] in ('atom', 'text'):
                if not _small(r1, 0.6).intersects(_small(r2, 0.6)):
                    continue
            elif not _small(r1).intersects(_small(r2)):
                continue
            out.append(Issue('glyph-glyph', o1, o2, f'"{t1}" on "{t2}"'))
    # letters / numbers on lines (a label's own bonds are cut back at it)
    for o, t, r in gl:
        box = _small(r)
        for bo, region in regions:
            if not region.intersects(box):
                continue
            if bo[0] == 'bond' and o[0] == 'atom':
                bd = doc.bonds[bo[1]]
                if owner_atom(o) in (bd.a, bd.b) and o[2] == 0:   # the symbol itself: the bond ends at its edge
                    continue
            out.append(Issue('glyph-bond', o, bo, f'"{t}" on {bo[0]} {bo[1]}'))
    # lone pairs
    for o, r in dots(rend):
        for bo, region in regions:
            if region.intersects(_small(r, 0.1)):
                out.append(Issue('dot-bond', o, bo, f'lone pair of atom {o[1]} on {bo[0]} {bo[1]}'))
        for go, t, gr in gl:
            if _small(gr, 0.1).intersects(_small(r, 0.1)):
                out.append(Issue('dot-glyph', o, go, f'lone pair of atom {o[1]} on "{t}"'))
    # atoms on atoms (a layout clash: two atoms where one should be)
    ids = list(doc.atoms)
    for k in range(len(ids)):
        a = doc.atoms[ids[k]]
        for j in range(k + 1, len(ids)):
            b = doc.atoms[ids[j]]
            if (a.x - b.x) ** 2 + (a.y - b.y) ** 2 < 0.3 ** 2 and doc.bond_between(ids[k], ids[j]) is None:
                out.append(Issue('clash', ('atom', ids[k]), ('atom', ids[j]), 'two atoms on one place'))
    return out


def summary(issues):
    from collections import Counter
    return Counter(i.kind for i in issues)

