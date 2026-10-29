"""
cdxml.py - ChemDraw's XML format (.cdxml), both ways, so drawings move between ChemDraw and this editor.

    doc, notes = cdxml.read(path)       # notes: what could not be brought over (stoichiometry grids ...)
    cdxml.write(doc, path)

Read: structures (elements, charges, nicknames such as OMe / Boc, bond orders, wedges, hashes, bold, dashed, wavy,
double-bond position), captions (size, bold, italic, formula style, justification, colour), arrows (plain, dashed,
equilibrium, resonance, retrosynthetic, curved and fish-hook), rectangles, ovals, lines, brackets; the document style
(bond length, line widths, font size). Written: the same. Coordinates in CDXML are points with y down, so a model
unit is one bond length (the file's BondLength).
"""
import math
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape, quoteattr

import model

ASCENT = 0.905                                              # Arial ascent / size: CDXML text sits on its baseline
FACE_BOLD, FACE_ITALIC, FACE_FORMULA = 1, 2, 96
NICK_TYPES = ('Nickname', 'Fragment', 'GenericNickname', 'Unspecified', 'MultiAttachment', 'Unknown')


def _reversed_nicknames():
    out = {}
    for k in model.NICKNAMES:
        toks, _anchor = model.label_order(k, True)
        out[''.join(toks)] = k
    return out


def _floats(s):
    return [float(v) for v in (s or '').split()]


def _text_of(t):
    """The text of a <t>: its <s> runs joined (line breaks kept)."""
    if t is None:
        return ''
    return ''.join(s.text or '' for s in t.iter('s')).replace('\r\n', '\n').replace('\r', '\n')


# ------------------------------------------------------------------------------------------------ read
def read(path):
    root = ET.parse(path).getroot()
    doc = model.Doc()
    notes = []
    BL = float(root.get('BondLength', 14.4)) or 14.4
    style = {'bond_length': BL}
    for key, attr in (('line_width', 'LineWidth'), ('bold_width', 'BoldWidth'), ('hash_spacing', 'HashSpacing'),
                      ('margin', 'MarginWidth'), ('font_size', 'LabelSize')):
        if root.get(attr):
            style[key] = float(root.get(attr))
    if root.get('BondSpacing'):
        style['bond_spacing'] = float(root.get('BondSpacing')) / 100
    if root.get('ShowAtomStereo') == 'yes':
        style['show_stereo'] = True
    if root.get('ShowAtomNumber') == 'yes':
        style['show_atom_numbers'] = True
    fonts = {f.get('id'): f.get('name') for f in root.iter('font')}
    if fonts.get(root.get('LabelFont')):
        style['font'] = fonts[root.get('LabelFont')]
    doc.style = style
    table = root.find('colortable')
    colors = []
    if table is not None:
        for c in table.findall('color'):
            colors.append('#%02X%02X%02X' % tuple(round(float(c.get(k, 0)) * 255) for k in 'rgb'))

    def color(el):
        c = el.get('color')
        if c is None:
            return ''
        k = int(c)
        if k == 0:
            return ''
        if k == 1:
            return '#FFFFFF'
        hexc = colors[k - 2] if 0 <= k - 2 < len(colors) else ''
        return '' if hexc in ('#000000', '') else hexc

    rev = _reversed_nicknames()
    nodes = {}                                              # CDXML node id -> atom id
    frag_atoms = {}                                         # CDXML fragment id -> atom ids
    grids = []
    unknown = set()

    def xy(p):
        x, y = _floats(p)[:2]
        return x / BL, y / BL

    def walk(parent):
        for el in parent:
            tag = el.tag
            if tag in ('fragment', 'group'):
                if tag == 'fragment':
                    read_fragment(el)
                else:
                    walk(el)
            elif tag == 't':
                read_text(el)
            elif tag == 'arrow':
                read_arrow(el)
            elif tag == 'graphic':
                read_graphic(el)
            elif tag == 'stoichiometrygrid':
                grids.append(el)
            elif tag in ('page',):
                walk(el)

    def read_fragment(frag):
        for n in frag.findall('n'):
            if n.get('NodeType') == 'ExternalConnectionPoint' or not n.get('p'):
                continue
            x, y = xy(n.get('p'))
            nt = n.get('NodeType', 'Element')
            if nt in NICK_TYPES:
                text = _text_of(n.find('t')).strip()
                label = text if text in model.NICKNAMES else rev.get(text, '')
                if not label:
                    inner = [m for m in n.iter('n') if m is not n and m.get('NodeType') != 'ExternalConnectionPoint']
                    unknown.add(text or '?')
                    label = text
                a = doc.add_atom(x, y)
                doc.atoms[a].label = label
            else:
                z = int(n.get('Element', 6))
                try:
                    from rdkit import Chem
                    el = Chem.GetPeriodicTable().GetElementSymbol(z)
                except Exception:
                    el = 'C'
                a = doc.add_atom(x, y, el, int(n.get('Charge', 0)))
            doc.atoms[a].color = color(n)
            nodes[n.get('id')] = a
            frag_atoms.setdefault(frag.get('id'), []).append(a)
        for b in frag.findall('b'):
            a1, a2 = nodes.get(b.get('B')), nodes.get(b.get('E'))
            if a1 is None or a2 is None:
                continue
            order = {'2': 2, '3': 3, '1.5': 2, '2.5': 3}.get(b.get('Order', '1'), 1)
            disp = b.get('Display', '')
            i = doc.add_bond(a1, a2, order)
            if i is None:
                continue
            bd = doc.bonds[i]
            if disp in ('WedgeBegin', 'WedgeEnd', 'WedgedHashBegin', 'WedgedHashEnd', 'HollowWedgeBegin', 'HollowWedgeEnd'):
                bd.stereo = 'hash' if 'Hash' in disp else 'wedge'
                bd.order = 1
                if disp.endswith('End'):
                    bd.a, bd.b = bd.b, bd.a
            else:
                bd.style = {'Bold': 'bold', 'Dash': 'dashed', 'Wavy': 'wavy', 'Hash': 'hashed'}.get(disp, '')
                if order == 2 and b.get('Display2') == 'Dash':
                    bd.style = 'dashed'
            pos = b.get('DoublePosition')
            if order == 2 and pos:
                bd.side = {'Center': 0, 'Left': -1, 'Right': 1}.get(pos)
            bd.color = color(b)

    def read_text(t):
        if not t.get('p'):
            return
        text = _text_of(t)
        if not text.strip():
            return
        s0 = t.find('s')
        size = float(s0.get('size', style.get('font_size', 10))) if s0 is not None else style.get('font_size', 10)
        face = int(s0.get('face', 0)) if s0 is not None else 0
        x, y = _floats(t.get('p'))[:2]
        just = (t.get('Justification') or t.get('CaptionJustification') or 'Left').lower()
        align = 'center' if just.startswith('center') else 'right' if just.startswith('right') else 'left'
        i = doc.add_text(x / BL, (y - ASCENT * size) / BL, text, align=align,
                         size=size if abs(size - style.get('font_size', 10)) > 0.01 else 0.0,
                         bold=bool(face & FACE_BOLD), italic=bool(face & FACE_ITALIC),
                         formula=(face & FACE_FORMULA) == FACE_FORMULA)
        doc.texts[i].color = color(s0 if s0 is not None and s0.get('color') else t)

    def read_arrow(ar):
        head, tail = _floats(ar.get('Head3D'))[:2], _floats(ar.get('Tail3D'))[:2]
        if len(head) < 2 or len(tail) < 2:
            return
        hh, ht = ar.get('ArrowheadHead', 'Full'), ar.get('ArrowheadTail', '')
        kind, bend = 'forward', 0.0
        ang = float(ar.get('AngularSize', 0) or 0)
        if ang:
            kind = 'curved_half' if hh.startswith('Half') else 'curved'
            bend = math.tan(math.radians(abs(ang)) / 4) / 2
            c = _floats(ar.get('Center3D'))[:2]
            if len(c) == 2:                                  # the bulge is on the side away from the arc's centre
                dx, dy = head[0] - tail[0], head[1] - tail[1]
                mx, my = (head[0] + tail[0]) / 2, (head[1] + tail[1]) / 2
                if (c[0] - mx) * dy - (c[1] - my) * dx > 0:
                    bend = -bend
        elif ar.get('ArrowheadType') == 'Angle' and ar.get('ArrowShaftSpacing'):
            kind = 'retro'
        elif ht.startswith('Half') or (hh.startswith('Half') and ar.get('ArrowShaftSpacing')):
            kind = 'equilibrium'
        elif ht == 'Full':
            kind = 'resonance'
        elif ar.get('LineType') == 'Dashed':
            kind = 'dashed'
        i = doc.add_arrow(tail[0] / BL, tail[1] / BL, head[0] / BL, head[1] / BL, kind, bend)
        doc.arrows[i].color = color(ar)

    def read_graphic(g):
        if g.get('SupersededBy'):
            return
        gt = g.get('GraphicType', '')
        bb = _floats(g.get('BoundingBox'))
        if len(bb) < 4:
            return
        x1, y1, x2, y2 = (v / BL for v in bb[:4])
        if gt == 'Line':
            at = g.get('ArrowType', 'NoHead')
            if at == 'NoHead':
                i = doc.add_shape(x1, y1, x2, y2, 'line')
                doc.shapes[i].color = color(g)
                return
            kind = {'Equilibrium': 'equilibrium', 'Resonance': 'resonance', 'RetroSynthetic': 'retro',
                    'HalfHead': 'curved_half'}.get(at, 'forward')
            if kind == 'curved_half':
                kind = 'forward'
            i = doc.add_arrow(x2, y2, x1, y1, kind)         # a line graphic lists the head first
            doc.arrows[i].color = color(g)
        elif gt in ('Rectangle', 'Oval', 'Bracket'):
            if gt == 'Oval' and g.get('Center3D'):
                c = _floats(g.get('Center3D')); ma = _floats(g.get('MajorAxisEnd3D')); mi = _floats(g.get('MinorAxisEnd3D'))
                if len(c) >= 2 and len(ma) >= 2 and len(mi) >= 2:
                    rx = math.hypot(ma[0] - c[0], ma[1] - c[1]); ry = math.hypot(mi[0] - c[0], mi[1] - c[1])
                    x1, y1, x2, y2 = (c[0] - rx) / BL, (c[1] - ry) / BL, (c[0] + rx) / BL, (c[1] + ry) / BL
            if gt == 'Rectangle':
                kind = 'roundrect' if 'Round' in (g.get('RectangleType') or '') else 'rect'
            elif gt == 'Oval':
                kind = 'oval'
            else:
                bt = g.get('BracketType', 'SquarePair')
                kind = 'parens' if 'Round' in bt else 'braces' if 'Curly' in bt else 'brackets'
            i = doc.add_shape(x1, y1, x2, y2, kind)
            doc.shapes[i].color = color(g)

    walk(root)
    for g in grids:
        read_grid(g, doc, frag_atoms, BL)
    if unknown:
        notes.append('these labels, kept as text but counted as carbon: ' + ', '.join(sorted(unknown)))
    return doc, notes


# A grid datum's SGPropertyType -> our row key; values as ChemDraw keeps them (fractions for percentages)
SG_REACTANT = {4: 'limiting', 6: 'equiv', 7: 'mass', 8: 'pct', 9: 'volume', 10: 'molarity', 12: 'density', 13: 'moles'}
SG_PRODUCT = {6: 'equiv', 22: 'completion', 17: 'measured', 18: 'purity'}
SG_DEFAULTS = {'equiv': 1.0, 'pct': 1.0, 'completion': 1.0, 'purity': 1.0}


def read_grid(g, doc, frag_atoms, BL):
    """A ChemDraw stoichiometry grid: its columns (by the fragment each one describes) and the values typed in."""
    p = _floats(g.get('p'))
    grid_id = doc._id()
    grid = model.Grid(p[0] / BL if p else 0.0, (p[1] / BL if p else 0.0) + 1.2)
    for comp in g.findall('sgcomponent'):
        if comp.get('ComponentIsHeader') == 'yes':
            continue
        atoms = frag_atoms.get(comp.get('ComponentReferenceID'))
        if not atoms:
            continue
        reactant = comp.get('ComponentIsReactant') == 'yes'
        keys = SG_REACTANT if reactant else SG_PRODUCT
        col = {'atoms': list(atoms)}
        for d in comp.findall('sgdatum'):
            key = keys.get(int(d.get('SGPropertyType', 0)))
            if key is None:
                continue
            try:
                v = float(d.get('SGDataValue', ''))
            except ValueError:
                continue
            if key == 'limiting':
                col['limiting'] = v == 1
                continue
            if d.get('IsEdited') != 'yes' and (v == SG_DEFAULTS.get(key, 0.0)):
                continue                                    # an untouched default
            if key in ('pct', 'completion', 'purity'):
                v *= 100
            col[key] = v
        (grid.reactants if reactant else grid.products).append(col)
    if grid.reactants or grid.products:
        doc.grids[grid_id] = grid


# ------------------------------------------------------------------------------------------------ write
def write(doc, path, program='Structure editor'):
    import render
    from rdkit import Chem
    s = render.style_of(doc)
    r = render.Renderer(doc)
    BL = s['bond_length']
    ids = iter(range(1000, 10 ** 9))
    colors = []

    def cidx(c):
        """CDXML colour index of '#RRGGBB' ('' = the default, black)."""
        if not c:
            return None
        c = c.upper()
        if c not in colors:
            colors.append(c)
        return 4 + colors.index(c)                          # 2 = white, 3 = black, then ours

    def cattr(c):
        k = cidx(c)
        return f' color="{k}"' if k is not None else ''

    def f(v):
        return f'{v:.2f}'

    out = []
    pt = lambda x, y: f'{f(x * BL)} {f(y * BL)}'
    PT = Chem.GetPeriodicTable()
    for frag in doc.fragments():
        fid = next(ids)
        node_id = {}
        body = []
        for a in frag:
            node_id[a] = next(ids)
        for a in frag:
            at = doc.atoms[a]
            lab = r.labels.get(a)
            label_xml = ''
            if lab and lab['clip']:                          # a label: ChemDraw shows the text it is given
                text = ''.join(p[0] for p in lab['pieces']).replace('−', '-')
                x0, base = lab['pieces'][0][1], lab['pieces'][0][2]
                just = 'Right' if (text and not text.startswith(at.label[:1] if at.label else at.el)) else 'Left'
                label_xml = (f'<t p="{f(x0)} {f(base)}" LabelJustification="{just}"><s font="3" size="{f(s["font_size"])}" '
                             f'face="96"{cattr(at.color)}>{escape(text)}</s></t>')
            if at.label and at.label in model.NICKNAMES:
                body.append(_nickname_node(doc, a, node_id[a], ids, pt, label_xml, cattr(at.color)))
                continue
            attrs = f'id="{node_id[a]}" p="{pt(at.x, at.y)}"'
            if at.el != 'C':
                attrs += f' Element="{PT.GetAtomicNumber(at.el)}"'
            if at.charge:
                attrs += f' Charge="{at.charge}"'
            if label_xml or at.el != 'C':
                attrs += f' NumHydrogens="{doc.implicit_h(a)}" NeedsClean="yes"'
            attrs += cattr(at.color)
            body.append(f'<n {attrs}>{label_xml}</n>' if label_xml else f'<n {attrs}/>')
        for i, bd in doc.bonds.items():
            if bd.a not in node_id or bd.b not in node_id:
                continue
            attrs = f'id="{next(ids)}" B="{node_id[bd.a]}" E="{node_id[bd.b]}"'
            if bd.order != 1:
                attrs += f' Order="{bd.order}"'
            disp = {'wedge': 'WedgeBegin', 'hash': 'WedgedHashBegin'}.get(bd.stereo) if bd.order == 1 else None
            if disp is None and bd.order == 1:
                disp = {'bold': 'Bold', 'dashed': 'Dash', 'wavy': 'Wavy', 'hashed': 'Hash'}.get(bd.style)
            if disp:
                attrs += f' Display="{disp}"'
            if bd.order == 2:
                if bd.style == 'dashed':
                    attrs += ' Display2="Dash"'
                if bd.side is not None:
                    attrs += f' DoublePosition="{ {0: "Center", -1: "Left", 1: "Right"}[bd.side] }"'
            attrs += cattr(bd.color)
            body.append(f'<b {attrs}/>')
        out.append(f'<fragment id="{fid}">' + ''.join(body) + '</fragment>')
    for i, t in doc.texts.items():
        size = t.size or s['font_size']
        face = (FACE_BOLD if t.bold else 0) | (FACE_ITALIC if t.italic else 0) | (FACE_FORMULA if t.formula else 0)
        just = {'left': 'Left', 'center': 'Center', 'right': 'Right'}[t.align]
        out.append(f'<t id="{next(ids)}" p="{f(t.x * BL)} {f(t.y * BL + ASCENT * size)}" Justification="{just}" '
                   f'CaptionJustification="{just}" InterpretChemically="no"><s font="3" size="{f(size)}" face="{face}"'
                   f'{cattr(t.color)}>{escape(t.text)}</s></t>')
    for i, ar in doc.arrows.items():
        head, tail = pt(ar.x2, ar.y2) + ' 0', pt(ar.x1, ar.y1) + ' 0'
        attrs = (f'id="{next(ids)}" FillType="None" ArrowheadHead="Full" ArrowheadType="Solid" HeadSize="1000" '
                 f'ArrowheadCenterSize="875" ArrowheadWidth="250" Head3D="{head}" Tail3D="{tail}"')
        if ar.kind == 'equilibrium':
            attrs = attrs.replace('ArrowheadHead="Full"', 'ArrowheadHead="HalfLeft" ArrowheadTail="HalfLeft" ArrowShaftSpacing="300"')
        elif ar.kind == 'resonance':
            attrs = attrs.replace('ArrowheadHead="Full"', 'ArrowheadHead="Full" ArrowheadTail="Full"')
        elif ar.kind == 'retro':
            attrs = attrs.replace('ArrowheadType="Solid"', 'ArrowheadType="Angle" ArrowShaftSpacing="300"')
        elif ar.kind == 'dashed':
            attrs += ' LineType="Dashed"'
        elif ar.kind in ('curved', 'curved_half') and ar.bend:
            if ar.kind == 'curved_half':
                attrs = attrs.replace('ArrowheadHead="Full"', 'ArrowheadHead="HalfLeft"')
            ang = math.degrees(4 * math.atan(2 * abs(ar.bend)))
            L = math.hypot(ar.x2 - ar.x1, ar.y2 - ar.y1) * BL
            R = L / (2 * math.sin(math.radians(ang) / 2))
            mx, my = (ar.x1 + ar.x2) / 2 * BL, (ar.y1 + ar.y2) / 2 * BL
            cx, cy = model.curve_control(ar)
            ux, uy = mx - cx * BL, my - cy * BL              # from the bulge back through the middle
            ul = math.hypot(ux, uy) or 1.0
            d = R * math.cos(math.radians(ang) / 2)
            ctr = (mx + ux / ul * d, my + uy / ul * d)
            attrs += (f' AngularSize="{f(ang if ar.bend > 0 else -ang)}" Center3D="{f(ctr[0])} {f(ctr[1])} 0"'
                      f' MajorAxisEnd3D="{f(ctr[0] + R)} {f(ctr[1])} 0" MinorAxisEnd3D="{f(ctr[0])} {f(ctr[1] + R)} 0"')
        out.append(f'<arrow {attrs}{cattr(ar.color)}/>')
    for i, sh in doc.shapes.items():
        bb = f'{f(sh.x1 * BL)} {f(sh.y1 * BL)} {f(sh.x2 * BL)} {f(sh.y2 * BL)}'
        if sh.kind == 'line':
            g = f'GraphicType="Line" ArrowType="NoHead"'
        elif sh.kind in ('rect', 'roundrect'):
            g = 'GraphicType="Rectangle"' + (' RectangleType="RoundEdge"' if sh.kind == 'roundrect' else '')
        elif sh.kind == 'oval':
            g = 'GraphicType="Oval" OvalType="Circle"' if abs(abs(sh.x2 - sh.x1) - abs(sh.y2 - sh.y1)) < 1e-6 else 'GraphicType="Oval"'
        else:
            g = f'GraphicType="Bracket" BracketType="{ {"brackets": "SquarePair", "parens": "RoundPair", "braces": "CurlyPair"}[sh.kind] }"'
        out.append(f'<graphic id="{next(ids)}" BoundingBox="{bb}" {g}{cattr(sh.color)}/>')
    b = r.bounds()
    W, H = render.PAGES.get(s['page'], render.PAGES['A4'])
    bbox = f'{f(b.left())} {f(b.top())} {f(b.right())} {f(b.bottom())}' if not b.isNull() else '0 0 0 0'
    table = '<color r="1" g="1" b="1"/><color r="0" g="0" b="0"/>' + ''.join(
        f'<color r="{int(c[1:3], 16) / 255:.3f}" g="{int(c[3:5], 16) / 255:.3f}" b="{int(c[5:7], 16) / 255:.3f}"/>'
        for c in colors)
    head = (f'<?xml version="1.0" encoding="UTF-8" ?>\n'
            f'<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n'
            f'<CDXML CreationProgram={quoteattr(program)} BoundingBox="{bbox}" FractionalWidths="yes" '
            f'InterpretChemically="yes" ShowAtomStereo="{"yes" if s["show_stereo"] else "no"}" '
            f'ShowBondStereo="{"yes" if s["show_stereo"] else "no"}" '
            f'ShowAtomNumber="{"yes" if s["show_atom_numbers"] else "no"}" LabelFont="3" LabelSize="{f(s["font_size"])}" '
            f'LabelFace="96" CaptionFont="3" CaptionSize="{f(s["font_size"])}" HashSpacing="{f(s["hash_spacing"])}" '
            f'MarginWidth="{f(s["margin"])}" LineWidth="{f(s["line_width"])}" BoldWidth="{f(s["bold_width"])}" '
            f'BondLength="{f(BL)}" BondSpacing="{round(s["bond_spacing"] * 100)}" ChainAngle="120" color="0" bgcolor="1">'
            f'<colortable>{table}</colortable><fonttable><font id="3" charset="iso-8859-1" name={quoteattr(s["font"])}/>'
            f'</fonttable><page id="{next(ids)}" BoundingBox="0 0 {f(W)} {f(H)}" HeaderPosition="36" FooterPosition="36" '
            f'PrintTrimMarks="yes" HeightPages="1" WidthPages="1">')
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(head + '\n'.join(out) + '</page></CDXML>\n')


def _nickname_node(doc, a, nid, ids, pt, label_xml, cattr):
    """A ChemDraw nickname: the outer node with the label, holding the group's atoms and a connection point for each
    bond to the rest of the structure."""
    from rdkit import Chem
    PT = Chem.GetPeriodicTable()
    at = doc.atoms[a]
    tmp = model.Doc()
    x = tmp.add_atom(at.x, at.y)
    outer = []
    for n, _b in doc.neighbors(a):
        o = doc.atoms[n]
        k = tmp.add_atom(o.x, o.y)
        tmp.add_bond(x, k)
        outer.append(k)
    tmp.put_group(x, model.NICKNAMES[at.label])
    inner_id = {k: next(ids) for k in tmp.atoms}
    parts = []
    for k, t in tmp.atoms.items():
        if k in outer:
            parts.append(f'<n id="{inner_id[k]}" p="{pt(t.x, t.y)}" NodeType="ExternalConnectionPoint"/>')
            continue
        attrs = f'id="{inner_id[k]}" p="{pt(t.x, t.y)}"'
        if t.el != 'C':
            attrs += f' Element="{PT.GetAtomicNumber(t.el)}"'
        if t.charge:
            attrs += f' Charge="{t.charge}"'
        parts.append(f'<n {attrs}/>')
    for bd in tmp.bonds.values():
        order = f' Order="{bd.order}"' if bd.order != 1 else ''
        parts.append(f'<b id="{next(ids)}" B="{inner_id[bd.a]}" E="{inner_id[bd.b]}"{order}/>')
    return (f'<n id="{nid}" p="{pt(at.x, at.y)}" NodeType="Nickname" NeedsClean="yes"{cattr}>'
            f'<fragment id="{next(ids)}">{"".join(parts)}</fragment>{label_xml}</n>')
