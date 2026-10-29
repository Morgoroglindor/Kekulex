"""
canvas.py - the drawing area: the tools, the mouse and keyboard, selection, undo / redo.

    c = Canvas(); c.set_tool('bond1'); c.doc (the Doc being edited); c.changed (signal after every edit)

Tools: select, bond1 / bond2 / bond3 / wedge / hash, chain, ring3 ... ring8 / benzene, atom (with c.element),
charge+ / charge-, erase, arrow / equilibrium, text.
Mouse: click or drag with a bond tool (from an atom, a bond or empty space); ring tools fuse onto a bond or an atom;
right / middle drag pans, Ctrl + wheel zooms, the wheel scrolls.
Keyboard over an atom: C N O S P F H I, L = Cl, B = Br, + / - charge, Delete; over a bond: 1 2 3 = order.
"""
import math

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter, QPen, QBrush, QPainterPath
from PyQt6.QtWidgets import QInputDialog, QWidget

import model
import render
import theme
from model import Doc

RING_TOOLS = {'ring3': 3, 'ring4': 4, 'ring5': 5, 'ring6': 6, 'ring7': 7, 'ring8': 8, 'benzene': 6}
BOND_TOOLS = {'bond1': (1, ''), 'bond2': (2, ''), 'bond3': (3, ''), 'wedge': (1, 'wedge'), 'hash': (1, 'hash')}
KEY_ELEMENTS = {Qt.Key.Key_C: 'C', Qt.Key.Key_N: 'N', Qt.Key.Key_O: 'O', Qt.Key.Key_S: 'S', Qt.Key.Key_P: 'P',
                Qt.Key.Key_F: 'F', Qt.Key.Key_H: 'H', Qt.Key.Key_I: 'I', Qt.Key.Key_L: 'Cl', Qt.Key.Key_B: 'Br'}
DRAG_PX = 4                                                 # the mouse must move this far for a press to be a drag
HINTS = {
    'select': 'Click or drag a box to select (Shift adds); drag the selection to move it; double-click an atom to type its label.',
    'bond1': 'Click an atom to add a bond, drag to place it; click a bond to make it double, then triple.',
    'bond2': 'Click or drag to draw a double bond; click a double bond to move its second line.',
    'bond3': 'Click or drag to draw a triple bond.',
    'wedge': 'Wedge (toward you): the narrow end is where you start. Click a wedge again to reverse it.',
    'hash': 'Hashed bond (away from you): the narrow end is where you start. Click it again to reverse it.',
    'chain': 'Drag to draw a zigzag chain; the number of bonds follows the mouse.',
    'ring': 'Click empty space for a ring, a bond to fuse a ring onto it, an atom to attach a ring there.',
    'atom': 'Click an atom to change its element, or empty space to place one. Keys C N O S P F H I L(Cl) B(Br) work with any tool.',
    'charge+': 'Click an atom to add a positive charge (+ key with any tool).',
    'charge-': 'Click an atom to add a negative charge (- key with any tool).',
    'erase': 'Click (or drag over) atoms, bonds, arrows and text to delete them.',
    'arrow': 'Drag to draw a reaction arrow.',
    'equilibrium': 'Drag to draw an equilibrium arrow.',
    'text': 'Click to add text (captions, compound numbers); click text to edit it.',
}


class Canvas(QWidget):
    changed = pyqtSignal()                                  # the drawing was edited
    selection_changed = pyqtSignal()
    zoom_changed = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.doc = Doc()
        self.tool = 'bond1'
        self.element = 'N'
        self.zoom = 2.0                                     # 1.0 = real size on the screen
        self.origin = QPointF(200, 200)                     # where model (0, 0) is on the widget, px
        self.sel = set()                                    # {(kind, id)}
        self.hover = None                                   # (kind, id)
        self.undo_stack, self.redo_stack = [], []
        self.press = None                                   # dict while a mouse button is down
        self.preview = None                                 # what a drag would add: ('bond' | 'chain' | 'arrow' ...)
        self.mouse = QPointF()
        self._renderer = None
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(300, 200)
        self.setAutoFillBackground(False)

    # ------------------------------------------------------------------------------------------- coordinates
    def scale(self):
        """Pixels per bond length."""
        return self.zoom * self.logicalDpiX() / 72.0 * render.style_of(self.doc)['bond_length']

    def px_per_pt(self):
        return self.zoom * self.logicalDpiX() / 72.0

    def to_model(self, pos):
        s = self.scale()
        return (pos.x() - self.origin.x()) / s, (pos.y() - self.origin.y()) / s

    def to_px(self, x, y):
        s = self.scale()
        return QPointF(self.origin.x() + x * s, self.origin.y() + y * s)

    def renderer(self):
        if self._renderer is None:
            self._renderer = render.Renderer(self.doc)
        return self._renderer

    # ------------------------------------------------------------------------------------------- editing
    def set_doc(self, doc):
        self.doc = doc
        self.sel.clear(); self.hover = None
        self.undo_stack.clear(); self.redo_stack.clear()
        self.edited(record=False)
        self.fit()

    def snapshot(self):
        self.undo_stack.append(self.doc.copy())
        del self.undo_stack[:-200]
        self.redo_stack.clear()

    def edited(self, record=True):
        """Call after the drawing changed (snapshot() first)."""
        self._renderer = None
        self.sel = {(k, i) for k, i in self.sel if i in self.doc.table(k)}
        if self.hover and self.hover[1] not in self.doc.table(self.hover[0]):
            self.hover = None
        self.update()
        self.changed.emit()
        self.selection_changed.emit()

    def undo(self):
        if self.undo_stack:
            self.redo_stack.append(self.doc)
            self.doc = self.undo_stack.pop()
            self.edited()

    def redo(self):
        if self.redo_stack:
            self.undo_stack.append(self.doc)
            self.doc = self.redo_stack.pop()
            self.edited()

    def set_tool(self, tool, element=None):
        self.tool = tool
        if element:
            self.element = element
        self.preview = None
        self.setCursor(Qt.CursorShape.ArrowCursor if tool == 'select' else Qt.CursorShape.CrossCursor)
        self.update()

    def hint(self):
        return HINTS['ring' if self.tool in RING_TOOLS else self.tool]

    def selected_atoms(self):
        return self.doc.atoms_of(self.sel)

    def set_selection(self, items):
        self.sel = set(items)
        self.update()
        self.selection_changed.emit()

    def select_all(self):
        self.set_selection(self.doc.all_items())

    def delete_selection(self):
        if not self.sel:
            return
        self.snapshot()
        for kind, i in sorted(self.sel, key=lambda t: t[0] != 'bond'):   # bonds first, then atoms and the rest
            if i in self.doc.table(kind):
                self.doc.remove(kind, i)
        self.sel.clear()
        self.edited()

    def insert_items(self, items_dict=None, mols=()):
        """Paste: a to_dict() result and / or RDKit molecules, placed right of what is there; selects the new items."""
        self.snapshot()
        r = render.Renderer(self.doc).bounds()
        BL = render.style_of(self.doc)['bond_length']
        x = (r.right() / BL + 3) if not r.isNull() else 0.0
        y = (r.center().y() / BL) if not r.isNull() else 0.0
        new = []
        if items_dict:
            tmp = Doc.from_dict(items_dict)
            b = render.Renderer(tmp).bounds()
            if not b.isNull():
                dx, dy = x - b.left() / BL, y - b.center().y() / BL
            else:
                dx = dy = 0.0
            new += self.doc.insert_dict(items_dict, dx, dy)
            x += (b.width() / BL + 2) if not b.isNull() else 0
        for mol in mols:
            items = self.doc.insert_rdkit(mol, (0, 0))
            pts = self.doc.points_of(items)
            minx = min(p[0] for p in pts) if pts else 0
            self.doc.transform(items, lambda px, py, dx=x - minx: (px + dx, py + y))
            x = max(p[0] for p in self.doc.points_of(items)) + 2 if pts else x
            new += items
        self.sel = set(new)
        self.edited()
        self.ensure_visible()

    def transform_selection(self, fn_name):
        """'flip_h' | 'flip_v' | 'rot+' | 'rot-' on the selection (everything when nothing is selected)."""
        items = list(self.sel) or self.doc.all_items()
        pts = self.doc.points_of(items)
        if not pts:
            return
        cx = (min(p[0] for p in pts) + max(p[0] for p in pts)) / 2
        cy = (min(p[1] for p in pts) + max(p[1] for p in pts)) / 2
        self.snapshot()
        if fn_name == 'flip_h':
            fn = lambda x, y: (2 * cx - x, y)
        elif fn_name == 'flip_v':
            fn = lambda x, y: (x, 2 * cy - y)
        else:
            a = math.radians(15 if fn_name == 'rot+' else -15)
            c, s = math.cos(a), math.sin(a)
            fn = lambda x, y: (cx + (x - cx) * c - (y - cy) * s, cy + (x - cx) * s + (y - cy) * c)
        self.doc.transform(items, fn)
        if fn_name in ('flip_h', 'flip_v'):                 # wedges and hashes swap, so the flipped drawing is
            for kind, i in items:                           # still the same molecule (not its mirror image)
                if kind == 'bond' and self.doc.bonds[i].stereo:
                    b = self.doc.bonds[i]
                    b.stereo = 'hash' if b.stereo == 'wedge' else 'wedge'
        self.edited()

    def clean_up(self):
        """Lay out the selected structures (all, when nothing is selected) again with RDKit, keeping their place."""
        import chem
        from rdkit.Chem import rdDepictor
        from rdkit import Chem
        atoms = self.selected_atoms() or set(self.doc.atoms)
        frags = [set(c) for c in self.doc.fragments() if atoms & set(c)]   # whole structures, even if half selected
        if not frags:
            return
        self.snapshot()
        for frag in frags:
            mol, order, _err = chem.to_mol(self.doc, frag)
            if mol is None:
                continue
            ox = sum(self.doc.atoms[a].x for a in order) / len(order)
            oy = sum(self.doc.atoms[a].y for a in order) / len(order)
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(mol)
            Chem.WedgeMolBonds(mol, mol.GetConformer())
            conf = mol.GetConformer()
            pos = [conf.GetAtomPosition(k) for k in range(mol.GetNumAtoms())]
            lens = sorted(math.hypot(pos[b.GetBeginAtomIdx()].x - pos[b.GetEndAtomIdx()].x,
                                     pos[b.GetBeginAtomIdx()].y - pos[b.GetEndAtomIdx()].y) for b in mol.GetBonds())
            L = (lens[len(lens) // 2] if lens else 1.5) or 1.5
            cx = sum(p.x for p in pos) / len(pos); cy = sum(p.y for p in pos) / len(pos)
            for k, a in enumerate(order):
                self.doc.atoms[a].x = ox + (pos[k].x - cx) / L
                self.doc.atoms[a].y = oy - (pos[k].y - cy) / L
            for b in mol.GetBonds():                         # wedges where RDKit puts them (same stereo)
                i = self.doc.bond_between(order[b.GetBeginAtomIdx()], order[b.GetEndAtomIdx()])
                if i is None:
                    continue
                bd = self.doc.bonds[i]
                d = {Chem.BondDir.BEGINWEDGE: 'wedge', Chem.BondDir.BEGINDASH: 'hash'}.get(b.GetBondDir(), '')
                bd.stereo = d if bd.order == 1 else ''
                if d:
                    bd.a, bd.b = order[b.GetBeginAtomIdx()], order[b.GetEndAtomIdx()]
                bd.side = None
        self.edited()

    # ------------------------------------------------------------------------------------------- view
    def set_zoom(self, z, anchor=None):
        z = min(max(z, 0.25), 12.0)
        anchor = anchor if anchor is not None else QPointF(self.width() / 2, self.height() / 2)
        mx, my = self.to_model(anchor)
        self.zoom = z
        p = self.to_px(mx, my)
        self.origin += anchor - p
        self.update()
        self.zoom_changed.emit(z)

    def fit(self):
        """Show the whole drawing (at most 3x; an empty drawing: 2x with the origin near the middle)."""
        b = render.Renderer(self.doc).bounds()
        if b.isNull():
            self.zoom = 2.0
            self.origin = QPointF(self.width() / 2, self.height() / 2)
        else:
            k = self.logicalDpiX() / 72.0
            z = min((self.width() - 60) / max(b.width() * k, 1), (self.height() - 60) / max(b.height() * k, 1), 3.0)
            self.zoom = max(z, 0.25)
            c = b.center()
            self.origin = QPointF(self.width() / 2 - c.x() * self.px_per_pt(), self.height() / 2 - c.y() * self.px_per_pt())
        self.update()
        self.zoom_changed.emit(self.zoom)

    def ensure_visible(self):
        b = render.Renderer(self.doc).bounds()
        if b.isNull():
            return
        k = self.px_per_pt()
        r = QRectF(self.origin.x() + b.left() * k, self.origin.y() + b.top() * k, b.width() * k, b.height() * k)
        if not QRectF(self.rect()).contains(r):
            self.fit()

    # ------------------------------------------------------------------------------------------- hit testing
    def item_at(self, pos):
        x, y = self.to_model(pos)
        s = self.scale()
        r_atom = max(0.22, 7 / s)
        a = self.doc.find_atom(x, y, r_atom)
        if a is not None:
            return ('atom', a)
        rend = self.renderer()
        BL = rend.BL
        for i in self.doc.texts:
            if rend.text_rect(i).adjusted(-1, -1, 1, 1).contains(QPointF(x * BL, y * BL)):
                return ('text', i)
        for a2 in self.doc.atoms:                           # a label wider than the atom circle (NH2, OH)
            r = rend.label_rect(a2)
            if r is not None and r.adjusted(-1, -1, 1, 1).contains(QPointF(x * BL, y * BL)):
                return ('atom', a2)
        b = self.doc.find_bond(x, y, max(0.14, 5 / s))
        if b is not None:
            return ('bond', b)
        ar = self.doc.find_arrow(x, y, max(0.15, 6 / s))
        if ar is not None:
            return ('arrow', ar)
        return None

    # ------------------------------------------------------------------------------------------- painting
    def paintEvent(self, _e):
        p = QPainter(self)
        p.fillRect(self.rect(), QColor('#FFFFFF'))
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.translate(self.origin)
        k = self.px_per_pt()
        p.scale(k, k)
        rend = self.renderer()
        BL = rend.BL
        accent = QColor(theme.current_accent)
        soft = QColor(accent); soft.setAlpha(60)
        selc = QColor(accent); selc.setAlpha(90)
        for kind, i in self.sel:
            self._mark(p, rend, kind, i, selc, BL)
        if self.hover and self.hover not in self.sel:
            self._mark(p, rend, *self.hover, soft, BL)
        rend.paint(p, errors=True)
        if self.preview:
            self._paint_preview(p, BL, accent)
        if self.press and self.press.get('band'):
            p.resetTransform()
            pen = QPen(accent, 1, Qt.PenStyle.DashLine)
            p.setPen(pen)
            fill = QColor(accent); fill.setAlpha(25)
            p.setBrush(fill)
            p.drawRect(QRectF(self.press['pos'], self.mouse).normalized())
        p.end()

    def _mark(self, p, rend, kind, i, color, BL):
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color))
        if kind == 'atom' and i in self.doc.atoms:
            at = self.doc.atoms[i]
            r = rend.label_rect(i)
            if r is not None and rend.labels[i]['clip']:
                p.drawRoundedRect(r.adjusted(-1.5, -1.5, 1.5, 1.5), 2, 2)
            else:
                p.drawEllipse(QPointF(at.x * BL, at.y * BL), 0.2 * BL, 0.2 * BL)
        elif kind == 'bond' and i in self.doc.bonds:
            b = self.doc.bonds[i]
            pa, pb = self.doc.atoms[b.a], self.doc.atoms[b.b]
            p.setPen(QPen(color, 0.32 * BL, cap=Qt.PenCapStyle.RoundCap)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawLine(QPointF(pa.x * BL, pa.y * BL), QPointF(pb.x * BL, pb.y * BL))
        elif kind == 'arrow' and i in self.doc.arrows:
            a = self.doc.arrows[i]
            p.setPen(QPen(color, 0.3 * BL, cap=Qt.PenCapStyle.RoundCap)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawLine(QPointF(a.x1 * BL, a.y1 * BL), QPointF(a.x2 * BL, a.y2 * BL))
        elif kind == 'text' and i in self.doc.texts:
            p.drawRoundedRect(rend.text_rect(i).adjusted(-2, -1, 2, 1), 2, 2)

    def _paint_preview(self, p, BL, accent):
        kind, pts = self.preview
        pen = QPen(accent, 1.2, cap=Qt.PenCapStyle.RoundCap, join=Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
        path = QPainterPath()
        if kind in ('bond', 'chain', 'arrow', 'equilibrium', 'ring'):
            path.moveTo(pts[0][0] * BL, pts[0][1] * BL)
            for x, y in pts[1:]:
                path.lineTo(x * BL, y * BL)
            if kind == 'ring':
                path.closeSubpath()
        p.drawPath(path)
        if kind == 'chain':
            p.resetTransform()
            p.setPen(accent)
            p.drawText(self.mouse + QPointF(12, -8), str(len(pts) - 1))

    # ------------------------------------------------------------------------------------------- mouse
    def mousePressEvent(self, e):
        self.setFocus()
        if e.button() in (Qt.MouseButton.MiddleButton, Qt.MouseButton.RightButton):
            self.press = dict(pan=True, pos=e.position(), origin=QPointF(self.origin))
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            return
        if e.button() != Qt.MouseButton.LeftButton:
            return
        pos = e.position()
        item = self.item_at(pos)
        self.press = dict(pos=pos, item=item, moved=False, shift=bool(e.modifiers() & Qt.KeyboardModifier.ShiftModifier))
        if self.tool == 'select':
            if item and item in self.sel:
                self.press['move'] = list(self.sel)
            elif item:
                if not self.press['shift']:
                    self.sel = set()
                self.sel.add(item)
                self.press['move'] = list(self.sel)
                self.selection_changed.emit()
            else:
                self.press['band'] = True
                if not self.press['shift']:
                    self.sel = set()
                    self.selection_changed.emit()
            if 'move' in self.press:
                self.press['start'] = self.to_model(pos)
                self.press['snap'] = None
            self.update()
        elif self.tool == 'erase':
            self.snapshot()
            self.press['erased'] = False
            self._erase(item)

    def _erase(self, item):
        if item and item[1] in self.doc.table(item[0]):
            self.doc.remove(*item)
            self.press['erased'] = True
            self.sel.discard(item)
            self.hover = None
            self.edited()

    def mouseMoveEvent(self, e):
        pos = e.position()
        self.mouse = pos
        pr = self.press
        if pr and pr.get('pan'):
            self.origin = pr['origin'] + (pos - pr['pos'])
            self.update()
            return
        if pr is None:
            item = self.item_at(pos)
            if item != self.hover:
                self.hover = item
                self.update()
            return
        if not pr['moved'] and (pos - pr['pos']).manhattanLength() < DRAG_PX:
            return
        pr['moved'] = True
        x, y = self.to_model(pos)
        t = self.tool
        if t == 'select':
            if 'move' in pr:
                if not pr.get('snapped'):
                    pr['snapped'] = True
                    self.snapshot()
                sx, sy = pr['start']
                dx, dy = x - sx, y - sy
                self.doc.transform(pr['move'], lambda px, py: (px + dx, py + dy))
                pr['start'] = (x, y)
                self._renderer = None
                moved_atoms = self.doc.atoms_of(pr['move'])
                self.hover = None
                if len(moved_atoms) == 1:                   # a single atom dragged onto another merges on release
                    a = next(iter(moved_atoms))
                    at = self.doc.atoms[a]
                    tgt = self.doc.find_atom(at.x, at.y, model.MERGE, exclude={a} | {n for n, _b in self.doc.neighbors(a)})
                    self.hover = ('atom', tgt) if tgt is not None else None
                    pr['merge'] = (a, tgt) if tgt is not None else None
                self.update()
            elif pr.get('band'):
                self._band_select(pr['pos'], pos, pr['shift'])
            return
        if t == 'erase':
            self._erase(self.item_at(pos))
            return
        if t in BOND_TOOLS or t == 'chain':
            start = self._start_point(pr['item'])
            if start is None:
                return
            sx, sy = start
            if t == 'chain':
                self.preview = ('chain', [(sx, sy)] + self.doc.chain_points(sx, sy, x, y))
            else:
                ex, ey = self._bond_end(sx, sy, x, y, pr['item'])
                self.preview = ('bond', [(sx, sy), (ex, ey)])
            self.update()
        elif t in ('arrow', 'equilibrium'):
            sx, sy = self.to_model(pr['pos'])
            ang = model.snap_angle(model.angle_of(x - sx, y - sy), 15)
            L = math.hypot(x - sx, y - sy)
            self.preview = (t, [(sx, sy), (sx + L * math.cos(ang), sy + L * math.sin(ang))])
            self.update()

    def _start_point(self, item):
        if item and item[0] == 'atom':
            at = self.doc.atoms[item[1]]
            return at.x, at.y
        if item is None:
            return self.to_model(self.press['pos'])
        return None

    def _bond_end(self, sx, sy, x, y, item):
        """Where a dragged bond ends: on an atom under the mouse, else one bond length away at a 15 degree angle."""
        exclude = {item[1]} if item and item[0] == 'atom' else set()
        hit = self.doc.find_atom(x, y, 0.3, exclude)
        if hit is not None:
            return self.doc.atoms[hit].x, self.doc.atoms[hit].y
        ang = model.snap_angle(model.angle_of(x - sx, y - sy), 15)
        return sx + math.cos(ang), sy + math.sin(ang)

    def _band_select(self, p0, p1, add):
        r = QRectF(p0, p1).normalized()
        items = set()
        for i, at in self.doc.atoms.items():
            if r.contains(self.to_px(at.x, at.y)):
                items.add(('atom', i))
        for i, b in self.doc.bonds.items():
            if ('atom', b.a) in items and ('atom', b.b) in items:
                items.add(('bond', i))
        for i, a in self.doc.arrows.items():
            if r.contains(self.to_px(a.x1, a.y1)) and r.contains(self.to_px(a.x2, a.y2)):
                items.add(('arrow', i))
        k = self.px_per_pt()
        for i in self.doc.texts:
            tr = self.renderer().text_rect(i)
            if r.intersects(QRectF(self.origin.x() + tr.left() * k, self.origin.y() + tr.top() * k, tr.width() * k, tr.height() * k)):
                items.add(('text', i))
        self.sel = (self.sel | items) if add else items
        self.update()
        self.selection_changed.emit()

    def mouseReleaseEvent(self, e):
        pr, self.press = self.press, None
        if pr is None:
            return
        if pr.get('pan'):
            self.set_tool(self.tool)
            return
        pos = e.position()
        x, y = self.to_model(pos)
        t, item = self.tool, pr['item']
        preview, self.preview = self.preview, None
        if t == 'select':
            if pr.get('merge') and pr.get('snapped'):
                a, tgt = pr['merge']
                self.doc.merge_atom(a, tgt)
                self.sel.discard(('atom', a))
                self.edited()
            elif pr.get('snapped'):
                self.edited()
            elif not pr['moved'] and item is None and not pr['shift']:
                self.set_selection([])
            self.update()
            return
        if t == 'erase':
            if not pr.get('erased'):
                self.undo_stack.pop()
            return
        if pr['moved'] and preview:
            self.snapshot()
            kind, pts = preview
            if kind == 'bond':
                a = item[1] if item and item[0] == 'atom' else self.doc.atom_at_or_new(*pts[0])
                b = self.doc.atom_at_or_new(*pts[1], exclude={a})
                order, stereo = BOND_TOOLS[t]
                bid = self.doc.add_bond(a, b, order, stereo)
                if bid is not None and self.doc.bonds[bid].a != a and stereo:
                    bd = self.doc.bonds[bid]; bd.a, bd.b = bd.b, bd.a
            elif kind == 'chain':
                prev = item[1] if item and item[0] == 'atom' else self.doc.atom_at_or_new(*pts[0])
                for px, py in pts[1:]:
                    nxt = self.doc.atom_at_or_new(px, py, exclude={prev})
                    self.doc.add_bond(prev, nxt)
                    prev = nxt
            elif kind in ('arrow', 'equilibrium'):
                (x1, y1), (x2, y2) = pts
                if math.hypot(x2 - x1, y2 - y1) > 0.3:
                    self.doc.add_arrow(x1, y1, x2, y2, kind)
            self.edited()
            return
        if not pr['moved']:
            self.click(item, x, y)

    def click(self, item, x, y):
        """A click (no drag) with the current tool on item (None: empty space at model x, y)."""
        t, doc = self.tool, self.doc
        if t in ('arrow', 'equilibrium'):
            self.snapshot(); doc.add_arrow(x, y, x + 2.5, y, t); self.edited(); return
        if t == 'text':
            self.edit_text(item[1] if item and item[0] == 'text' else None, x, y)
            return
        if item and item[0] == 'text':
            return
        if t in BOND_TOOLS or t == 'chain':
            order, stereo = BOND_TOOLS.get(t, (1, ''))
            self.snapshot()
            if item is None:
                a = doc.add_atom(x - math.cos(math.radians(30)) / 2, y + 0.25)
                b = doc.atom_at_or_new(doc.atoms[a].x + math.cos(math.radians(30)), doc.atoms[a].y - 0.5, exclude={a})
                doc.add_bond(a, b, order, stereo)
            elif item[0] == 'atom':
                a = item[1]
                ang = doc.next_bond_angle(a)
                at = doc.atoms[a]
                b = doc.atom_at_or_new(at.x + math.cos(ang), at.y + math.sin(ang), exclude={a})
                bid = doc.add_bond(a, b, order, stereo)
                if bid is not None and stereo and doc.bonds[bid].a != a:
                    bd = doc.bonds[bid]; bd.a, bd.b = bd.b, bd.a
            elif item[0] == 'bond':
                bd = doc.bonds[item[1]]
                if stereo:
                    if bd.stereo == stereo:
                        bd.a, bd.b = bd.b, bd.a
                    bd.stereo, bd.order = stereo, 1
                elif t == 'bond1' or t == 'chain':
                    bd.stereo = ''
                    bd.order = 1 if bd.order == 3 else bd.order + 1
                    bd.side = None
                elif t == 'bond2':
                    if bd.order == 2 and not bd.stereo:     # cycle the second line: auto -> other side -> centred -> auto
                        cur = doc.double_side(item[1])
                        if bd.side is None:
                            bd.side = -cur if cur else 1
                        else:
                            bd.side = 0 if bd.side else None
                    bd.order, bd.stereo = 2, ''
                else:
                    bd.order, bd.stereo = 3, ''
            else:
                self.undo_stack.pop(); return
            self.edited()
        elif t in RING_TOOLS:
            n = RING_TOOLS[t]
            self.snapshot()
            if item is None:
                pts = doc.ring_points(n, (x, y), math.radians(-90))
            elif item[0] == 'bond':
                pts = doc.ring_on_bond(item[1], n)
            elif item[0] == 'atom':
                pts = doc.ring_on_atom(item[1], n)
            else:
                self.undo_stack.pop(); return
            doc.add_ring(pts, aromatic=(t == 'benzene'))
            self.edited()
        elif t == 'atom':
            self.snapshot()
            if item is None:
                doc.add_atom(x, y, self.element)
            elif item[0] == 'atom':
                doc.atoms[item[1]].el = self.element
            else:
                self.undo_stack.pop(); return
            self.edited()
        elif t in ('charge+', 'charge-') and item and item[0] == 'atom':
            self.snapshot()
            doc.atoms[item[1]].charge += 1 if t == 'charge+' else -1
            self.edited()

    def mouseDoubleClickEvent(self, e):
        item = self.item_at(e.position())
        if item and item[0] == 'atom':
            self.edit_label(item[1])
        elif item and item[0] == 'text':
            self.edit_text(item[1])

    def wheelEvent(self, e):
        d = e.angleDelta()
        if e.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.set_zoom(self.zoom * (1.15 ** (d.y() / 120)), e.position())
        elif e.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.origin += QPointF(d.y() / 2, 0); self.update()
        else:
            self.origin += QPointF(d.x() / 2, d.y() / 2); self.update()

    def leaveEvent(self, _e):
        if self.hover:
            self.hover = None
            self.update()

    # ------------------------------------------------------------------------------------------- keyboard
    def keyPressEvent(self, e):
        key = e.key()
        mods = e.modifiers() & ~Qt.KeyboardModifier.KeypadModifier
        h = self.hover
        if key == Qt.Key.Key_Escape:
            self.press = None; self.preview = None
            self.set_selection([])
            return
        if key in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            if h and h not in self.sel:
                self.snapshot(); self.doc.remove(*h); self.hover = None; self.edited()
            else:
                self.delete_selection()
            return
        if key in (Qt.Key.Key_Left, Qt.Key.Key_Right, Qt.Key.Key_Up, Qt.Key.Key_Down) and self.sel:
            step = 0.5 if mods & Qt.KeyboardModifier.ShiftModifier else 0.05
            dx = {Qt.Key.Key_Left: -step, Qt.Key.Key_Right: step}.get(key, 0)
            dy = {Qt.Key.Key_Up: -step, Qt.Key.Key_Down: step}.get(key, 0)
            self.snapshot()
            self.doc.transform(list(self.sel), lambda x, y: (x + dx, y + dy))
            self.edited()
            return
        if mods & ~Qt.KeyboardModifier.ShiftModifier:
            return super().keyPressEvent(e)
        if h and h[0] == 'atom':
            at = self.doc.atoms[h[1]]
            if key in KEY_ELEMENTS:
                el = KEY_ELEMENTS[key]
                self.snapshot(); at.el = el; self.edited(); return
            if key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal, Qt.Key.Key_Minus, Qt.Key.Key_Underscore):
                self.snapshot()
                at.charge += 1 if key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal) else -1
                self.edited(); return
        if h and h[0] == 'bond' and key in (Qt.Key.Key_1, Qt.Key.Key_2, Qt.Key.Key_3):
            self.snapshot()
            bd = self.doc.bonds[h[1]]
            bd.order, bd.stereo, bd.side = {Qt.Key.Key_1: 1, Qt.Key.Key_2: 2, Qt.Key.Key_3: 3}[key], '', None
            self.edited(); return
        super().keyPressEvent(e)

    # ------------------------------------------------------------------------------------------- dialogs
    def edit_label(self, a):
        at = self.doc.atoms[a]
        cur = at.el + (render.Renderer._charge_text(at.charge) if at.charge else '')
        text, ok = QInputDialog.getText(self, 'Atom', 'Element (e.g. N, Cl, Na+, O-):', text=cur)
        if not ok or not text.strip():
            return
        el, charge = parse_label(text.strip())
        if el is None:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, 'Atom', f'"{text}" is not an element symbol.')
            return
        self.snapshot()
        at.el, at.charge = el, charge
        self.edited()

    def edit_text(self, i=None, x=0.0, y=0.0):
        cur = self.doc.texts[i].text if i is not None else ''
        text, ok = QInputDialog.getMultiLineText(self, 'Text', 'Text:', cur)
        if not ok:
            return
        self.snapshot()
        if i is None:
            if text.strip():
                self.doc.add_text(x, y, text)
        elif text.strip():
            self.doc.texts[i].text = text
        else:
            del self.doc.texts[i]
        self.edited()


def parse_label(text):
    """'Na+' -> ('Na', 1), 'o' -> ('O', 0), 'N2+' -> ('N', 2); (None, 0) when not an element."""
    import re
    from rdkit import Chem
    m = re.fullmatch(r'([A-Za-z]{1,2})\s*(\d*)\s*([+\-−]*)', text)
    if not m:
        return None, 0
    sym = m.group(1)[0].upper() + m.group(1)[1:].lower()
    try:
        if Chem.GetPeriodicTable().GetAtomicNumber(sym) <= 0:
            return None, 0
    except RuntimeError:
        return None, 0
    signs = m.group(3).replace('−', '-')
    n = int(m.group(2)) if m.group(2) else (len(signs) or 0)
    if not signs:
        return sym, 0
    return sym, n if signs[0] == '+' else -n
