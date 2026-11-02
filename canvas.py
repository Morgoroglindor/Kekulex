"""
canvas.py - the drawing area: the page, the tools, the mouse and keyboard, typing labels and captions in place,
selection, arranging (align, distribute, rotate, scale, flip, colour), undo / redo.

    c = Canvas(); c.set_tool('bond1'); c.doc (the Doc being edited); c.changed (signal after every edit)

Tools: select, bond tools (BOND_TOOLS), chain, ring3 ... ring8 / benzene / cyclopentadiene, atom (with c.element),
charge+ / charge-, erase, arrow tools (ARROW_TOOLS), shape tools (SHAPE_TOOLS), text.
Mouse: click or drag with a bond tool (from an atom, a bond or empty space); ring tools fuse onto a bond or an atom;
right / middle drag pans, Ctrl + wheel zooms, the wheel scrolls. Double-click: an atom's label or a caption is typed in
place; a bond selects its whole structure.
Keyboard: ChemDraw's hotkeys (ATOM_KEYS, SPROUT_KEYS, FUSE_KEYS, BOND_STYLE_KEYS, TOOL_KEYS below), + / - charge, Delete.
"""
import math

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QLineEdit, QPlainTextEdit, QWidget

import model
import render
import theme
from model import Doc

RING_TOOLS = {'ring3': (3, False), 'ring4': (4, False), 'ring5': (5, False), 'ring6': (6, False), 'ring7': (7, False),
              'ring8': (8, False), 'benzene': (6, True), 'cyclopentadiene': (5, True)}
BOND_TOOLS = {'bond1': (1, '', ''), 'bond2': (2, '', ''), 'bond3': (3, '', ''), 'wedge': (1, 'wedge', ''),
              'hash': (1, 'hash', ''), 'bold': (1, '', 'bold'), 'dashed': (1, '', 'dashed'), 'wavy': (1, '', 'wavy'),
              'hashed': (1, '', 'hashed')}
ARROW_TOOLS = {'arrow': 'forward', 'equilibrium': 'equilibrium', 'retro': 'retro', 'resonance': 'resonance',
               'arrow_dashed': 'dashed', 'curved': 'curved', 'curved_half': 'curved_half'}
# shape-like tools: name -> (Shape kind, style, params). Drawing elements are 'kind:style' (circle:shaded,
# line:wavy, arc:120:dashed), orbitals 'orb:type:look' (orb:p:shaded)
SHAPE_TOOLS = {k: (k, '', {}) for k in model.SHAPE_KINDS + model.BIO_KINDS}
ELEMENT_TOOLS = []                                          # in the palette's order: 8 columns, as ChemDraw's
for _row, (_kind, _line, _ang) in enumerate((('circle', 'line', 270), ('oval', 'line:dashed', 180),
                                             ('roundrect', 'line:wavy', 120), ('rect', 'line:bold', 90))):
    for _st in ('', 'dashed', 'shaded', 'filled', 'shadow'):
        _name = f'{_kind}:{_st}' if _st else _kind
        SHAPE_TOOLS[_name] = (_kind, _st, {})
        ELEMENT_TOOLS.append(_name)
    SHAPE_TOOLS[_line] = ('line', _line.partition(':')[2], {})
    ELEMENT_TOOLS.append(_line)
    for _st in ('', 'dashed'):
        _name = f'arc:{_ang}' + (':dashed' if _st else '')
        SHAPE_TOOLS[_name] = ('arc', _st, {'angle': _ang})
        ELEMENT_TOOLS.append(_name)
ORBITAL_TOOLS = []                                          # 3 columns: plain, shaded, filled
for _typ in ('s', 'sigma', 'lobe', 'p', 'hybrid', 'dxy', 'dz2'):
    for _look in ('plain', 'shaded', 'filled'):
        SHAPE_TOOLS[f'orb:{_typ}:{_look}'] = ('orbital', '' if _look == 'plain' else _look, {'type': _typ})
        ORBITAL_TOOLS.append(f'orb:{_typ}:{_look}')
CUT_TOOLS = ('massfrag', 'retrosyn', 'synthesis')           # drag a line across bonds
CURVE = 0.3                                                 # a new curved arrow's bend
# Hotkeys laid out as in ChemDraw (case matters: Shift+key is the capital), so the same fingers work in both.
# Over an atom - the atom becomes this element or nickname:
ATOM_KEYS = {'c': 'C', 'C': 'Cl', 'n': 'N', 'w': 'N', 'o': 'O', 'q': 'O', 's': 'S', 'S': 'Si', 'p': 'P', 'f': 'F',
             'h': 'H', 'i': 'I', 'l': 'Cl', 'b': 'Br', 'B': 'B', 'L': 'Li', 'm': 'Me', 'e': 'Et', 'A': 'Ac',
             'E': 'CO2Me', 'F': 'CF3', 'H': 'Cbz', 'M': 'MgBr', 'N': 'NO2', 'O': 'OMe', 'Q': 'Fmoc', 'P': 'Ph',
             'y': 'Boc', 'Z': 'N3'}
# Over an atom - a new bond from it to (group, bond order, stereo); '9' adds two methyls
SPROUT_KEYS = {'1': ('C', 1, ''), '0': ('C', 1, ''), '2': ('O', 2, ''), '3': ('Ph', 1, ''), 'a': ('Ph', 1, ''),
               'J': ('Ph', 1, ''), '4': ('C', 1, 'wedge'), '5': ('C', 1, 'hash'), '6': ('Cy', 1, ''),
               '7': ('cPent', 1, ''), '8': ('C', 2, ''), 'z': ('ethynyl', 1, ''), 'v': ('cPr', 1, ''),
               'u': ('cBu', 1, ''), 'K': ('tBu', 1, ''), 'k': ('SO2', 1, ''), 'j': ('Cp', 1, '')}
# Over a bond - a fused ring (size, aromatic); a drawing style
FUSE_KEYS = {'v': (3, False), '4': (4, False), '5': (5, False), '6': (6, False), '7': (7, False), '8': (8, False),
             'a': (6, True), 'z': (5, True)}
BOND_STYLE_KEYS = {'b': 'bold', 'B': 'bold', 'd': 'dashed', 'D': 'dashed', 'y': 'wavy', 'H': 'hashed'}
# Anywhere else - a tool
TOOL_KEYS = {'g': 'tlc', 'G': 'orb:p:shaded', 'j': 'benzene', 'J': 'cyclopentadiene', 'x': 'bond1', 'X': 'chain', 'e': 'arrow', 't': 'text',
             ' ': 'select'}
DRAG_PX = 4                                                 # the mouse must move this far for a press to be a drag
PAGE_BG, PAGE_SHADOW = QColor('#C8CCD1'), QColor(0, 0, 0, 40)
HINTS = {
    'lasso': 'Drag around what to select, free-hand (Shift adds); drag the selection to move it. Right-click: object menu.',
    'select':'Click or drag a box to select (Shift adds); drag the selection to move it; double-click an atom or text '
              'to type, a bond to select its structure.',
    'bond': 'Click an atom to add a bond, drag to place it; click a bond to change it (single -> double -> triple).',
    'chain': 'Drag to draw a zigzag chain; the number of bonds follows the mouse.',
    'ring': 'Click empty space for a ring, a bond to fuse a ring onto it, an atom to attach a ring there.',
    'atom': 'Click an atom to change its element, or empty space to place one. ChemDraw hotkeys (Ctrl+/) work with any tool.',
    'charge+': 'Click an atom to add a positive charge (+ key with any tool).',
    'charge-': 'Click an atom to add a negative charge (- key with any tool).',
    'erase': 'Click (or drag over) atoms, bonds, arrows, shapes and text to delete them.',
    'arrow': 'Drag to draw the arrow. Click a curved arrow with a curved-arrow tool to bend it the other way.',
    'shape': 'Drag to draw the shape (brackets come as a pair around what you frame).',
    'sequence': 'Click where the sequence goes, then type or paste it (FASTA, H-Ala-Gly-OH or HELM); click a sequence to edit it.',
    'orbital': 'Click an atom to put the orbital on it, or drag from the atom to set its size and direction.',
    'arc': 'Drag from end to end: left to right bulges up, right to left down.',
    'massfrag': 'Drag a line across bonds: the formula and exact mass of the fragments appear (Alt keeps the line).',
    'retrosyn': 'Drag a line across bonds of a product: its retrosynthesis (product, ⇒, the pieces) is drawn.',
    'synthesis': 'Drag a line across bonds of a product: the pieces, an arrow and the product are drawn.',
    'tlc': 'Drag a TLC plate (its width sets the lanes) or click; drag a spot to change its Rf (Shift: resize, '
           'Ctrl: copy to another lane); right-click for lanes and spots, double-click to type the Rf values.',
    'bio': 'Click for the default size, or drag (membranes, DNA, helices and mRNA grow more subunits as they get longer).',
    'text': 'Click empty space to type a caption, click an atom to type its label (N, OH, OMe, Boc, CO2H ...).',
}


class _LabelEdit(QLineEdit):
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.returnPressed.connect(lambda: canvas.commit_editor())

    def keyPressEvent(self, e):
        if e.key() == Qt.Key.Key_Escape:
            self.canvas.commit_editor(cancel=True)
            return
        super().keyPressEvent(e)

    def focusOutEvent(self, e):
        super().focusOutEvent(e)
        self.canvas.commit_editor()


class _TextEdit(QPlainTextEdit):
    def __init__(self, canvas):
        super().__init__(canvas)
        self.canvas = canvas
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.textChanged.connect(self.fit)

    def fit(self):
        fm = self.fontMetrics()
        lines = self.toPlainText().split('\n')
        w = max([fm.horizontalAdvance(line) for line in lines] + [60]) + 24
        h = fm.lineSpacing() * max(1, len(lines)) + 14
        self.resize(w, h)

    def keyPressEvent(self, e):
        if e.key() == Qt.Key.Key_Escape or (e.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter)
                                            and e.modifiers() & Qt.KeyboardModifier.ControlModifier):
            self.canvas.commit_editor()
            return
        super().keyPressEvent(e)

    def focusOutEvent(self, e):
        super().focusOutEvent(e)
        self.canvas.commit_editor()


MIN_ZOOM = 0.1                                              # low enough to see a whole A0 page


class Canvas(QWidget):
    changed = pyqtSignal()                                  # the drawing was edited
    selection_changed = pyqtSignal()
    zoom_changed = pyqtSignal(float)
    tool_key = pyqtSignal(str)                              # a tool hotkey was pressed (the window switches tools)
    message = pyqtSignal(str)                               # something to show in the status bar
    context_menu = pyqtSignal(object, object)               # right-click: (global QPoint, item under the mouse or None)
    edit_object = pyqtSignal(str, int)                      # double-click on a sequence or a plasmid map: (kind, id)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.doc = Doc()
        self.tool = 'bond1'
        self.element = 'N'
        self.zoom = 1.5                                     # 1.0 = real size on the screen
        self.origin = QPointF(40, 40)                       # where the page's top-left corner is on the widget, px
        self.sel = set()                                    # {(kind, id)}
        self.hover = None                                   # (kind, id)
        self.undo_stack, self.redo_stack = [], []
        self.press = None                                   # dict while a mouse button is down
        self.preview = None                                 # what a drag would add: (kind, points)
        self.mouse = QPointF()
        self.editor = None                                  # (widget, 'label' | 'text', id, x, y) while typing
        self.text_props = dict(align='left', size=0.0, bold=False, italic=False, formula=False, color='')
        self.show_reaction = False                          # View > Show reaction interpretation
        self.show_rulers = False
        self.show_warnings = True                           # red boxes on atoms over their valence
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

    def rect_px(self, r):
        """A QRectF in pt (model x bond length) as widget pixels."""
        k = self.px_per_pt()
        return QRectF(self.origin.x() + r.left() * k, self.origin.y() + r.top() * k, r.width() * k, r.height() * k)

    def page_size(self):
        return render.page_size(render.style_of(self.doc))

    def renderer(self):
        if self._renderer is None:
            self._renderer = render.Renderer(self.doc)
        return self._renderer

    # ------------------------------------------------------------------------------------------- editing
    def set_doc(self, doc):
        self.commit_editor(cancel=True)
        self.doc = doc
        self.place_on_page()
        self.sel.clear(); self.hover = None
        self.undo_stack.clear(); self.redo_stack.clear()
        self.edited(record=False)
        self.fit()

    def place_on_page(self):
        """Move a drawing that sticks out past the page's top or left edge (from a MOL / SMILES file) onto the page."""
        b = render.Renderer(self.doc).bounds()
        if b.isNull() or (b.left() >= 0 and b.top() >= 0):
            return
        BL = render.style_of(self.doc)['bond_length']
        dx, dy = (72 - b.left()) / BL, (72 - b.top()) / BL
        self.doc.transform(self.doc.all_items(), lambda x, y: (x + dx, y + dy))

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
        self.commit_editor()
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
        self.commit_editor()
        self.tool = tool
        if element:
            self.element = element
        self.preview = None
        self.setCursor(Qt.CursorShape.ArrowCursor if tool in ('select', 'lasso') else
                       Qt.CursorShape.IBeamCursor if tool == 'text' else Qt.CursorShape.CrossCursor)
        self.update()

    def hint(self):
        t = self.tool
        key = ('ring' if t in RING_TOOLS else 'bond' if t in BOND_TOOLS else 'arrow' if t in ARROW_TOOLS else
               'bio' if t in model.BIO_KINDS else 'orbital' if t.startswith('orb:') else
               'arc' if t.startswith('arc') else 'shape' if t in SHAPE_TOOLS else t)
        return HINTS.get(key, '')

    def selected_atoms(self):
        return self.doc.atoms_of(self.sel)

    def selected(self, kind):
        return [i for k, i in self.sel if k == kind and i in self.doc.table(kind)]

    def set_selection(self, items):
        self.sel = set(items)
        self.update()
        self.selection_changed.emit()

    def select_all(self):
        self.set_selection(self.doc.all_items())

    def select_fragment(self, atom):
        for comp in self.doc.fragments():
            if atom in comp:
                atoms = set(comp)
                items = [('atom', a) for a in atoms] + [('bond', i) for i, b in self.doc.bonds.items() if b.a in atoms]
                self.set_selection(items)
                return

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
        """Paste: a to_dict() result and / or RDKit molecules, placed right of what is there (on an empty page: near
        its top left); selects the new items."""
        self.snapshot()
        r = render.Renderer(self.doc).bounds()
        BL = render.style_of(self.doc)['bond_length']
        x = (r.right() / BL + 3) if not r.isNull() else 108 / BL
        y = (r.center().y() / BL) if not r.isNull() else 144 / BL
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

    # ------------------------------------------------------------------------------------------- biopolymers
    def add_sequence(self, values, pos=None):
        """A new sequence (model.Sequence fields in values) at pos (bond lengths), else right of what is drawn."""
        self.snapshot()
        if pos is None:
            r = render.Renderer(self.doc).bounds()
            BL = render.style_of(self.doc)['bond_length']
            pos = ((r.right() / BL + 3) if not r.isNull() else 108 / BL, (r.top() / BL) if not r.isNull() else 144 / BL)
        v = dict(values)
        i = self.doc.add_sequence(pos[0], pos[1], v.pop('polymer'), v.pop('residues'), **v)
        self.sel = {('seq', i)}
        self.edited()
        self.ensure_visible()
        return i

    def update_sequence(self, i, values):
        self.snapshot()
        for k, v in values.items():
            setattr(self.doc.seqs[i], k, v)
        self.edited()

    def set_sequence_props(self, **props):
        """Set fields on every selected sequence (one undo step), e.g. three_letter=True."""
        seqs = self.selected('seq')
        if not seqs:
            return
        self.snapshot()
        for i in seqs:
            for k, v in props.items():
                setattr(self.doc.seqs[i], k, v)
        self.edited()

    def expand_sequences(self):
        """Replace the selected sequences by their full structures (Structure > Expand label in ChemDraw)."""
        import biopolymer
        seqs = self.selected('seq')
        if not seqs:
            return
        mols = []
        for i in seqs:
            mol, err = biopolymer.to_mol(self.doc.seqs[i])
            if mol is None:
                self.message.emit(err)
                return
            mols.append((i, mol))
        self.snapshot()
        new = []
        for i, mol in mols:
            q = self.doc.seqs.pop(i)
            items = self.doc.insert_rdkit(mol, (0, 0))
            pts = self.doc.points_of(items)
            minx, miny = min(p[0] for p in pts), min(p[1] for p in pts)
            self.doc.transform(items, lambda px, py, dx=q.x - minx, dy=q.y - miny: (px + dx, py + dy))
            new += items
        self.sel = set(new)
        self.edited()
        self.ensure_visible()

    # ------------------------------------------------------------------------------------------- arranging
    def _targets(self):
        return list(self.sel) or self.doc.all_items()

    def _center(self, items):
        pts = self.doc.points_of(items)
        return ((min(p[0] for p in pts) + max(p[0] for p in pts)) / 2,
                (min(p[1] for p in pts) + max(p[1] for p in pts)) / 2) if pts else None

    def transform_selection(self, fn_name):
        """On the selection (everything when nothing is selected):
        'flip_h' / 'flip_v'      Rotate 180° Horizontal / Vertical: turned over, still the same molecule (wedges and
                                 hashes swap, as when a model is turned round)
        'mirror_h' / 'mirror_v'  Flip Horizontal / Vertical: the mirror image (stereocentres invert)
        'rot+' / 'rot-'          rotate 15° clockwise / counter-clockwise"""
        if fn_name in ('rot+', 'rot-'):
            self.rotate(15 if fn_name == 'rot+' else -15)
            return
        items = self._targets()
        c = self._center(items)
        if c is None:
            return
        cx, cy = c
        self.snapshot()
        fn = (lambda x, y: (2 * cx - x, y)) if fn_name.endswith('_h') else (lambda x, y: (x, 2 * cy - y))
        self.doc.transform(items, fn)
        for kind, i in items:
            if kind == 'bond' and self.doc.bonds[i].stereo and fn_name.startswith('flip'):
                b = self.doc.bonds[i]
                b.stereo = 'hash' if b.stereo == 'wedge' else 'wedge'
            elif kind == 'arrow' and self.doc.arrows[i].bend:
                self.doc.arrows[i].bend *= -1
        self.edited()

    def rotate(self, degrees):
        items = self._targets()
        c = self._center(items)
        if c is None:
            return
        cx, cy = c
        a = math.radians(degrees)
        co, si = math.cos(a), math.sin(a)
        self.snapshot()
        self.doc.transform(items, lambda x, y: (cx + (x - cx) * co - (y - cy) * si, cy + (x - cx) * si + (y - cy) * co))
        self.edited()

    def scale_selection(self, factor):
        items = self._targets()
        c = self._center(items)
        if c is None or factor <= 0:
            return
        cx, cy = c
        self.snapshot()
        self.doc.transform(items, lambda x, y: (cx + (x - cx) * factor, cy + (y - cy) * factor))
        self.edited()

    def object_groups(self):
        """The selected objects for Align / Distribute: each structure (its selected atoms and bonds), each arrow,
        text and shape - as [(items, QRectF in pt)]."""
        rend = self.renderer()
        atoms = self.selected_atoms()
        groups = []
        for comp in self.doc.fragments(atoms):
            cs = set(comp)
            items = [('atom', a) for a in cs] + [('bond', i) for i, b in self.doc.bonds.items() if b.a in cs and b.b in cs]
            rect = QRectF()
            for it in items:
                r = rend.item_rect(*it)
                rect = r if rect.isNull() else rect.united(r)
            groups.append((items, rect))
        for kind in ('arrow', 'text', 'shape'):
            for i in self.selected(kind):
                groups.append(([(kind, i)], rend.item_rect(kind, i)))
        return groups

    def align(self, how):
        """'left' | 'hcenter' | 'right' | 'top' | 'vcenter' | 'bottom' the selected objects to each other."""
        groups = self.object_groups()
        if len(groups) < 2:
            self.message.emit('Select two or more objects (structures, arrows, text, shapes) to align.')
            return
        BL = self.renderer().BL
        rects = [r for _i, r in groups]
        target = {'left': min(r.left() for r in rects), 'right': max(r.right() for r in rects),
                  'top': min(r.top() for r in rects), 'bottom': max(r.bottom() for r in rects),
                  'hcenter': sum(r.center().x() for r in rects) / len(rects),
                  'vcenter': sum(r.center().y() for r in rects) / len(rects)}[how]
        self.snapshot()
        for items, r in groups:
            cur = {'left': r.left(), 'right': r.right(), 'top': r.top(), 'bottom': r.bottom(),
                   'hcenter': r.center().x(), 'vcenter': r.center().y()}[how]
            d = (target - cur) / BL
            dx, dy = (d, 0.0) if how in ('left', 'right', 'hcenter') else (0.0, d)
            self.doc.transform(items, lambda x, y: (x + dx, y + dy))
        self.edited()

    def distribute(self, axis):
        """Equal gaps between the selected objects, left to right ('h') or top to bottom ('v')."""
        groups = self.object_groups()
        if len(groups) < 3:
            self.message.emit('Select three or more objects to distribute.')
            return
        BL = self.renderer().BL
        key = (lambda g: g[1].center().x()) if axis == 'h' else (lambda g: g[1].center().y())
        groups.sort(key=key)
        size = (lambda r: r.width()) if axis == 'h' else (lambda r: r.height())
        start = groups[0][1].left() if axis == 'h' else groups[0][1].top()
        end = groups[-1][1].right() if axis == 'h' else groups[-1][1].bottom()
        gap = (end - start - sum(size(r) for _i, r in groups)) / (len(groups) - 1)
        self.snapshot()
        pos = start
        for items, r in groups:
            cur = r.left() if axis == 'h' else r.top()
            d = (pos - cur) / BL
            dx, dy = (d, 0.0) if axis == 'h' else (0.0, d)
            self.doc.transform(items, lambda x, y: (x + dx, y + dy))
            pos += size(r) + gap
        self.edited()

    def edit_selected(self, kind, **props):
        """Set properties on the selected atoms or bonds (the atoms / bonds of the selection), one undo step."""
        if kind == 'atom':
            targets = [self.doc.atoms[a] for a in self.selected_atoms()]
        else:
            atoms = self.selected_atoms()
            targets = [b for i, b in self.doc.bonds.items() if ('bond', i) in self.sel or (b.a in atoms and b.b in atoms)]
        if not targets:
            return
        self.snapshot()
        for t in targets:
            for k, v in props.items():
                setattr(t, k, v)
        self.edited()

    def reverse_bonds(self):
        """Swap the ends of the selected bonds (a wedge then points the other way)."""
        bonds = [self.doc.bonds[i] for k, i in self.sel if k == 'bond']
        if bonds:
            self.snapshot()
            for b in bonds:
                b.a, b.b = b.b, b.a
            self.edited()

    def invert_selection(self):
        self.set_selection(set(self.doc.all_items()) - self.sel)

    def set_highlight(self, color):
        """Highlight the selected atoms and bonds with a colour band behind them ('' removes it)."""
        atoms = self.selected_atoms()
        bonds = [i for i, b in self.doc.bonds.items() if ('bond', i) in self.sel or (b.a in atoms and b.b in atoms)]
        if not atoms and not bonds:
            self.message.emit('Select atoms or bonds to highlight first.')
            return
        self.snapshot()
        for a in atoms:
            self.doc.atoms[a].highlight = color
        for i in bonds:
            self.doc.bonds[i].highlight = color
        self.edited()

    def set_ring_fill(self, color):
        """Fill the rings of the selection (all its atoms selected) with a colour ('' removes the fill)."""
        atoms = self.selected_atoms()
        rings = self.doc.rings(atoms) if atoms else []
        if not rings:
            self.message.emit('Select a ring (or a structure with rings) to fill first.')
            return
        self.snapshot()
        for ring in rings:
            old = self.doc.fill_of_ring(ring)
            if not color:
                if old is not None:
                    del self.doc.fills[old]
            elif old is not None:
                self.doc.fills[old].color = color
            else:
                self.doc.fills[self.doc._id()] = model.Fill(list(ring), color)
        self.edited()

    def add_frame(self, kind):
        """A rectangle / rounded rectangle / brackets around the selection (Add Frame)."""
        rend = self.renderer()
        r = QRectF()
        for it in self.sel:
            if it[1] in self.doc.table(it[0]):
                ir = rend.item_rect(*it)
                r = ir if r.isNull() else r.united(ir)
        if r.isNull():
            self.message.emit('Select what to frame first.')
            return
        pad = render.style_of(self.doc)['font_size'] * 0.8
        r = r.adjusted(-pad, -pad, pad, pad)
        BL = rend.BL
        self.snapshot()
        i = self.doc.add_shape(r.left() / BL, r.top() / BL, r.right() / BL, r.bottom() / BL, kind)
        self.sel.add(('shape', i))
        self.edited()

    def set_color(self, color):
        """Colour the selection ('' = back to the drawing's colour)."""
        if not self.sel:
            self.message.emit('Select what to colour first.')
            return
        self.snapshot()
        atoms = self.selected_atoms()
        for a in atoms:
            self.doc.atoms[a].color = color
        for i, b in self.doc.bonds.items():
            if ('bond', i) in self.sel or (b.a in atoms and b.b in atoms):
                b.color = color
        for kind in ('arrow', 'text', 'shape'):
            for i in self.selected(kind):
                self.doc.table(kind)[i].color = color
        self.edited()

    def set_text_props(self, **props):
        """Change the selected captions (and what new captions get)."""
        self.text_props.update(props)
        texts = self.selected('text')
        if texts:
            self.snapshot()
            for i in texts:
                for k, v in props.items():
                    setattr(self.doc.texts[i], k, v)
            self.edited()

    def expand_labels(self):
        """Write the selected nicknames (all, when nothing is selected) out as atoms."""
        atoms = [a for a in (self.selected_atoms() or self.doc.atoms) if self.doc.atoms[a].label in model.NICKNAMES]
        if not atoms:
            self.message.emit('No nickname (OMe, Boc, Ph ...) to expand.')
            return
        self.snapshot()
        for a in atoms:
            self.doc.put_group(a, model.NICKNAMES[self.doc.atoms[a].label])
        self.edited()

    def add_analysis(self):
        """A caption under each selected structure (all, when nothing is selected): formula, exact mass, m/z,
        molecular weight and elemental analysis - like a property label."""
        import chem
        atoms = self.selected_atoms() or set(self.doc.atoms)
        frags = [c for c in self.doc.fragments() if atoms & set(c)]
        if not frags:
            self.message.emit('Draw a structure first.')
            return
        self.snapshot()
        for comp in frags:
            info = chem.analyse(self.doc, comp)
            peaks = chem.isotope_pattern(info.counts, 0)[:3]
            mz = ', '.join(f'{m:.2f} ({p:.1f}%)' for m, p in peaks)
            text = (f'Chemical Formula: {info.formula}{chem.charge_text(info.charge)}\nExact Mass: {info.exact:.2f}\n'
                    f'm/z: {mz}\nMolecular Weight: {info.average:.2f}\nElemental Analysis: {info.analysis}')
            self._caption_under(comp, text, True)
        self.edited()

    def add_caption_below(self, atoms, text, formula=False):
        """A centred caption under the structure(s) holding these atoms, or under a sequence when atoms is
        ('seq', id) (the analysis window's Paste)."""
        if isinstance(atoms, tuple):
            if atoms[1] not in self.doc.seqs:
                return
        else:
            atoms = [a for a in atoms if a in self.doc.atoms]
            if not atoms:
                return
        self.snapshot()
        self._caption_under(atoms, text, formula)
        self.edited()

    def _caption_under(self, atoms, text, formula):
        """Add text centred under atoms (or under the sequence ('seq', id)), below any caption or shape already
        under them."""
        rend = self.renderer()
        BL = rend.BL
        rect = QRectF()
        if isinstance(atoms, tuple):
            rect = rend.item_rect(*atoms)
            atoms = []
        for a in atoms:
            r = rend.item_rect('atom', a)
            rect = r if rect.isNull() else rect.united(r)
        below = rect.bottom()
        for kind in ('text', 'shape', 'arrow'):
            for i in self.doc.table(kind):
                r = rend.item_rect(kind, i)
                if r.right() > rect.left() and r.left() < rect.right() and rect.top() < r.top() < below + 3 * BL:
                    below = max(below, r.bottom())
        self.doc.add_text(rect.center().x() / BL, (below + 6) / BL, text, align='center', formula=formula)
        self._renderer = None

    def analyze_stoichiometry(self):
        """Structure > Analyze stoichiometry: a grid under the reaction (the selected arrow's, else the first
        straight arrow's), with a column for each reactant, reagent and product."""
        import stoich
        arrows = [i for k, i in self.sel if k == 'arrow' and not self.doc.arrows[i].kind.startswith('curved')]
        cols = stoich.new_grid_columns(self.doc, arrows[0] if arrows else None)
        if cols is None:
            self.message.emit('Draw a reaction first: structures left and right of a reaction arrow.')
            return None
        rend = self.renderer()
        BL = rend.BL
        b = rend.bounds()
        self.snapshot()
        gid = self.doc._id()
        self.doc.grids[gid] = model.Grid(b.left() / BL, (b.bottom() + 2.2 * render.style_of(self.doc)['font_size']) / BL,
                                         cols[0], cols[1])
        self.sel = {('grid', gid)}
        self.edited()
        self.ensure_visible()
        return gid

    def edit_grid(self, gid):
        from stoich_dialog import GridDialog
        GridDialog(self, gid, self.window()).exec()

    def check_structure(self):
        """Structure > Check structure: the problems in the selected structures (all, when nothing is selected)."""
        import chem
        atoms = self.selected_atoms() or set(self.doc.atoms)
        problems = []
        for comp in self.doc.fragments(atoms):
            info = chem.analyse(self.doc, comp)
            for e in info.errors:
                problems.append(f'{info.formula or "?"}: {e}')
        unknown = {self.doc.atoms[a].label for a in atoms if self.doc.atoms[a].label
                   and self.doc.atoms[a].label not in model.NICKNAMES}
        problems += [f'Label "{u}" is not a known nickname: it counts as {"C"}.' for u in sorted(unknown)]
        return problems

    def clean_up(self):
        """Lay out the selected structures (all, when nothing is selected) again with RDKit, keeping their place."""
        import chem
        import layout2d
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
            layout2d.compute(mol)
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
        z = min(max(z, MIN_ZOOM), 12.0)
        anchor = anchor if anchor is not None else QPointF(self.width() / 2, self.height() / 2)
        mx, my = self.to_model(anchor)
        self.zoom = z
        p = self.to_px(mx, my)
        self.origin += anchor - p
        self.update()
        self.zoom_changed.emit(z)

    def fit(self):
        """Show the whole drawing (at most 3x); an empty drawing: the top of the page at 150 %."""
        b = render.Renderer(self.doc).bounds()
        if b.isNull():
            self.zoom = 1.5
            W, _H = self.page_size()
            k = self.px_per_pt()
            self.origin = QPointF(max(20.0, (self.width() - W * k) / 2), 20.0)
        else:
            k = self.logicalDpiX() / 72.0
            z = min((self.width() - 60) / max(b.width() * k, 1), (self.height() - 60) / max(b.height() * k, 1), 3.0)
            self.zoom = max(z, MIN_ZOOM)
            c = b.center()
            self.origin = QPointF(self.width() / 2 - c.x() * self.px_per_pt(), self.height() / 2 - c.y() * self.px_per_pt())
        self.update()
        self.zoom_changed.emit(self.zoom)

    def show_page(self):
        """The whole page in the window."""
        W, H = self.page_size()
        k = self.logicalDpiX() / 72.0
        self.zoom = max(MIN_ZOOM, min((self.width() - 40) / (W * k), (self.height() - 40) / (H * k)))
        kk = self.px_per_pt()
        self.origin = QPointF((self.width() - W * kk) / 2, (self.height() - H * kk) / 2)
        self.update()
        self.zoom_changed.emit(self.zoom)

    def ensure_visible(self):
        b = render.Renderer(self.doc).bounds()
        if not b.isNull() and not QRectF(self.rect()).contains(self.rect_px(b)):
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
        pt = QPointF(x * BL, y * BL)
        for i in self.doc.texts:
            if rend.text_rect(i).adjusted(-1, -1, 1, 1).contains(pt):
                return ('text', i)
        for a2 in self.doc.atoms:                           # a label wider than the atom circle (NH2, OMe)
            r = rend.label_rect(a2)
            if r is not None and r.adjusted(-1, -1, 1, 1).contains(pt):
                return ('atom', a2)
        b = self.doc.find_bond(x, y, max(0.14, 5 / s))
        if b is not None:
            return ('bond', b)
        ar = self.doc.find_arrow(x, y, max(0.15, 6 / s))
        if ar is not None:
            return ('arrow', ar)
        sh = self.doc.find_shape(x, y, max(0.12, 5 / s))
        if sh is not None:
            return ('shape', sh)
        for i in self.doc.grids:
            if rend.grid_rect(i).contains(pt):
                return ('grid', i)
        for i in self.doc.seqs:
            if rend.item_rect('seq', i).adjusted(-2, -2, 2, 2).contains(pt):
                return ('seq', i)
        for i in self.doc.tlcs:
            if rend.item_rect('tlc', i).adjusted(-3, -3, 3, 3).contains(pt):
                return ('tlc', i)
        return None

    # ------------------------------------------------------------------------------------------- painting
    def paintEvent(self, _e):
        p = QPainter(self)
        p.fillRect(self.rect(), PAGE_BG)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        W, H = self.page_size()
        page = self.rect_px(QRectF(0, 0, W, H))
        p.fillRect(page.translated(3, 3), PAGE_SHADOW)
        p.fillRect(page, QColor('#FFFFFF'))
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
        if self.show_reaction:
            self._paint_reaction(p, rend, BL)
        rend.paint(p, errors=self.show_warnings)
        if self.sel:                                        # the selection's frame, as a dashed box
            r = self._sel_box(rend)
            if not r.isNull():
                p.setPen(QPen(accent, 0.6 / max(self.zoom, 0.3), Qt.PenStyle.DashLine)); p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(r)
        if self.preview:
            self._paint_preview(p, BL, accent)
        if self.press and self.press.get('lasso'):
            p.resetTransform()
            p.setPen(QPen(accent, 1, Qt.PenStyle.DashLine))
            fill = QColor(accent); fill.setAlpha(25)
            p.setBrush(fill)
            from PyQt6.QtGui import QPolygonF
            p.drawPolygon(QPolygonF(self.press['lasso']))
        elif self.press and self.press.get('band'):
            p.resetTransform()
            p.setPen(QPen(accent, 1, Qt.PenStyle.DashLine))
            fill = QColor(accent); fill.setAlpha(25)
            p.setBrush(fill)
            p.drawRect(QRectF(self.press['pos'], self.mouse).normalized())
        if self.show_rulers:
            self._paint_rulers(p)
        p.end()

    def _paint_rulers(self, p):
        """Centimetre rulers along the top and left edges, measured from the page's top-left corner."""
        p.resetTransform()
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        size = 18
        bg, fg = QColor('#F4F5F7'), QColor('#5A6270')
        p.fillRect(QRectF(0, 0, self.width(), size), bg)
        p.fillRect(QRectF(0, 0, size, self.height()), bg)
        p.setPen(QPen(QColor('#B8BEC6'), 1))
        p.drawLine(QPointF(0, size), QPointF(self.width(), size))
        p.drawLine(QPointF(size, 0), QPointF(size, self.height()))
        cm = 72 / 2.54 * self.px_per_pt()                   # pixels per centimetre
        step = 1 if cm > 25 else 2 if cm > 12 else 5
        f = QFont(self.font()); f.setPixelSize(9); p.setFont(f)
        for axis in ('x', 'y'):
            o = self.origin.x() if axis == 'x' else self.origin.y()
            length = self.width() if axis == 'x' else self.height()
            k0 = int((size - o) / cm) - 1
            k1 = int((length - o) / cm) + 1
            for k in range(k0, k1 + 1):
                for sub in range(10 if cm > 60 else 2 if cm > 20 else 1):
                    pos = o + (k + sub / (10 if cm > 60 else 2 if cm > 20 else 1)) * cm
                    if pos < size or pos > length:
                        continue
                    major = sub == 0
                    tick = 8 if major else 4
                    p.setPen(QPen(fg, 1))
                    if axis == 'x':
                        p.drawLine(QPointF(pos, size - tick), QPointF(pos, size))
                        if major and k % step == 0:
                            p.drawText(QPointF(pos + 2, 10), str(k))
                    else:
                        p.drawLine(QPointF(size - tick, pos), QPointF(size, pos))
                        if major and k % step == 0:
                            p.save(); p.translate(10, pos + 2); p.rotate(-90); p.drawText(QPointF(-12, 0), str(k)); p.restore()
        p.fillRect(QRectF(0, 0, size, size), bg)

    def _mark(self, p, rend, kind, i, color, BL):
        if i not in self.doc.table(kind):
            return
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QBrush(color))
        if kind == 'atom':
            at = self.doc.atoms[i]
            r = rend.label_rect(i)
            if r is not None and rend.labels[i]['clip']:
                p.drawRoundedRect(r.adjusted(-1.5, -1.5, 1.5, 1.5), 2, 2)
            else:
                p.drawEllipse(QPointF(at.x * BL, at.y * BL), 0.2 * BL, 0.2 * BL)
        elif kind == 'bond':
            b = self.doc.bonds[i]
            pa, pb = self.doc.atoms[b.a], self.doc.atoms[b.b]
            p.setPen(QPen(color, 0.32 * BL, cap=Qt.PenCapStyle.RoundCap)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawLine(QPointF(pa.x * BL, pa.y * BL), QPointF(pb.x * BL, pb.y * BL))
        elif kind == 'arrow':
            path = QPainterPath()
            pts = model.arrow_points(self.doc.arrows[i])
            path.moveTo(pts[0][0] * BL, pts[0][1] * BL)
            for x, y in pts[1:]:
                path.lineTo(x * BL, y * BL)
            p.setPen(QPen(color, 0.3 * BL, cap=Qt.PenCapStyle.RoundCap)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(path)
        elif kind == 'shape':
            p.setPen(QPen(color, 0.25 * BL)); p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(rend.shape_paths(i))
        elif kind == 'text':
            p.drawRoundedRect(rend.text_rect(i).adjusted(-2, -1, 2, 1), 2, 2)
        elif kind == 'grid':
            p.drawRoundedRect(rend.grid_rect(i).adjusted(-2, -2, 2, 2), 3, 3)
        elif kind in ('seq', 'tlc'):                        # a light tint and a frame, so the content stays readable
            light = QColor(color); light.setAlpha(min(color.alpha(), 28))
            p.setBrush(QBrush(light)); p.setPen(QPen(color, 1.2))
            p.drawRoundedRect(rend.item_rect(kind, i).adjusted(-3, -3, 3, 3), 4, 4)

    def _paint_reaction(self, p, rend, BL):
        """View > Show reaction interpretation: reactants red, products magenta, above the arrow cyan, below blue
        (the colours ChemDraw uses), as soft halos behind the structures."""
        import stoich
        found = stoich.reaction_columns(self.doc)
        if not found:
            return
        reac, prod = found
        ar = next(a for a in self.doc.arrows.values() if not a.kind.startswith('curved'))
        lo, hi = sorted((ar.x1, ar.x2))
        ay = (ar.y1 + ar.y2) / 2
        for frags, role in ((reac, 'reactant'), (prod, 'product')):
            for frag in frags:
                cx = sum(self.doc.atoms[a].x for a in frag) / len(frag)
                cy = sum(self.doc.atoms[a].y for a in frag) / len(frag)
                col = QColor({'reactant': '#E53935', 'product': '#D81B9A', 'above': '#00ACC1', 'below': '#1E63E9'}[
                    role if role == 'product' or not lo <= cx <= hi else ('above' if cy < ay else 'below')])
                col.setAlpha(70)
                p.setPen(QPen(col, 0.38 * BL, cap=Qt.PenCapStyle.RoundCap)); p.setBrush(Qt.BrushStyle.NoBrush)
                fs = set(frag)
                for b in self.doc.bonds.values():
                    if b.a in fs:
                        pa, pb = self.doc.atoms[b.a], self.doc.atoms[b.b]
                        p.drawLine(QPointF(pa.x * BL, pa.y * BL), QPointF(pb.x * BL, pb.y * BL))
                for a in frag:
                    at = self.doc.atoms[a]
                    p.drawPoint(QPointF(at.x * BL, at.y * BL))

    def _paint_preview(self, p, BL, accent):
        kind, pts = self.preview
        pen = QPen(accent, 1.2, cap=Qt.PenCapStyle.RoundCap, join=Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)
        path = QPainterPath()
        if kind in SHAPE_TOOLS:
            sk, st, params = SHAPE_TOOLS[kind]
            tmp = Doc()
            tmp.style = self.doc.style
            j = tmp.add_shape(*pts[0], *pts[1], sk)
            tmp.shapes[j].style, tmp.shapes[j].params = st, dict(params)
            if sk in model.BIO_KINDS or sk in model.ELEMENT_KINDS or st:   # the object itself, half transparent
                p.save(); p.setOpacity(0.6)
                render.Renderer(tmp).paint(p)
                p.restore()
                return
            p.drawPath(render.Renderer(tmp).shape_paths(j))
            return
        if kind in CUT_TOOLS:                               # the cut line
            pen.setStyle(Qt.PenStyle.DashLine); p.setPen(pen)
            p.drawLine(QPointF(pts[0][0] * BL, pts[0][1] * BL), QPointF(pts[1][0] * BL, pts[1][1] * BL))
            return
        if kind == 'tlc':
            p.drawRect(QRectF(QPointF(pts[0][0] * BL, pts[0][1] * BL), QPointF(pts[1][0] * BL, pts[1][1] * BL)))
            return
        if kind in ('curved', 'curved_half'):
            ar = model.Arrow(*pts[0], *pts[1], kind, CURVE)
            pts = model.arrow_points(ar)
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

    # ------------------------------------------------------------------------------------------- typing in place
    def _editor_font(self, size_pt, bold=False, italic=False):
        f = QFont(render.style_of(self.doc)['font'])
        f.setPixelSize(max(9, round(size_pt * self.px_per_pt())))
        f.setBold(bold); f.setItalic(italic)
        return f

    def edit_label(self, a):
        """Type an atom's label in place (an element, a charge, OH / NH2, or a nickname such as OMe, Boc, CO2H)."""
        self.commit_editor()
        at = self.doc.atoms[a]
        w = _LabelEdit(self)
        w.setFont(self._editor_font(render.style_of(self.doc)['font_size']))
        w.setText(at.label or (at.el + (render.Renderer._charge_text(at.charge).replace('−', '-') if at.charge else '')))
        w.selectAll()
        width = max(70, w.fontMetrics().horizontalAdvance('MMMMM'))
        c = self.to_px(at.x, at.y)
        w.setGeometry(round(c.x() - width / 2), round(c.y() - w.sizeHint().height() / 2), width, w.sizeHint().height())
        self.editor = (w, 'label', a, at.x, at.y)
        w.show(); w.setFocus()

    def edit_text(self, i=None, x=0.0, y=0.0):
        """Type a caption in place: a new one at (x, y), or caption i. Esc or a click elsewhere ends typing."""
        self.commit_editor()
        props = dict(self.text_props) if i is None else {k: getattr(self.doc.texts[i], k) for k in self.text_props}
        size = props['size'] or render.style_of(self.doc)['caption_size']
        w = _TextEdit(self)
        w.setFont(self._editor_font(size, props['bold'], props['italic']))
        if i is not None:
            t = self.doc.texts[i]
            w.setPlainText(t.text)
            r = self.rect_px(self.renderer().text_rect(i))
            pos = r.topLeft()
        else:
            pos = self.to_px(x, y)
        w.fit()
        w.move(round(pos.x()) - 6, round(pos.y()) - 6)
        self.editor = (w, 'text', i, x, y)
        w.show(); w.setFocus()
        w.selectAll() if i is not None else None

    def commit_editor(self, cancel=False):
        if not self.editor:
            return
        w, kind, i, x, y = self.editor
        self.editor = None
        text = w.text() if kind == 'label' else w.toPlainText()
        w.hide(); w.deleteLater()
        self.setFocus()
        if cancel:
            return
        if kind == 'label':
            if i not in self.doc.atoms or not text.strip():
                return
            what = model.parse_label(text)
            if what is None:
                self.message.emit(f'"{text}" is not an element or a known nickname (OMe, Boc, CO2H, Ph ...).')
                return
            self.snapshot()
            at = self.doc.atoms[i]
            if what[0] == 'nickname':
                self.doc.set_nickname(i, what[1])
            else:
                at.label, at.el, at.charge = '', what[1], what[2]
            self.edited()
            return
        if i is None:
            if text.strip():
                self.snapshot()
                self.doc.add_text(x, y, text, **self.text_props)
                self.edited()
        elif i in self.doc.texts and text != self.doc.texts[i].text:
            self.snapshot()
            if text.strip():
                self.doc.texts[i].text = text
            else:
                del self.doc.texts[i]
            self.edited()

    # ------------------------------------------------------------------------------------------- mouse
    def mousePressEvent(self, e):
        if self.editor:                                     # a click elsewhere ends typing (and does nothing else)
            self.commit_editor()
            return
        self.setFocus()
        if e.button() in (Qt.MouseButton.MiddleButton, Qt.MouseButton.RightButton):
            self.press = dict(pan=True, pos=e.position(), origin=QPointF(self.origin), moved=False,
                              right=e.button() == Qt.MouseButton.RightButton, item=self.item_at(e.position()),
                              global_pos=e.globalPosition().toPoint())
            return
        if e.button() != Qt.MouseButton.LeftButton:
            return
        pos = e.position()
        item = self.item_at(pos)
        self.press = dict(pos=pos, item=item, moved=False, shift=bool(e.modifiers() & Qt.KeyboardModifier.ShiftModifier))
        if self.tool in ('select', 'lasso'):
            grabbed = self._charge_at(pos)
            if grabbed is not None:                         # the charge sign itself, not its atom
                self.press['charge'] = grabbed
                return
            if item and item[0] == 'tlc':                   # a spot on a plate: drag it (Rf), resize, copy
                import tlc
                rend = self.renderer()
                x, y = self.to_model(pos)
                hit = tlc.spot_at(rend, item[1], x * rend.BL, y * rend.BL)
                if hit is not None:
                    k, j = hit
                    s = self.doc.tlcs[item[1]].lanes[k][j]
                    self.press['spot'] = (item[1], k, j, dict(s))
                    self.press['ctrl'] = bool(e.modifiers() & Qt.KeyboardModifier.ControlModifier)
                    self.sel = {item}
                    self.selection_changed.emit()
                    return
        if self.tool in ('select', 'lasso'):
            if item and item in self.sel:
                self.press['move'] = list(self.sel)
            elif not item and not self.press['shift'] and self._in_sel_box(pos):
                self.press['move'] = list(self.sel)         # anywhere inside the dashed box moves the selection
            elif item:
                if not self.press['shift']:
                    self.sel = set()
                self.sel.add(item)
                self.press['move'] = list(self.sel)
                self.selection_changed.emit()
            else:
                self.press['band'] = True
                if self.tool == 'lasso':
                    self.press['lasso'] = [pos]
                if not self.press['shift']:
                    self.sel = set()
                    self.selection_changed.emit()
            if 'move' in self.press:
                self.press['start'] = self.to_model(pos)
            self.update()
        elif self.tool == 'erase':
            self.snapshot()
            self.press['erased'] = False
            self._erase(item)

    def _sel_box(self, rend):
        """The dashed frame around the selection, in pt (null when nothing is selected)."""
        r = QRectF()
        for it in self.sel:
            if it[1] in self.doc.table(it[0]):
                ir = rend.item_rect(*it)
                r = ir if r.isNull() else r.united(ir)
        return r if r.isNull() else r.adjusted(-3, -3, 3, 3)

    def _in_sel_box(self, pos):
        rend = self.renderer()
        x, y = self.to_model(pos)
        return self._sel_box(rend).contains(QPointF(x * rend.BL, y * rend.BL))

    def _charge_at(self, pos):
        """The atom whose charge sign is under the mouse, or None."""
        x, y = self.to_model(pos)
        rend = self.renderer()
        pt = QPointF(x * rend.BL, y * rend.BL)
        for a, at in self.doc.atoms.items():
            if at.charge:
                r = rend.charge_rect(a)
                if r is not None and r.adjusted(-1.2, -1.2, 1.2, 1.2).contains(pt):
                    return a
        return None

    def _explain(self, item, e):
        """Pointing at an atom with a bond-count problem says what is wrong (tooltip and status bar)."""
        from PyQt6.QtWidgets import QToolTip
        problem = self.doc.valence_problem(item[1]) if item and item[0] == 'atom' and self.show_warnings else ''
        if item and item[0] == 'seq':                       # pointing at a residue: what it is and its number
            import biodraw
            import biopolymer
            rend = self.renderer()
            x, y = self.to_model(e.position())
            k = biodraw.residue_at(rend, item[1], x * rend.BL, y * rend.BL)
            if k is not None:
                q = self.doc.seqs[item[1]]
                problem = f'{k + 1}: {biopolymer.residue_name(q.polymer, q.residues[k])}'
        elif item and item[0] == 'shape' and self.doc.shapes[item[1]].kind in model.BIO_KINDS:
            import biodraw
            problem = biodraw.NAMES.get(self.doc.shapes[item[1]].kind, '')
        if problem:
            QToolTip.showText(e.globalPosition().toPoint(), problem, self)
            self.message.emit(problem)
        else:
            QToolTip.hideText()

    def reset_charge_positions(self):
        atoms = [a for a in self.selected_atoms() if self.doc.atoms[a].charge_pos]
        if atoms:
            self.snapshot()
            for a in atoms:
                self.doc.atoms[a].charge_pos = None
            self.edited()

    def _erase(self, item):
        if item and item[0] == 'tlc' and item[1] in self.doc.tlcs:   # a spot (its lane goes with its last spot)
            import tlc
            rend = self.renderer()
            x, y = self.to_model(self.press['pos'])
            hit = tlc.spot_at(rend, item[1], x * rend.BL, y * rend.BL)
            if hit is not None:
                plate = self.doc.tlcs[item[1]]
                del plate.lanes[hit[0]][hit[1]]
                if not plate.lanes[hit[0]]:
                    del plate.lanes[hit[0]]
                if not plate.lanes:
                    self.doc.tlcs.pop(item[1])
                self.press['erased'] = True
                self.edited()
                return
        if item and item[1] in self.doc.table(item[0]):
            if item[0] == 'atom' and self.doc.strip_atom(item[1]):  # N, OMe, C+ ... become a plain carbon first
                self.press['erased'] = True
                self.edited()
                return
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
            if not pr['moved'] and (pos - pr['pos']).manhattanLength() < DRAG_PX:
                return
            if not pr['moved']:
                pr['moved'] = True
                self.setCursor(Qt.CursorShape.ClosedHandCursor)
            self.origin = pr['origin'] + (pos - pr['pos'])
            self.update()
            return
        if pr is None:
            item = self.item_at(pos)
            if item != self.hover:
                self.hover = item
                self.update()
                self._explain(item, e)
            elif item and item[0] == 'seq':                 # residue by residue along a sequence
                self._explain(item, e)
            if self.tool in ('select', 'lasso'):            # a charge sign can be dragged
                self.setCursor(Qt.CursorShape.SizeAllCursor if self._charge_at(pos) is not None
                               else Qt.CursorShape.ArrowCursor)
            return
        if not pr['moved'] and (pos - pr['pos']).manhattanLength() < DRAG_PX:
            return
        pr['moved'] = True
        x, y = self.to_model(pos)
        t = self.tool
        if 'spot' in pr:                                    # a TLC spot: Rf, size (Shift) or a copy (Ctrl)
            import tlc
            i, k, j, orig = pr['spot']
            if not pr.get('snapped'):
                pr['snapped'] = True
                self.snapshot()
            rend = self.renderer()
            spot = self.doc.tlcs[i].lanes[k][j]
            sx, sy = self.to_model(pr['pos'])
            if pr['shift']:
                spot['w'] = round(max(0.3, orig.get('w', 1.0) + (x - sx) * 3), 2)
                spot['h'] = round(max(0.3, orig.get('h', 1.0) + (y - sy) * 3), 2)
                self.message.emit(f"Spot size {spot['w']:.2f} × {spot['h']:.2f}")
            elif pr.get('ctrl'):
                self.message.emit(f'Copy to lane {tlc.lane_at(rend, i, x * rend.BL) + 1}, '
                                  f'Rf {tlc.rf_at(rend, i, y * rend.BL):.2f}')
            else:
                spot['rf'] = round(tlc.rf_at(rend, i, y * rend.BL), 3)
                self.message.emit(f"Rf {spot['rf']:.2f}")
            self._renderer = None
            self.update()
            return
        if 'charge' in pr:                                  # dragging a charge sign: it stays with its atom
            if not pr.get('snapped'):
                pr['snapped'] = True
                self.snapshot()
            at = self.doc.atoms[pr['charge']]
            at.charge_pos = [round(x - at.x, 3), round(y - at.y, 3)]
            self._renderer = None
            self.update()
            return
        if t in ('select', 'lasso'):
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
            elif 'lasso' in pr:
                pr['lasso'].append(pos)
                self._lasso_select(pr['lasso'], pr['shift'])
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
        elif t in ARROW_TOOLS:
            sx, sy = self.to_model(pr['pos'])
            ang = model.snap_angle(model.angle_of(x - sx, y - sy), 15)
            L = math.hypot(x - sx, y - sy)
            end = (x, y) if ARROW_TOOLS[t].startswith('curved') else (sx + L * math.cos(ang), sy + L * math.sin(ang))
            self.preview = (ARROW_TOOLS[t], [(sx, sy), end])
            self.update()
        elif t in SHAPE_TOOLS:
            start = self.to_model(pr['pos'])
            if t.startswith('orb:') and pr['item'] and pr['item'][0] == 'atom':   # orbitals start on their atom
                at = self.doc.atoms[pr['item'][1]]
                start = (at.x, at.y)
            self.preview = (t, [start, (x, y)])
            self.update()
        elif t in CUT_TOOLS or t == 'tlc':
            self.preview = (t, [self.to_model(pr['pos']), (x, y)])
            self.update()

    def _cut(self, tool, pts, keep_line=False):
        """Mass fragmentation / retrosynthesis / synthesis along the dragged line; False when it crosses no bond."""
        import fragment
        (x1, y1), (x2, y2) = pts
        line = (x1, y1, x2, y2)
        if tool == 'massfrag':
            frags = fragment.mass_fragments(self.doc, line)
            if not frags:
                self.message.emit('The line crosses no bond.')
                return False
            for formula, mass, (tx, ty) in frags:
                self.doc.add_text(tx, ty, f'{formula}\n{mass:.4f}', formula=True, align='center')
            if keep_line:
                j = self.doc.add_shape(x1, y1, x2, y2, 'line')
                self.doc.shapes[j].style = 'dashed'
            return True
        new = fragment.disconnect(self.doc, line, retro=(tool == 'retrosyn'))
        if not new:
            self.message.emit('The line crosses no bond.')
            return False
        self.sel = set(new)
        return True

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

    def _lasso_select(self, points, add):
        """Select what lies inside a free-hand outline (widget pixels): atoms by their position, bonds between two such
        atoms, arrows and shapes by both ends, captions and grids by their middle."""
        from PyQt6.QtGui import QPolygonF
        poly = QPolygonF(points)
        inside = lambda pt: poly.containsPoint(pt, Qt.FillRule.OddEvenFill)
        items = set()
        for i, at in self.doc.atoms.items():
            if inside(self.to_px(at.x, at.y)):
                items.add(('atom', i))
        for i, b in self.doc.bonds.items():
            if ('atom', b.a) in items and ('atom', b.b) in items:
                items.add(('bond', i))
        for kind, table in (('arrow', self.doc.arrows), ('shape', self.doc.shapes)):
            for i, a in table.items():
                if inside(self.to_px(a.x1, a.y1)) and inside(self.to_px(a.x2, a.y2)):
                    items.add((kind, i))
        rend = self.renderer()
        for i in self.doc.texts:
            if inside(self.rect_px(rend.text_rect(i)).center()):
                items.add(('text', i))
        for i in self.doc.grids:
            if inside(self.rect_px(rend.grid_rect(i)).center()):
                items.add(('grid', i))
        self.sel = (self.sel | items) if add else items
        self.update()
        self.selection_changed.emit()

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
        for i, s in self.doc.shapes.items():
            if r.contains(self.to_px(s.x1, s.y1)) and r.contains(self.to_px(s.x2, s.y2)):
                items.add(('shape', i))
        for i in self.doc.texts:
            if r.intersects(self.rect_px(self.renderer().text_rect(i))):
                items.add(('text', i))
        for i in self.doc.grids:
            if r.contains(self.rect_px(self.renderer().grid_rect(i))):
                items.add(('grid', i))
        self.sel = (self.sel | items) if add else items
        self.update()
        self.selection_changed.emit()

    def mouseReleaseEvent(self, e):
        pr, self.press = self.press, None
        if pr is None:
            return
        if pr.get('pan'):
            self.set_tool(self.tool)
            if pr.get('right') and not pr['moved']:        # a right click (not a drag): the object menu
                item = pr['item']
                if item and item not in self.sel:
                    self.set_selection([item])
                self.context_menu.emit(pr['global_pos'], item)
            return
        pos = e.position()
        x, y = self.to_model(pos)
        t, item = self.tool, pr['item']
        preview, self.preview = self.preview, None
        if 'charge' in pr:                                  # a charge sign was dragged (or just clicked)
            if pr.get('snapped'):
                self.edited()
            return
        if 'spot' in pr:
            if pr.get('snapped'):
                if pr.get('ctrl'):                          # Ctrl+drag: a copy in the lane under the mouse
                    import tlc
                    i, k, j, orig = pr['spot']
                    rend = self.renderer()
                    lane = tlc.lane_at(rend, i, x * rend.BL)
                    self.doc.tlcs[i].lanes[lane].append(dict(orig, rf=round(tlc.rf_at(rend, i, y * rend.BL), 3)))
                self.edited()
            return
        if t in ('select', 'lasso'):
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
                order, stereo, style = BOND_TOOLS[t]
                bid = self.doc.add_bond(a, b, order, stereo)
                if bid is not None:
                    bd = self.doc.bonds[bid]
                    if stereo and bd.a != a:
                        bd.a, bd.b = bd.b, bd.a
                    bd.style = style
            elif kind == 'chain':
                prev = item[1] if item and item[0] == 'atom' else self.doc.atom_at_or_new(*pts[0])
                for px, py in pts[1:]:
                    nxt = self.doc.atom_at_or_new(px, py, exclude={prev})
                    self.doc.add_bond(prev, nxt)
                    prev = nxt
            elif kind in model.ARROW_KINDS:
                (x1, y1), (x2, y2) = pts
                if math.hypot(x2 - x1, y2 - y1) > 0.3:
                    self.doc.add_arrow(x1, y1, x2, y2, kind, CURVE if kind.startswith('curved') else 0.0)
            elif kind in SHAPE_TOOLS:
                (x1, y1), (x2, y2) = pts
                if math.hypot(x2 - x1, y2 - y1) > 0.2:
                    sk, st, params = SHAPE_TOOLS[kind]
                    j = self.doc.add_shape(x1, y1, x2, y2, sk)
                    self.doc.shapes[j].style, self.doc.shapes[j].params = st, dict(params)
            elif kind in CUT_TOOLS:
                if not self._cut(kind, pts, bool(e.modifiers() & Qt.KeyboardModifier.AltModifier)):
                    self.undo_stack.pop()
                    return
            elif kind == 'tlc':
                import tlc
                (x1, y1), (x2, y2) = pts
                i = tlc.new_plate(self.doc, x1, y1, x2 - x1, y2 - y1)
                self.sel = {('tlc', i)}
            self.edited()
            return
        if not pr['moved']:
            self.click(item, x, y)

    def click(self, item, x, y):
        """A click (no drag) with the current tool on item (None: empty space at model x, y)."""
        t, doc = self.tool, self.doc
        if t in ARROW_TOOLS:
            kind = ARROW_TOOLS[t]
            if item and item[0] == 'arrow' and kind.startswith('curved') and doc.arrows[item[1]].kind.startswith('curved'):
                self.snapshot(); doc.arrows[item[1]].bend *= -1; self.edited()
                return
            self.snapshot()
            doc.add_arrow(x, y, x + 2.5, y, kind, CURVE if kind.startswith('curved') else 0.0)
            self.edited()
            return
        if t == 'sequence':                                 # the window asks for the sequence, then puts it here
            if item and item[0] == 'seq':
                self.edit_object.emit('seq', item[1])
            else:
                self.pending_pos = (x, y)
                self.edit_object.emit('seq', -1)
            return
        if t in model.BIO_KINDS:                            # a click: the object at its default size
            import biodraw
            self.snapshot()
            i = doc.add_shape(*biodraw.default_shape(t, x, y), t)
            if t == 'plasmid':
                doc.shapes[i].params = {'bp': 5000, 'name': '', 'regions': [], 'markers': []}
            self.sel = {('shape', i)}
            self.edited()
            if t == 'plasmid':
                self.edit_object.emit('shape', i)
            return
        if t.startswith('orb:'):                            # an orbital on the atom clicked (or where clicked)
            sk, st, params = SHAPE_TOOLS[t]
            cx, cy = (doc.atoms[item[1]].x, doc.atoms[item[1]].y) if item and item[0] == 'atom' else (x, y)
            size = 0.5 if params['type'] in ('s', 'sigma') else 0.85
            self.snapshot()
            j = doc.add_shape(cx, cy, cx, cy - size, 'orbital')
            doc.shapes[j].style, doc.shapes[j].params = st, dict(params)
            self.edited()
            return
        if t == 'tlc':
            import tlc
            self.snapshot()
            i = tlc.new_plate(doc, x - 1.5, y - 2.0, 3.0, 4.0)
            self.sel = {('tlc', i)}
            self.edited()
            return
        if t in SHAPE_TOOLS or t in CUT_TOOLS:
            return
        if t == 'text':
            if item and item[0] == 'atom':
                self.edit_label(item[1])
            else:
                self.edit_text(item[1] if item and item[0] == 'text' else None, x, y)
            return
        if item and item[0] in ('text', 'shape', 'arrow', 'grid'):
            return
        if t in BOND_TOOLS or t == 'chain':
            order, stereo, style = BOND_TOOLS.get(t, (1, '', ''))
            self.snapshot()
            if item is None:
                a = doc.add_atom(x - math.cos(math.radians(30)) / 2, y + 0.25)
                b = doc.atom_at_or_new(doc.atoms[a].x + math.cos(math.radians(30)), doc.atoms[a].y - 0.5, exclude={a})
                bid = doc.add_bond(a, b, order, stereo)
                doc.bonds[bid].style = style
            elif item[0] == 'atom':
                a = item[1]
                ang = doc.next_bond_angle(a)
                at = doc.atoms[a]
                b = doc.atom_at_or_new(at.x + math.cos(ang), at.y + math.sin(ang), exclude={a})
                bid = doc.add_bond(a, b, order, stereo)
                if bid is not None:
                    bd = doc.bonds[bid]
                    if stereo and bd.a != a:
                        bd.a, bd.b = bd.b, bd.a
                    bd.style = style
            else:
                bd = doc.bonds[item[1]]
                if stereo:
                    if bd.stereo == stereo:
                        bd.a, bd.b = bd.b, bd.a
                    bd.stereo, bd.order, bd.style = stereo, 1, ''
                elif style:
                    if style == 'dashed' and bd.order == 2:  # a dashed second line (aromatic / partial bond)
                        bd.style = '' if bd.style == 'dashed' else 'dashed'
                    else:
                        bd.style = '' if bd.style == style and bd.order == 1 else style
                        bd.order, bd.stereo = 1, ''
                elif t in ('bond1', 'chain'):
                    bd.stereo = bd.style = ''
                    bd.order = 1 if bd.order == 3 else bd.order + 1
                    bd.side = None
                elif t == 'bond2':
                    if bd.order == 2 and not bd.stereo:     # cycle the second line: auto -> other side -> centred -> auto
                        cur = doc.double_side(item[1])
                        if bd.side is None:
                            bd.side = -cur if cur else 1
                        else:
                            bd.side = 0 if bd.side else None
                    bd.order, bd.stereo, bd.style = 2, '', ''
                else:
                    bd.order, bd.stereo, bd.style = 3, '', ''
            self.edited()
        elif t in RING_TOOLS:
            n, aromatic = RING_TOOLS[t]
            self.snapshot()
            if item is None:
                pts = doc.ring_points(n, (x, y), math.radians(-90))
            elif item[0] == 'bond':
                pts = doc.ring_on_bond(item[1], n)
            else:
                pts = doc.ring_on_atom(item[1], n)
            doc.add_ring(pts, aromatic=aromatic)
            self.edited()
        elif t == 'atom':
            self.snapshot()
            if item is None:
                doc.add_atom(x, y, self.element)
            elif item[0] == 'atom':
                at = doc.atoms[item[1]]
                at.el, at.label = self.element, ''
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
        elif item and item[0] == 'grid':
            self.edit_grid(item[1])
        elif item and item[0] == 'seq':
            self.edit_object.emit('seq', item[1])
        elif item and item[0] == 'tlc':
            self.edit_object.emit('tlc', item[1])
        elif item and item[0] == 'shape' and self.doc.shapes[item[1]].kind == 'plasmid':
            self.edit_object.emit('shape', item[1])
        elif item and item[0] == 'bond' and self.tool == 'select':
            self.select_fragment(self.doc.bonds[item[1]].a)

    def wheelEvent(self, e):
        d = e.angleDelta()
        if e.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.set_zoom(self.zoom * (1.15 ** (d.y() / 120)), e.position())
        elif e.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.origin += QPointF(d.y() / 2, 0); self.update()
        else:
            self.origin += QPointF(d.x() / 2, d.y() / 2); self.update()
        if self.editor:
            self.commit_editor()

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
                self.snapshot()
                if not (h[0] == 'atom' and self.doc.strip_atom(h[1])):   # a label goes first, as in ChemDraw
                    self.doc.remove(*h)
                    self.hover = None
                self.edited()
            else:
                self.delete_selection()
            return
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and h and h[0] == 'atom':
            self.edit_label(h[1])                           # Enter over an atom: type its label (as in ChemDraw)
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
        if self.hotkey(e.text(), h):
            return
        super().keyPressEvent(e)

    def hotkey(self, ch, h):
        """A ChemDraw-style hotkey ch over the hovered item h (None: a tool key). True if it did something."""
        doc = self.doc
        if not ch:
            return False
        if h and h[0] == 'atom':
            a = h[1]
            if ch in ('+', '=', '-', '_'):
                self.snapshot(); doc.atoms[a].charge += 1 if ch in ('+', '=') else -1; self.edited()
                return True
            if ch in ATOM_KEYS:
                name = ATOM_KEYS[ch]
                self.snapshot()
                if name in model.NICKNAMES:
                    doc.set_nickname(a, name)
                else:
                    doc.atoms[a].el, doc.atoms[a].label = name, ''
                self.edited()
                return True
            if ch in SPROUT_KEYS:
                group, order, stereo = SPROUT_KEYS[ch]
                self.snapshot()
                doc.sprout(a, model.GROUPS.get(group, group), order, stereo)
                self.edited()
                return True
            if ch == '9':
                self.snapshot(); doc.sprout(a); doc.sprout(a); self.edited()
                return True
            return False
        if h and h[0] == 'bond':
            bd = doc.bonds[h[1]]
            if ch in ('1', '2', '3'):
                self.snapshot(); bd.order, bd.stereo, bd.style, bd.side = int(ch), '', '', None; self.edited()
                return True
            if ch in ('l', 'r', 'c') and bd.order == 2:     # second line left / right of the bond (on screen), centred
                self.snapshot(); bd.side = {'l': -1, 'r': 1, 'c': 0}[ch]; self.edited()
                return True
            if ch in FUSE_KEYS:
                n, aromatic = FUSE_KEYS[ch]
                self.snapshot(); doc.add_ring(doc.ring_on_bond(h[1], n), aromatic); self.edited()
                return True
            if ch in ('w', 'h', 'W'):
                stereo = 'wedge' if ch == 'w' else 'hash'
                self.snapshot()
                if bd.stereo == stereo:
                    bd.a, bd.b = bd.b, bd.a
                bd.stereo, bd.order, bd.style = stereo, 1, ''
                self.edited()
                return True
            if ch in BOND_STYLE_KEYS:
                self.snapshot()
                style = BOND_STYLE_KEYS[ch]
                if style == 'dashed' and bd.order == 2:
                    bd.style = 'dashed'
                else:
                    bd.style, bd.order, bd.stereo = style, 1, ''
                self.edited()
                return True
            return False
        if h is None and ch in TOOL_KEYS:
            self.tool_key.emit(TOOL_KEYS[ch])
            return True
        return False
