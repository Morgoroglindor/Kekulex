"""
kekulex.py - Kekulex: draw chemical structures for figures and slides and read off their formula,
exact mass and the m/z of the usual LC-MS adducts.

    python kekulex.py [drawing.sed | drawing.cdxml | file.cdx | file.mol | file.sdf | file.smi]

Left: the tool palette (drop-down groups for bonds, rings, arrows and shapes) and the elements. Top: the file / edit
bar and the text style bar. Middle: the page. Right: formula, masses, elemental analysis, adduct m/z with isotope
pattern, SMILES, InChIKey. Drawings are saved as .sed (JSON) or ChemDraw .cdxml; structures go out as MOL / SMILES,
pictures as PNG, SVG and PDF or to the clipboard (Ctrl+Shift+C) for PowerPoint / Word.
"""
import json
import os
import sys
import time

from PyQt6 import QtCore, QtGui, QtWidgets

import applog
import autosave
import biodraw
import biopolymer
import canvas as canvas_mod
import cdxml
import chem
import elements
import features
import icons
import model
import nmr
import render
import templates
import theme
import winkeys
from canvas import Canvas
from info_panel import InfoPanel
from model import Doc
from version import __version__

APP_NAME = 'Kekulex'
MIME = 'application/x-structure-editor'
OPEN_FILTER = ('Structures and drawings (*.sed *.cdxml *.cdx *.mol *.sdf *.sd *.smi *.smiles);;'
               'Kekulex drawings (*.sed);;ChemDraw files (*.cdxml *.cdx);;MOL / SDF files (*.mol *.sdf *.sd);;'
               'SMILES files (*.smi *.smiles *.txt);;All files (*)')
SAVE_FILTER = 'Kekulex drawing (*.sed);;ChemDraw XML (*.cdxml)'
_STYLE_NAMES = {'': '', 'dashed': 'dashed', 'shaded': 'shaded', 'filled': 'filled', 'shadow': 'with shadow',
                'wavy': 'wavy', 'bold': 'bold'}
_KIND_NAMES = {'circle': 'Circle', 'oval': 'Oval', 'roundrect': 'Rounded rectangle', 'rect': 'Rectangle', 'line': 'Line',
               'arc': 'Arc'}


def _element_tip(name):
    kind, style, params = canvas_mod.SHAPE_TOOLS[name]
    if kind == 'arc':
        return f"Arc {params['angle']}°" + (', dashed' if style else '') + ' (drag left to right: convex)'
    return (_KIND_NAMES[kind] + (', ' + _STYLE_NAMES[style] if style else '')
            + ('  (drag out from the centre)' if kind == 'circle' else ''))


ELEMENT_TIPS = [(n, _element_tip(n)) for n in canvas_mod.ELEMENT_TOOLS]
ORBITAL_TIPS = [(n, f"{elements.ORBITALS[n.split(':')[1]]}, {n.split(':')[2]}  (Shift+G)") for n in
                canvas_mod.ORBITAL_TOOLS]
# The palette: single tools, and drop-down groups (the button shows the last one used; click its arrow for the rest)
TOOLS = (('lasso', 'Lasso: select free-hand'), ('select', 'Select / move  (Space)'),
         ('erase', 'Eraser  (Delete over an item)'), None,
         ('bond1', 'Single bond  (x; 1 over a bond)'), ('bond2', 'Double bond  (2 over a bond)'),
         ('group', 'bonds', [('wedge', 'Wedged bond (toward you)  (w over a bond)'),
                             ('hash', 'Hashed wedged bond (away from you)  (h over a bond)'),
                             ('bold', 'Bold bond  (b over a bond)'), ('hashed', 'Hashed bond  (Shift+H over a bond)'),
                             ('dashed', 'Dashed bond / dashed second line  (d over a bond)'),
                             ('wavy', 'Wavy bond  (y over a bond)'), ('bond3', 'Triple bond  (3 over a bond)')]),
         ('chain', 'Chain  (Shift+X)'), None,
         ('benzene', 'Benzene  (j)'),
         ('group', 'rings', [('ring6', 'Cyclohexane'), ('ring5', 'Cyclopentane'), ('ring3', 'Cyclopropane'),
                             ('ring4', 'Cyclobutane'), ('ring7', 'Cycloheptane'), ('ring8', 'Cyclooctane'),
                             ('cyclopentadiene', 'Cyclopentadiene  (Shift+J)')]), None,
         ('group', 'charges', [('charge+', 'Positive charge  (+ over an atom)'),
                               ('charge-', 'Negative charge  (- over an atom)'),
                               ('radical', 'Radical •  (. over an atom)'),
                               ('diradical', 'Diradical / carbene ••  (. twice over an atom)'),
                               ('radical+', 'Radical cation •+'), ('radical-', 'Radical anion •−')], 2), None,
         ('group', 'arrows', [('arrow', 'Reaction arrow  (e)'), ('equilibrium', 'Equilibrium arrow'),
                              ('retro', 'Retrosynthetic arrow'), ('resonance', 'Resonance arrow'),
                              ('arrow_dashed', 'Dashed arrow'), ('curved', 'Curved arrow (electron pair)'),
                              ('curved_half', 'Fish-hook arrow (one electron)')]),
         ('group', 'shapes', [('brackets', 'Brackets [ ]'), ('parens', 'Parentheses ( )'),
                              ('braces', 'Braces { }')]),
         ('group', 'elements', ELEMENT_TIPS, 8),
         ('group', 'orbitals', ORBITAL_TIPS, 3),
         ('group', 'fragmentation', [('massfrag', 'Mass fragmentation: drag across bonds'),
                                     ('retrosyn', 'Retrosynthesis: drag across bonds of a product'),
                                     ('synthesis', 'Synthesis: drag across bonds of a product')]),
         ('tlc', 'TLC plate  (g)'),
         ('text', 'Text: captions, and atom labels typed in place  (t)'), None,
         ('sequence', 'Biopolymer sequence: click where it goes, then type or paste it (FASTA, H-Ala-Gly-OH, HELM)'),
         ('group', 'bio', list(biodraw.NAMES.items())))
RADIOLABELS = {'H': (3,), 'C': (11, 14), 'N': (13,), 'O': (15,), 'F': (18,), 'P': (32, 33), 'S': (35,),
               'I': (123, 124, 125, 131), 'Cl': (36,), 'Br': (76, 82)}
chem_sup = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')


def isotope_choices(el):
    """Mass numbers to offer for an element: its stable isotopes other than the usual one, then the usual
    radiolabels (14C, 3H, 18F, 125I ...)."""
    import chem as _chem
    try:
        usual = round(_chem.PT.GetMostCommonIsotopeMass(el))
        stable = [round(m) for m, _ab in _chem.isotopes(el)]
    except RuntimeError:
        return []
    out = [n for n in stable if n != usual] + [n for n in RADIOLABELS.get(el, ()) if n not in stable]
    return [n for n in out if model.isotope_mass(el, n)]


COLORS = (('Black', ''), ('Red', '#E01E1E'), ('Blue', '#1E4FE0'), ('Green', '#1E8C3A'), ('Orange', '#E07B00'),
          ('Magenta', '#C0189A'), ('Cyan', '#0096B4'), ('Grey', '#808080'))
HIGHLIGHT_COLORS = ('#2EE62E', '#2EE6F0', '#FFEA2E', '#F2994A', '#FF5252', '#B5E62E', '#8AD8FF', '#E8D44D',
                    '#B4B8FF', '#C27BFF', '#6BE04A', '#8CF5CF', '#FFBE80', '#FF9CF0', '#C4C4C4', '#FFC21A',
                    '#E08A8A', '#FF8A5C', '#C9B98A', '#FFD8C2')
RING_FILL_COLORS = ('#E5F7B5', '#C8F5D2', '#7FF28C', '#B6EB3B', '#33D633', '#FFF7B0', '#FFE14D', '#FFF23B',
                    '#FFC07A', '#FF9F40', '#CFF7F7', '#C9DDF7', '#3BF2F2', '#6FB8EB', '#1FCFCF', '#FFD0D0',
                    '#FBD7F7', '#FF8ADF', '#FF8A80', '#FF4F4F', '#ECECEC', '#E4D9F5', '#CBA8EF', '#B8B8B8',
                    '#9B6BE0')
ELEMENTS = ('C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I')
STYLE_FIELDS = (('bond_length', 'Bond length (pt)', 4, 60, 0.1), ('line_width', 'Line width (pt)', 0.1, 5, 0.05),
                ('bold_width', 'Bold / wedge width (pt)', 0.3, 10, 0.1), ('hash_spacing', 'Hash spacing (pt)', 0.5, 10, 0.1),
                ('bond_spacing', 'Double bond spacing (× bond length)', 0.05, 0.5, 0.01),
                ('margin', 'Margin around labels (pt)', 0, 6, 0.1), ('font_size', 'Atom label size (pt)', 4, 48, 0.5),
                ('caption_size', 'Caption size (pt)', 4, 72, 0.5))


class StyleDialog(QtWidgets.QDialog):
    """File > Document settings: a style sheet to start from (ACS 1996, RSC, Wiley, Nature, slides ...), then bond
    length, line widths, fonts, page size and orientation."""

    def __init__(self, style, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Document settings')
        theme.register_window(self)
        lay = QtWidgets.QFormLayout(self)
        self.sheet = QtWidgets.QComboBox()
        self.sheet.addItem('(as set below)')
        self.sheet.addItems(list(render.STYLE_SHEETS))
        self.sheet.setToolTip('Load the settings of a journal or document style; they can still be changed below')
        self.sheet.textActivated.connect(self._sheet_chosen)
        lay.addRow('Style sheet:', self.sheet)
        self.spins = {}
        for key, label, lo, hi, step in STYLE_FIELDS:
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(2)
            lay.addRow(label + ':', s)
            self.spins[key] = s
        self.font = QtWidgets.QFontComboBox()
        lay.addRow('Atom label font:', self.font)
        self.label_bold = QtWidgets.QCheckBox('Bold atom labels')
        lay.addRow('', self.label_bold)
        self.caption_font = QtWidgets.QFontComboBox()
        lay.addRow('Caption font:', self.caption_font)
        self.terminal = QtWidgets.QCheckBox('Show terminal carbon labels (CH3)')
        lay.addRow('', self.terminal)
        self.color_atoms = QtWidgets.QCheckBox('Colour heteroatoms (N blue, O red ...)')
        lay.addRow('', self.color_atoms)
        self.page = QtWidgets.QComboBox(); self.page.addItems(list(render.PAGES))
        self.orientation = QtWidgets.QComboBox(); self.orientation.addItems(['Portrait', 'Landscape'])
        row = QtWidgets.QHBoxLayout(); row.addWidget(self.page, 1); row.addWidget(self.orientation, 1)
        lay.addRow('Page:', row)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                         QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.accepted.connect(self.accept); box.rejected.connect(self.reject)
        lay.addRow(box)
        self.load(style)

    def _sheet_chosen(self, name):
        if name in render.STYLE_SHEETS:
            s = dict(render.DEFAULT_STYLE)
            s.update(render.STYLE_SHEETS[name])
            s['page'], s['orientation'] = self.page.currentText(), self.orientation.currentText().lower()
            self.load(s, keep_sheet=True)

    def load(self, style, keep_sheet=False):
        for key, s in self.spins.items():
            s.setValue(float(style[key]))
        self.font.setCurrentFont(QtGui.QFont(style['font']))
        self.caption_font.setCurrentFont(QtGui.QFont(style['caption_font']))
        self.label_bold.setChecked(bool(style['label_bold']))
        self.terminal.setChecked(bool(style['show_terminal_carbons']))
        self.color_atoms.setChecked(bool(style['color_atoms']))
        self.page.setCurrentText(style.get('page', 'A4'))
        self.orientation.setCurrentText(str(style.get('orientation', 'portrait')).capitalize())
        if not keep_sheet:
            self.sheet.setCurrentIndex(0)

    def style(self):
        """The settings shown (the document's other options, such as stereo labels, are left as they are)."""
        s = {k: sp.value() for k, sp in self.spins.items()}
        s['font'] = self.font.currentFont().family()
        s['caption_font'] = self.caption_font.currentFont().family()
        s['label_bold'] = self.label_bold.isChecked()
        s['show_terminal_carbons'] = self.terminal.isChecked()
        s['color_atoms'] = self.color_atoms.isChecked()
        s['page'] = self.page.currentText()
        s['orientation'] = self.orientation.currentText().lower()
        return s



def claim_keys(window, text, action):
    """Keeps a key that another program registers system-wide (Notion takes Ctrl+Shift+K) for this program while it is
    in front: the main window's action runs, in whichever main window is active."""
    app = QtWidgets.QApplication.instance()
    keys = getattr(app, 'kekulex_keys', None)
    if keys is None:
        keys = app.kekulex_keys = winkeys.KeyClaim(app)
        keys.windows = []
    keys.windows.append(window)

    def run():
        live = [w for w in keys.windows if not sip_deleted(w)]
        active = QtWidgets.QApplication.activeWindow()
        target = next((w for w in live if w is active or (active is not None and w.isAncestorOf(active))), None)
        target = target or (live[-1] if live else None)
        if target is not None:
            getattr(target, action).trigger()
    keys.claim(text, run)


def sip_deleted(obj):
    from PyQt6 import sip
    return sip.isdeleted(obj)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.path = None
        self.dirty = False
        self.settings = theme.settings()
        self.canvas = Canvas(self)
        self.autosave = autosave.Slot()                     # a copy of unsaved work, for crash recovery
        self._revision = self._autosaved = 0
        self.autosave_timer = QtCore.QTimer(self, interval=autosave.INTERVAL_S * 1000)
        self.autosave_timer.timeout.connect(self.autosave_now)
        self.autosave_timer.start()
        self.panel3d = None
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(self.canvas, 'Drawing')
        if features.installed(features.VIEW_3D):            # chosen in the installer
            from view3d import Panel3D
            self.panel3d = Panel3D(self)
            self.tabs.addTab(self.panel3d, '3D')
        else:
            self.tabs.tabBar().hide()
        self.tabs.currentChanged.connect(self._tab_changed)
        self.setCentralWidget(self.tabs)
        self.info = InfoPanel(self.canvas)
        dock = QtWidgets.QDockWidget('Formula and masses', self)
        dock.setObjectName('info')
        dock.setWidget(self.info)
        self.info_dock = dock
        dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetMovable |
                         QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.canvas.changed.connect(self._on_change)
        self.canvas.tool_key.connect(self.set_tool)
        self.canvas.message.connect(lambda m: self.statusBar().showMessage(m, 6000))
        self.canvas.selection_changed.connect(self._sync_style_bar)
        self.canvas.context_menu.connect(self.object_menu)
        self.canvas.edit_object.connect(self.edit_object)
        self.canvas.zoom_changed.connect(lambda z: self.zoom_label.setText(f'{z * 100:.0f} %'))
        self._build_tools()
        self._build_color_bars()
        self._build_menus()
        self.zoom_label = QtWidgets.QLabel('200 %')
        self.statusBar().addPermanentWidget(self.zoom_label)
        self.set_tool('bond1')
        self._drawing_bars = []
        self.setAcceptDrops(True)
        theme.register_window(self)
        theme.listeners.append(self._retheme)
        self.resize(1200, 760)
        theme.restore_geometry(self, 'main')
        self._title()

    # ------------------------------------------------------------------------------------------- building
    def _act(self, text, slot, shortcut=None, icon=None, tip=None):
        a = QtGui.QAction(text, self)
        if shortcut:
            a.setShortcut(QtGui.QKeySequence(shortcut))
        if icon:
            a.setIcon(icons.icon(icon)); a.icon_name = icon
        if tip:
            a.setToolTip(tip); a.setStatusTip(tip)
        a.triggered.connect(lambda _c=False: slot())
        return a

    def _build_tools(self):
        self.tool_actions = {}
        tb = QtWidgets.QToolBar('Tools', self)
        tb.setObjectName('tools')
        tb.setIconSize(QtCore.QSize(24, 24))
        tb.setMovable(False)
        self.addToolBar(QtCore.Qt.ToolBarArea.LeftToolBarArea, tb)
        group = QtGui.QActionGroup(self)
        group.setExclusionPolicy(QtGui.QActionGroup.ExclusionPolicy.ExclusiveOptional)   # "…" element: none checked
        palette =QtWidgets.QWidget()                       # two columns of tool buttons, then the elements
        g = QtWidgets.QGridLayout(palette)
        g.setContentsMargins(3, 3, 3, 3); g.setSpacing(2)
        row = col = 0

        def place(w):
            nonlocal row, col
            g.addWidget(w, row, col)
            col += 1
            if col == 2:
                row, col = row + 1, 0

        def gap():
            nonlocal row, col
            if col:
                row, col = row + 1, 0
            line = QtWidgets.QFrame(); line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
            g.addWidget(line, row, 0, 1, 2)
            row += 1
        def tool_action(name, tip):
            a = QtGui.QAction(icons.icon(name), tip, self, checkable=True)
            a.icon_name = name
            a.triggered.connect(lambda _c=False, n=name: self.set_tool(n))
            group.addAction(a)
            self.tool_actions[name] = a
            return a
        self.group_buttons = {}
        for t in TOOLS:
            if t is None:
                gap(); continue
            b = QtWidgets.QToolButton(); b.setIconSize(QtCore.QSize(24, 24)); b.setAutoRaise(True)
            if t[0] == 'group':                             # a drop-down: shows the last tool picked from it
                _g, gname, members = t[:3]
                cols = t[3] if len(t) > 3 else 0
                menu = QtWidgets.QMenu(b)
                grid = None
                if cols:                                    # a grid of buttons, like ChemDraw's palettes
                    holder = QtWidgets.QWidget()
                    grid = QtWidgets.QGridLayout(holder)
                    grid.setContentsMargins(4, 4, 4, 4); grid.setSpacing(1)
                    wa = QtWidgets.QWidgetAction(menu); wa.setDefaultWidget(holder)
                    menu.addAction(wa)
                for n, (name, tip) in enumerate(members):
                    a = tool_action(name, tip)
                    if grid is not None:
                        tb_ = QtWidgets.QToolButton(); tb_.setDefaultAction(a); tb_.setAutoRaise(True)
                        tb_.setIconSize(QtCore.QSize(24, 24)); tb_.setFixedSize(32, 30)
                        grid.addWidget(tb_, n // cols, n % cols)
                        a.triggered.connect(lambda _c=False, m=menu: m.hide())
                    else:
                        menu.addAction(a)
                    a.triggered.connect(lambda _c=False, a=a, b=b: b.setDefaultAction(a))
                b.setMenu(menu)
                b.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.MenuButtonPopup)
                b.setDefaultAction(self.tool_actions[members[0][0]])
                self.group_buttons[gname] = b
            else:
                b.setDefaultAction(tool_action(*t))
            b.setFixedSize(40, 32)
            place(b)
        gap()
        for el in ELEMENTS + ('…',):
            b = QtWidgets.QToolButton()
            b.setText(el); b.setCheckable(el != '…'); b.setFixedSize(40, 26); b.setAutoRaise(True)
            f = b.font(); f.setBold(True); b.setFont(f)
            b.setToolTip('Other element…' if el == '…' else f'{el}: click atoms to make them {el}')
            if el == '…':
                b.clicked.connect(self.other_element)
            else:
                act = QtGui.QAction(el, self, checkable=True)
                act.triggered.connect(lambda _c=False, e=el: self.set_tool('atom', e))
                group.addAction(act)
                self.tool_actions['atom:' + el] = act
                b.setDefaultAction(act)
                b.setToolTip(f'{el}: click atoms to make them {el}')
            place(b)
        g.setRowStretch(row + 1, 1)
        tb.addWidget(palette)

        top = self.addToolBar('File')
        self.file_bar = top
        top.setObjectName('file')
        top.setIconSize(QtCore.QSize(20, 20))
        top.setMovable(False)
        self.a_new = self._act('&New', self.new, 'Ctrl+N', 'new', 'New drawing')
        self.a_open = self._act('&Open…', self.open, 'Ctrl+O', 'open', 'Open a drawing or a MOL / SDF / SMILES file')
        self.a_save = self._act('&Save', self.save, 'Ctrl+S', 'save', 'Save the drawing (.sed)')
        self.a_undo = self._act('&Undo', self.canvas.undo, 'Ctrl+Z', 'undo', 'Undo')
        self.a_redo = self._act('&Redo', self.canvas.redo, 'Ctrl+Y', 'redo', 'Redo')
        self.a_copy_img = self._act('Copy as &picture', self.copy_picture, 'Ctrl+Shift+C', 'image',
                                    'Copy the selection (or everything) as a picture, for PowerPoint / Word')
        self.a_clean = self._act('C&lean up structure', self.canvas.clean_up, 'Ctrl+Shift+K', 'clean',
                                 'Lay out the selected structures (or all) again')
        self.a_clean.setShortcuts([QtGui.QKeySequence('Ctrl+Shift+K'), QtGui.QKeySequence('Ctrl+Alt+K')])
        self.a_clean.setShortcutContext(QtCore.Qt.ShortcutContext.ApplicationShortcut)   # also from the program's other windows
        claim_keys(self, 'Ctrl+Shift+K', 'a_clean')
        self.a_fit = self._act('&Fit to window', self.canvas.fit, 'Ctrl+0', 'fit', 'Show the whole drawing')
        self.a_send3d = self._act('Send to &3D', self.send_to_3d, 'Ctrl+3', 'send3d',
                                  'Show the selected molecule in the 3D tab')
        for a in (self.a_new, self.a_open, self.a_save, None, self.a_undo, self.a_redo, None, self.a_copy_img,
                  self.a_clean, self.a_fit) + ((self.a_send3d,) if self.panel3d else ()):
            top.addSeparator() if a is None else top.addAction(a)
        top.addSeparator()
        top.addWidget(QtWidgets.QLabel(' SMILES: '))
        self.smiles_box = QtWidgets.QLineEdit()
        self.smiles_box.setPlaceholderText('type or paste SMILES and press Enter')
        self.smiles_box.setMaximumWidth(360)
        self.smiles_box.returnPressed.connect(self.insert_smiles_box)
        top.addWidget(self.smiles_box)
        self._build_style_bar()

    def _build_color_bars(self):
        """The Highlight Colors and Ring Fill Colors palettes (View menu, or right-click): click a swatch to colour the
        selection's band / its rings; Custom… picks any colour, Remove takes it off."""
        self.highlight_bar = self._swatch_bar('Highlight Colors', 'highlight', HIGHLIGHT_COLORS, 4,
                                              self.canvas.set_highlight)
        self.fill_bar = self._swatch_bar('Ring Fill Colors', 'ringfill', RING_FILL_COLORS, 5, self.canvas.set_ring_fill)

    def _swatch_bar(self, title, key, colors, columns, apply):
        bar = QtWidgets.QToolBar(title, self)
        bar.setObjectName(key)
        bar.setToolTip(title)
        self.addToolBar(QtCore.Qt.ToolBarArea.LeftToolBarArea, bar)
        box = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(box)
        g.setContentsMargins(3, 3, 3, 3); g.setSpacing(3)
        g.addWidget(QtWidgets.QLabel(title.replace(' Colors', '')), 0, 0, 1, columns)
        for k, c in enumerate(colors):
            b = QtWidgets.QToolButton()
            b.setFixedSize(18, 18)
            b.setToolTip(c)
            b.setStyleSheet(f'QToolButton {{ background: {c}; border: 1px solid #8A9099; border-radius: 2px; }}'
                            f'QToolButton:hover {{ border: 2px solid #1B6CA8; }}')
            b.clicked.connect(lambda _c=False, c=c: apply(c))
            g.addWidget(b, 1 + k // columns, k % columns)
        row = 2 + (len(colors) - 1) // columns
        custom = QtWidgets.QPushButton('Custom…')
        custom.clicked.connect(lambda: self._custom_color(apply))
        remove = QtWidgets.QPushButton('Remove')
        remove.clicked.connect(lambda: apply(''))
        g.addWidget(custom, row, 0, 1, columns)
        g.addWidget(remove, row + 1, 0, 1, columns)
        bar.addWidget(box)
        bar.toggleViewAction().setText(f'Show {title} Toolbar')
        bar.setVisible(str(self.settings.value('show_' + key, 'false')).lower() == 'true')
        bar.visibilityChanged.connect(lambda v, k=key: self.settings.setValue('show_' + k, v))
        return bar

    def _custom_color(self, apply):
        c = QtWidgets.QColorDialog.getColor(QtGui.QColor('#FFF59D'), self, 'Colour')
        if c.isValid():
            apply(c.name().upper())

    def _transform_actions(self, m):
        """Flip / Rotate 180° / Rotate… / Scale…, as in ChemDraw's Object menu (also in the right-click menu)."""
        c = self.canvas
        if not hasattr(self, '_tf_actions'):
            self._tf_actions = [
                self._act('Flip &Horizontal', lambda: c.transform_selection('mirror_h'),
                          tip='The mirror image, left to right (stereocentres invert)'),
                self._act('Flip &Vertical', lambda: c.transform_selection('mirror_v'),
                          tip='The mirror image, top to bottom (stereocentres invert)'),
                self._act('Rotate 180° Hori&zontal', lambda: c.transform_selection('flip_h'), 'Ctrl+Shift+H',
                          tip='Turned over left to right: still the same molecule'),
                self._act('Rotate 180° V&ertical', lambda: c.transform_selection('flip_v'), 'Ctrl+Shift+V',
                          tip='Turned over top to bottom: still the same molecule'),
                self._act('R&otate…', self.rotate_dialog, 'Ctrl+R'),
                self._act('Rotate 15° clockwise', lambda: c.transform_selection('rot+'), 'Ctrl+]'),
                self._act('Rotate 15° counter-clockwise', lambda: c.transform_selection('rot-'), 'Ctrl+['),
                self._act('&Scale…', self.scale_dialog, 'Ctrl+K')]
        for a in self._tf_actions:
            m.addAction(a)

    def object_menu(self, global_pos, item):
        """The right-click menu, laid out like ChemDraw's: edit, arrange, colours, then a submenu per kind of object
        in the selection (Atom, Bond, Text, Bracket, Table, Molecule, Analysis)."""
        c = self.canvas
        m = QtWidgets.QMenu(self)
        theme.register_window(m)
        sel_kinds = {k for k, _i in c.sel}
        if not c.sel:                                       # empty page
            m.addAction(self.a_undo); m.addAction(self.a_redo)
            m.addSeparator()
            m.addAction('&Paste', self.paste)
            self._paste_special_menu(m.addMenu('Paste S&pecial'))
            m.addAction('Convert &Name to Structure…', self.name_to_structure)
            m.addAction('Select &All', c.select_all)
            m.addSeparator()
            m.addAction('Document &Settings…', self.style_dialog)
            m.addAction(self.highlight_bar.toggleViewAction())
            m.addAction(self.fill_bar.toggleViewAction())
            m.exec(global_pos)
            return
        m.addAction('Cu&t', self.cut); m.addAction('&Copy', self.copy); m.addAction('&Paste', self.paste)
        m.addAction('Clea&r', c.delete_selection)
        m.addSeparator()
        m.addAction(self.a_warnings)
        al = m.addMenu('&Align')
        for text, how in (('&Left edges', 'left'), ('Horizontal &centres', 'hcenter'), ('&Right edges', 'right'),
                          ('&Top edges', 'top'), ('&Vertical centres', 'vcenter'), ('&Bottom edges', 'bottom')):
            al.addAction(text, lambda h=how: c.align(h))
        di = m.addMenu('&Distribute')
        di.addAction('&Horizontally', lambda: c.distribute('h'))
        di.addAction('&Vertically', lambda: c.distribute('v'))
        if any(c.doc.group_of(it) is not None for it in c.sel):
            m.addAction(self.a_ungroup)
        else:
            m.addAction(self.a_group)
        fr = m.addMenu('Add &Frame')
        for text, kind in (('&Rectangle', 'rect'), ('R&ounded rectangle', 'roundrect'), ('&Brackets', 'brackets'),
                           ('&Parentheses', 'parens'), ('Bra&ces', 'braces')):
            fr.addAction(text, lambda k=kind: c.add_frame(k))
        self._transform_actions(m)
        m.addAction('Object &Settings…', self.style_dialog)
        m.addMenu(self._color_menu(m))
        m.addAction(self.highlight_bar.toggleViewAction())
        m.addAction(self.fill_bar.toggleViewAction())
        m.addSeparator()
        atoms = c.selected_atoms()
        if atoms:
            am = m.addMenu('A&tom')
            el = am.addMenu('&Element')
            for e in ELEMENTS:
                el.addAction(e, lambda e=e: c.edit_selected('atom', el=e, label='', iso=0))
            ch = am.addMenu('&Charge')
            for text, z in (('+1', 1), ('0', 0), ('−1', -1), ('+2', 2), ('−2', -2)):
                ch.addAction(text, lambda z=z: c.edit_selected('atom', charge=z))
            ch.addSeparator()
            for text, z, r in (('Radical cation •+', 1, 1), ('Radical anion •−', -1, 1)):
                ch.addAction(text, lambda z=z, r=r: c.edit_selected('atom', charge=z, radical=r))
            rm = am.addMenu('&Radical')
            for text, r in (('&None', 0), ('&Radical •', 1), ('&Carbene / diradical ••', 2)):
                a_ = rm.addAction(text, lambda r=r: c.edit_selected('atom', radical=r))
                a_.setCheckable(True)
                a_.setChecked(all(c.doc.atoms[a].radical == r for a in atoms))
            hm = am.addMenu('&Hydrogens')
            for text, h in (('&Automatic (from the valence)', None), ('&0', 0), ('&1', 1), ('&2', 2), ('&3', 3)):
                a_ = hm.addAction(text, lambda h=h: c.edit_selected('atom', hcount=h))
                a_.setCheckable(True)
                a_.setChecked(all(c.doc.atoms[a].hcount == h for a in atoms))
            els = {c.doc.atoms[a].el for a in atoms if not c.doc.atoms[a].label}
            if len(els) == 1:                               # 13C, 15N, D ...: the isotopes of that element
                el = els.pop()
                im = am.addMenu('&Isotope')
                for text, n in [('&Natural', 0)] + [(model.isotope_symbol(el, n) if el == 'H' else
                                                     str(n).translate(chem_sup) + el, n) for n in isotope_choices(el)]:
                    a_ = im.addAction(text, lambda n=n: c.edit_selected('atom', iso=n))
                    a_.setCheckable(True)
                    a_.setChecked(all(c.doc.atoms[a].iso == n for a in atoms))
            lm = am.addMenu('Lone &pairs')
            for text, n in (('&None', 0), ('&All it has', -1), ('&1', 1), ('&2', 2), ('&3', 3)):
                a_ = lm.addAction(text, lambda n=n: c.edit_selected('atom', lone_pairs=n))
                a_.setCheckable(True)
                a_.setChecked(all(c.doc.atoms[a].lone_pairs == n for a in atoms))
            if item and item[0] == 'atom':
                am.addAction('Edit &label…', lambda a=item[1]: c.edit_label(a))
            if any(c.doc.atoms[a].label for a in atoms):
                am.addAction('E&xpand label', c.expand_labels)
            if any(c.doc.atoms[a].charge_pos for a in atoms):
                am.addAction('&Reset charge position', c.reset_charge_positions)
            self._swatch_menu(am, '&Highlight', HIGHLIGHT_COLORS, c.set_highlight)
        bonds = [i for i, b in c.doc.bonds.items() if ('bond', i) in c.sel or (b.a in atoms and b.b in atoms)]
        if bonds:
            bm = m.addMenu('&Bond')
            om = bm.addMenu('&Order')
            for text, order in (('&Single', 1), ('&Double', 2), ('&Triple', 3)):
                om.addAction(text, lambda o=order: c.edit_selected('bond', order=o, stereo='', style='', side=None))
            sm = bm.addMenu('&Style')
            for text, stereo, style in (('&Plain', '', ''), ('&Wedged', 'wedge', ''), ('&Hashed wedged', 'hash', ''),
                                        ('&Bold', '', 'bold'), ('Hashe&d', '', 'hashed'), ('D&ashed', '', 'dashed'),
                                        ('Wa&vy', '', 'wavy')):
                sm.addAction(text, lambda st=stereo, sy=style: c.edit_selected(
                    'bond', stereo=st, style=sy, **({'order': 1} if st or sy not in ('', 'dashed') else {})))
            pm = bm.addMenu('Double bond &position')
            for text, side in (('&Automatic', None), ('&Left', -1), ('&Right', 1), ('&Centred', 0)):
                pm.addAction(text, lambda s=side: c.edit_selected('bond', side=s))
            bm.addAction('&Reverse direction', c.reverse_bonds)
            target = item[1] if item and item[0] == 'bond' else (bonds[0] if len(bonds) == 1 else None)
            if target is not None:
                self._cleave_menu(bm.addMenu('Clea&ve (MS fragments)'), target)
            self._swatch_menu(bm, '&Highlight', HIGHLIGHT_COLORS, c.set_highlight)
        if 'text' in sel_kinds:
            tm = m.addMenu('Te&xt')
            if item and item[0] == 'text':
                tm.addAction('&Edit…', lambda i=item[1]: c.edit_text(i))
            for text, key in (('&Bold', 'bold'), ('&Italic', 'italic'), ('&Formula (subscripts)', 'formula')):
                first = c.doc.texts[c.selected('text')[0]]
                a = tm.addAction(text, lambda k=key, v=not getattr(first, key): c.set_text_props(**{k: v}))
                a.setCheckable(True); a.setChecked(bool(getattr(first, key)))
            jm = tm.addMenu('&Justification')
            for text, al_ in (('&Left', 'left'), ('&Centre', 'center'), ('&Right', 'right')):
                jm.addAction(text, lambda a_=al_: c.set_text_props(align=a_))
        bio_shapes = [i for i in c.selected('shape') if c.doc.shapes[i].kind in model.BIO_KINDS]
        if bio_shapes:
            self._bio_object_menu(m, bio_shapes)
        if 'seq' in sel_kinds:
            self._biopolymer_menu(m.addMenu('Bio&polymer'))
        styled = [i for i in c.selected('shape') if c.doc.shapes[i].kind in
                  ('circle', 'oval', 'rect', 'roundrect', 'line', 'arc', 'orbital')]
        if styled:
            self._style_menu(m, styled)
        if 'tlc' in sel_kinds:
            self._tlc_menu(m, c.selected('tlc')[0], global_pos)
        frames = [i for i in c.selected('shape') if c.doc.shapes[i].kind in ('brackets', 'parens', 'braces')]
        if frames:
            km = m.addMenu('Brac&ket / shape')
            for text, kind in (('&Brackets', 'brackets'), ('&Parentheses', 'parens'), ('Bra&ces', 'braces'),
                               ('&Rectangle', 'rect'), ('R&ounded rectangle', 'roundrect'), ('O&val', 'oval')):
                km.addAction(text, lambda k=kind: self._set_shape_kind(k))
            km.addSeparator()
            km.addAction('Polymer &repeat count (n)…', lambda f=frames: self.set_repeat(f))
        if any(c.doc.arrows[i].kind.startswith('curved') for i in c.selected('arrow')):
            m.addAction('&Run the curved arrows (next step)', self.run_arrows)
        straight = [i for i in c.selected('arrow') if not c.doc.arrows[i].kind.startswith('curved')]
        if len(straight) == 1:
            m.addAction('Reaction &conditions (above / below the arrow)…',
                        lambda a=straight[0]: self.conditions_dialog(a))
        if 'grid' in sel_kinds:
            gm = m.addMenu('T&able')
            gid = c.selected('grid')[0]
            gm.addAction('&Edit values…', lambda: c.edit_grid(gid))
            gm.addAction('&Remove grid', c.delete_selection)
        if atoms:
            mm = m.addMenu('&Molecule')
            mm.addAction('&Select whole structure', lambda: c.select_fragment(next(iter(atoms))))
            mm.addAction(self.a_clean)
            mm.addAction('E&xpand labels', c.expand_labels)
            mm.addAction('&Check structure', self.check_structure)
            mm.addAction('Add &missing charges', c.fix_charges)
            if self.panel3d:
                mm.addAction(self.a_send3d)
            self._swatch_menu(mm, '&Ring fill', RING_FILL_COLORS, c.set_ring_fill)
            an = m.addMenu('A&nalysis')
            an.addAction('Show &analysis', c.add_analysis)
            an.addAction('Analyze &stoichiometry', c.analyze_stoichiometry)
            if nmr.available():
                an.addAction('Predict ¹&H NMR…', lambda: self.open_nmr('1H'))
                an.addAction('Predict ¹³&C NMR…', lambda: self.open_nmr('13C'))
            an.addAction('Convert Structure to &Name', self.structure_to_iupac)
            an.addAction('Common na&me (PubChem, online)', self.structure_to_name)
            an.addAction(self.a_stereo)
        m.exec(global_pos)

    # ------------------------------------------------------------------------------------------- styles, TLC
    def _style_menu(self, m, shapes):
        c = self.canvas
        kinds = {c.doc.shapes[i].kind for i in shapes}
        sm = m.addMenu('St&yle')
        if kinds & {'orbital'}:
            opts = (('&Plain', ''), ('&Shaded (phase)', 'shaded'), ('&Filled (phase)', 'filled'))
        elif kinds <= {'line', 'arc'}:
            opts = (('&Plain', ''), ('&Dashed', 'dashed'), ('&Bold', 'bold'), ('&Wavy (lines)', 'wavy'))
        else:
            opts = (('&Plain', ''), ('&Dashed', 'dashed'), ('&Shaded', 'shaded'), ('&Filled', 'filled'),
                    ('S&hadow', 'shadow'))
        for text, st in opts:
            sm.addAction(text, lambda st=st: self._set_shape_style(shapes, st))

    def _set_shape_style(self, shapes, style):
        c = self.canvas
        c.snapshot()
        for i in shapes:
            if i in c.doc.shapes:
                c.doc.shapes[i].style = style
        c.edited()

    def _tlc_menu(self, m, i, global_pos):
        """Right-click on a TLC plate: lanes, what the plate shows, the spots (ChemDraw's TLC menu)."""
        import tlc
        c = self.canvas
        plate = c.doc.tlcs[i]
        rend = c.renderer()
        x, y = c.to_model(QtCore.QPointF(c.mapFromGlobal(global_pos)))
        lane = tlc.lane_at(rend, i, x * rend.BL)
        spot = tlc.spot_at(rend, i, x * rend.BL, y * rend.BL)
        tm = m.addMenu('T&LC plate')

        def change(fn):
            c.snapshot(); fn(); c.edited()
        tm.addAction('&Edit Rf values…', lambda: self.edit_object('tlc', i))
        tm.addAction('&Add lane', lambda: change(lambda: plate.lanes.insert(lane + 1, [tlc.new_spot(0.5)])))
        tm.addAction('&Duplicate lane', lambda: change(lambda: plate.lanes.insert(
            lane + 1, [dict(s) for s in plate.lanes[lane]])))
        a = tm.addAction('De&lete lane', lambda: change(lambda: plate.lanes.pop(lane)))
        a.setEnabled(len(plate.lanes) > 1)
        tm.addSeparator()
        for text, attr in (('Show &origin', 'show_origin'), ('Show solvent &front', 'show_front'),
                           ('Show &border', 'show_border'), ('Show side &ticks', 'show_ticks')):
            a = tm.addAction(text, lambda at=attr: change(lambda: setattr(plate, at, not getattr(plate, at))))
            a.setCheckable(True); a.setChecked(getattr(plate, attr))
        sp = tm.addMenu('TLC &spots')
        every = [s for ln in plate.lanes for s in ln]
        a = sp.addAction('Show &Rf (all spots)', lambda v=not all(s['show_rf'] for s in every): change(
            lambda: [s.__setitem__('show_rf', v) for s in every]))
        a.setCheckable(True); a.setChecked(bool(every) and all(s['show_rf'] for s in every))
        for text, st in tlc.SPOT_STYLES:
            sp.addAction(f'All spots: {text}', lambda st=st: change(lambda: [s.__setitem__('style', st) for s in every]))
        if spot is not None:
            s = plate.lanes[spot[0]][spot[1]]
            sm = tm.addMenu(f'This spot (Rf {s["rf"]:.2f})')
            a = sm.addAction('Show &Rf', lambda: change(lambda: s.__setitem__('show_rf', not s['show_rf'])))
            a.setCheckable(True); a.setChecked(s['show_rf'])
            sm.addAction('&Set Rf…', lambda: self._set_spot_rf(s))
            for text, st in tlc.SPOT_STYLES:
                sm.addAction(text, lambda st=st: change(lambda: s.__setitem__('style', st)))
            self._swatch_menu(sm, '&Colour', [col for _n, col in COLORS if col] + list(HIGHLIGHT_COLORS[:6]),
                              lambda col: change(lambda: s.__setitem__('color', col)))
            sm.addAction('&Delete spot', lambda: change(lambda: (plate.lanes[spot[0]].pop(spot[1]),
                                                                 plate.lanes.pop(spot[0]) if not plate.lanes[spot[0]]
                                                                 and len(plate.lanes) > 1 else None)))

    def _set_spot_rf(self, spot):
        v, ok = QtWidgets.QInputDialog.getDouble(self, 'Rf', 'Rf (0 at the origin, 1 at the solvent front):',
                                                 spot['rf'], 0.0, 1.0, 2)
        if ok:
            self.canvas.snapshot(); spot['rf'] = round(v, 3); self.canvas.edited()

    # ------------------------------------------------------------------------------------------- biopolymers
    def _fill_bio_menu(self):
        """The Biopolymer menu's lower part follows the selection (rebuilt each time it opens)."""
        m = self.bio_menu
        for a in m.actions()[self._bio_menu_fixed:]:
            m.removeAction(a)
        self._biopolymer_menu(m)
        m.addSeparator()
        bm = m.addMenu('&BioDraw objects')
        for kind, name in biodraw.NAMES.items():
            bm.addAction(icons.icon(kind), name, lambda k=kind: self.set_tool(k))

    def _biopolymer_menu(self, m):
        """Sequence commands (menu bar and right-click): edit, style, numbers, expand, copy as, analysis."""
        c = self.canvas
        seqs = c.selected('seq')
        first = c.doc.seqs[seqs[0]] if seqs else None
        a = m.addAction('&Edit sequence…', lambda: self.edit_object('seq', seqs[0]) if seqs else None)
        a.setEnabled(bool(seqs))
        a = m.addAction('&Toggle standard amino acid style (A / Ala)',
                        lambda: c.set_sequence_props(three_letter=not first.three_letter))
        a.setEnabled(bool(first) and first.polymer == 'PEPTIDE')
        a = m.addAction('Show residue &numbers', lambda: c.set_sequence_props(show_ids=not first.show_ids))
        a.setCheckable(True); a.setChecked(bool(first and first.show_ids)); a.setEnabled(bool(seqs))
        a = m.addAction('&Residues per line…', self.sequence_per_line); a.setEnabled(bool(seqs))
        a = m.addAction('E&xpand to structure', c.expand_sequences); a.setEnabled(bool(seqs))
        a = m.addAction('Show &analysis (formula, mass, charge states)', self.sequence_analysis_caption)
        a.setEnabled(bool(seqs))
        ca = m.addMenu('&Copy As')
        ca.setEnabled(bool(seqs))
        ca.addAction('&HELM', lambda: self.copy_sequence('helm'))
        ca.addAction('HELM (&natural analog, FASTA)', lambda: self.copy_sequence('fasta'))
        ca.addAction('&Three-letter text', lambda: self.copy_sequence('three'))

    def _paste_special_menu(self, m):
        for text, polymer in (('FASTA &Peptide', 'PEPTIDE'), ('FASTA &DNA', 'DNA'), ('FASTA &RNA', 'RNA'),
                              ('&Biopolymer (H-Ala-Gly-Pro-NH2)', 'PEPTIDE'), ('&HELM', None)):
            m.addAction(text, lambda p=polymer: self.paste_sequence(p))

    def _bio_object_menu(self, m, shapes):
        c = self.canvas
        bm = m.addMenu('Bio &object')
        self._swatch_menu(bm, '&Fill colour', HIGHLIGHT_COLORS[:6] + RING_FILL_COLORS[:6], self.set_bio_fill)
        if any(c.doc.shapes[i].kind in model.BIO_LINEAR for i in shapes):
            bm.addAction('&Thicker', lambda: self.bio_thickness(1.2))
            bm.addAction('Th&inner', lambda: self.bio_thickness(1 / 1.2))
        plasmids = [i for i in shapes if c.doc.shapes[i].kind == 'plasmid']
        if plasmids:
            bm.addAction('&Plasmid map: size, regions, markers…', lambda: self.edit_object('shape', plasmids[0]))

    def set_bio_fill(self, color):
        c = self.canvas
        shapes = [i for i in c.selected('shape') if c.doc.shapes[i].kind in model.BIO_KINDS]
        if shapes:
            c.snapshot()
            for i in shapes:
                c.doc.shapes[i].fill = color
            c.edited()

    def bio_thickness(self, factor):
        c = self.canvas
        shapes = [i for i in c.selected('shape') if c.doc.shapes[i].kind in model.BIO_LINEAR]
        if shapes:
            c.snapshot()
            for i in shapes:
                sh = c.doc.shapes[i]
                sh.params = dict(sh.params, height=round(biodraw.thickness(sh) * factor, 3))
            c.edited()

    def edit_object(self, kind, i):
        """Double-click (or the sequence tool) on a sequence or a plasmid map: its editor."""
        from bio_dialog import PlasmidDialog, SequenceDialog
        c = self.canvas
        if kind == 'seq':
            seq = c.doc.seqs.get(i) if i >= 0 else None
            d = SequenceDialog(seq, self)
            if d.exec() != QtWidgets.QDialog.DialogCode.Accepted:
                return
            if seq is None:
                c.add_sequence(d.values(), getattr(c, 'pending_pos', None))
            else:
                c.update_sequence(i, d.values())
        elif kind == 'tlc' and i in c.doc.tlcs:
            from bio_dialog import TLCDialog
            d = TLCDialog(c.doc.tlcs[i], self)
            if d.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                c.snapshot()
                c.doc.tlcs[i].lanes = d.lanes
                c.edited()
        elif kind == 'shape' and i in c.doc.shapes:
            sh = c.doc.shapes[i]
            d = PlasmidDialog(sh.params, self)
            if d.exec() == QtWidgets.QDialog.DialogCode.Accepted:
                c.snapshot()
                sh.params = dict(sh.params, **d.params)
                c.edited()

    def insert_sequence(self):
        self.canvas.pending_pos = None
        self.edit_object('seq', -1)

    def paste_sequence(self, polymer, text=None):
        """Edit > Paste Special: the clipboard's text as a sequence (polymer None: HELM or detect)."""
        text = QtWidgets.QApplication.clipboard().text() if text is None else text
        try:
            pol, residues, extras = biopolymer.parse(text, polymer)
        except ValueError as e:
            QtWidgets.QMessageBox.information(self, APP_NAME, f'The clipboard does not hold a sequence: {e}')
            return
        self.canvas.add_sequence(dict(polymer=pol, residues=residues, **extras))

    def copy_sequence(self, how):
        """Edit > Copy As: HELM, FASTA (the natural analogs) or three-letter text of the selected sequences."""
        c = self.canvas
        seqs = c.selected('seq') or list(c.doc.seqs)
        if not seqs:
            self.statusBar().showMessage('Select a sequence first.', 4000)
            return
        out = []
        for i in seqs:
            q = c.doc.seqs[i]
            out.append(biopolymer.helm(q) if how == 'helm' else biopolymer.fasta(q) if how == 'fasta'
                       else biopolymer.three_letter(q))
        QtWidgets.QApplication.clipboard().setText('\n'.join(out))
        self.statusBar().showMessage('Copied: ' + out[0][:120], 5000)

    def sequence_per_line(self):
        c = self.canvas
        seqs = c.selected('seq')
        if not seqs:
            return
        n, ok = QtWidgets.QInputDialog.getInt(self, 'Residues per line', 'Residues per line:',
                                              c.doc.seqs[seqs[0]].per_line, 1, 500)
        if ok:
            c.set_sequence_props(per_line=n)

    def sequence_analysis_caption(self):
        """A caption under each selected sequence: formula, exact mass, MW and its charge states."""
        c = self.canvas
        for i in c.selected('seq'):
            q = c.doc.seqs[i]
            info = biopolymer.analyse(q)
            if info.errors:
                self.statusBar().showMessage(info.errors[0], 6000)
                continue
            sign = 1 if q.polymer == 'PEPTIDE' else -1
            ions = [ion for ion in info.ions if ion.z * sign > 0 and 'H]' in ion.label and '2M' not in ion.label
                    and 'H2O' not in ion.label]
            text = (f'{biopolymer.summary(q)}\nChemical Formula: {info.formula}\nExact Mass: {info.exact:.4f}\n'
                    f'Molecular Weight: {info.average:.2f}\n'
                    + ', '.join(f'{ion.label} {ion.mz:.4f}' for ion in ions[:6]))
            rend = c.renderer()
            rect = rend.item_rect('seq', i)
            c.snapshot()
            c.doc.add_text(rect.left() / rend.BL, rect.bottom() / rend.BL + 0.4, text)
            c.edited()

    def _swatch_menu(self, parent, title, colors, apply):
        sub = parent.addMenu(title)
        for c in colors[:12]:
            pm = QtGui.QPixmap(14, 14); pm.fill(QtGui.QColor(c))
            sub.addAction(QtGui.QIcon(pm), c, lambda c=c: apply(c))
        sub.addSeparator()
        sub.addAction('Custom…', lambda: self._custom_color(apply))
        sub.addAction('Remove', lambda: apply(''))

    def _cleave_menu(self, menu, b):
        """Bond > Cleave: each way of breaking bond b, named by its two fragments with their m/z."""
        import ms
        c = self.canvas
        opts = ms.options(c.doc, b)
        if not opts:
            menu.addAction('(a ring bond: breaking it leaves one piece)').setEnabled(False)
            return
        for k, (spec, how) in enumerate(opts):
            a, bb = (ms.describe(p) for p, _cut in ms.pieces(c.doc, b, spec))
            menu.addAction(f'{a}   +   {bb}     ({how.split(":")[0]})', lambda k=k: c.cleave_bond(b, k))

    def set_repeat(self, shapes):
        """Makes bracket pairs polymer brackets: n (or 3, x, m ...) at the foot of the right bracket; the atoms
        inside then count that many times in the formula and masses (a letter: shown as (unit)n, masses for 1)."""
        c = self.canvas
        shapes = [i for i in shapes if i in c.doc.shapes and c.doc.shapes[i].kind in model.BRACKET_KINDS]
        if not shapes:
            return
        cur = str(c.doc.shapes[shapes[0]].params.get('repeat', '')) or 'n'
        text, ok = QtWidgets.QInputDialog.getText(
            self, 'Polymer repeat count', 'Repeat count of the unit inside the brackets\n'
            '(a letter such as n, or a number; empty: plain brackets):', text=cur)
        if not ok:
            return
        text = text.strip()
        c.snapshot()
        for i in shapes:
            if text:
                c.doc.shapes[i].params['repeat'] = text
            else:
                c.doc.shapes[i].params.pop('repeat', None)
        c.edited()

    def _set_shape_kind(self, kind):
        c = self.canvas
        shapes = c.selected('shape')
        if shapes:
            c.snapshot()
            for i in shapes:
                c.doc.shapes[i].kind = kind
            c.edited()

    def _build_style_bar(self):
        """The text style bar: font size, bold, italic, formula (subscript numbers), alignment, colour - for the
        selected captions, or for the next ones typed."""
        self.addToolBarBreak()
        bar = self.addToolBar('Style')
        bar.setObjectName('style')
        bar.setMovable(False)
        bar.setIconSize(QtCore.QSize(18, 18))
        bar.addWidget(QtWidgets.QLabel(' Text: '))
        self.size_box = QtWidgets.QComboBox()
        self.size_box.setEditable(True)
        self.size_box.addItems(['Default'] + [str(v) for v in (6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 24, 36)])
        self.size_box.setToolTip('Caption size (pt); Default = the document font size')
        self.size_box.setFixedWidth(100)
        self.size_box.textActivated.connect(self._size_chosen)
        bar.addWidget(self.size_box)
        self.style_buttons = {}
        for key, text, tip in (('bold', 'B', 'Bold'), ('italic', 'I', 'Italic'),
                               ('formula', 'X₂', 'Formula: numbers after letters as subscripts (H2SO4)')):
            b = QtWidgets.QToolButton(); b.setText(text); b.setCheckable(True); b.setToolTip(tip); b.setAutoRaise(True)
            f = b.font(); f.setBold(key == 'bold'); f.setItalic(key == 'italic'); b.setFont(f)
            b.toggled.connect(lambda v, k=key: self._style_changed(k, v))
            bar.addWidget(b)
            self.style_buttons[key] = b
        bar.addSeparator()
        self.align_group = QtGui.QActionGroup(self)
        for key, tip in (('left', 'Align left'), ('center', 'Centre'), ('right', 'Align right')):
            a = QtGui.QAction(icons.icon('align_' + key), tip, self, checkable=True)
            a.icon_name = 'align_' + key
            a.setChecked(key == 'left')
            a.triggered.connect(lambda _c=False, k=key: self._style_changed('align', k))
            self.align_group.addAction(a); bar.addAction(a)
            self.style_buttons['align_' + key] = a
        bar.addSeparator()
        self.color_button = QtWidgets.QToolButton()
        self.color_button.setText('Colour'); self.color_button.setAutoRaise(True)
        self.color_button.setToolTip('Colour of the selected atoms, bonds, arrows, shapes and text')
        self.color_button.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.color_button.setMenu(self._color_menu())
        bar.addWidget(self.color_button)

    def _color_menu(self, parent=None):
        menu = QtWidgets.QMenu('&Color', parent or self)
        for name, c in COLORS:
            pm = QtGui.QPixmap(14, 14); pm.fill(QtGui.QColor(c or '#000000'))
            a = menu.addAction(QtGui.QIcon(pm), name)
            a.triggered.connect(lambda _c=False, c=c: self.canvas.set_color(c))
        menu.addSeparator()
        other = menu.addAction('Other…')
        other.triggered.connect(self._other_color)
        return menu

    def _other_color(self):
        c = QtWidgets.QColorDialog.getColor(QtGui.QColor('#000000'), self, 'Colour')
        if c.isValid():
            self.canvas.set_color(c.name().upper())

    def _size_chosen(self, text):
        try:
            size = 0.0 if text.lower().startswith('default') else float(text.replace(',', '.'))
        except ValueError:
            return
        self._style_changed('size', size)

    def _style_changed(self, key, value):
        if getattr(self, '_syncing', False):
            return
        self.canvas.set_text_props(**{key: value})

    def _sync_style_bar(self):
        """Show the style of the selected caption (or of the next one)."""
        texts = self.canvas.selected('text')
        props = ({k: getattr(self.canvas.doc.texts[texts[0]], k) for k in self.canvas.text_props} if texts
                 else self.canvas.text_props)
        self._syncing = True
        try:
            for k in ('bold', 'italic', 'formula'):
                self.style_buttons[k].setChecked(bool(props[k]))
            self.style_buttons['align_' + props['align']].setChecked(True)
            self.size_box.setCurrentText(f"{props['size']:g}" if props['size'] else 'Default')
        finally:
            self._syncing = False

    def _build_menus(self):
        mb = self.menuBar()
        m = mb.addMenu('&File')
        m.addAction(self.a_new)
        ns = m.addMenu('New from &Style Sheet')
        for name in render.STYLE_SHEETS:
            ns.addAction(name, lambda n=name: self.new_from_style_sheet(n))
        m.addAction(self.a_open)
        self.recent_menu = m.addMenu('Open &recent')
        self.recent_menu.aboutToShow.connect(self._fill_recent)
        m.addAction(self.a_save)
        m.addAction(self._act('Save &as…', self.save_as, 'Ctrl+Shift+S'))
        m.addSeparator()
        m.addAction(self._act('&Insert SMILES / MOL…', self.insert_dialog, 'Ctrl+I'))
        m.addAction(self._act('Insert &picture (PNG, JPEG …)…', self.insert_picture,
                              tip='A picture on the drawing; Ctrl+V also pastes a copied picture'))
        m.addAction(self._act('Export &picture (PNG, SVG, PDF)…', self.export_picture, 'Ctrl+E'))
        m.addAction(self._act('Export &MOL…', self.export_mol))
        m.addAction(self._act('Send the drawing to the lab &notebook…', self.send_to_notebook,
                              tip='Put this drawing into an entry of the lab notebook (its reaction and molecules)'))
        m.addAction(self._act('&Batch export (SDF / SMILES list to pictures)…', self.batch_export))
        m.addSeparator()
        m.addAction(self._act('&Document settings (page size, style)…', self.style_dialog))
        m.addAction(self._act('Print Pre&view…', self.print_preview))
        m.addAction(self._act('&Print…', self.print_document, 'Ctrl+P'))
        m.addSeparator()
        m.addAction(self._act('E&xit', self.close, 'Ctrl+Q'))

        m = mb.addMenu('&Edit')
        m.addAction(self.a_undo); m.addAction(self.a_redo)
        m.addSeparator()
        m.addAction(self._act('Cu&t', self.cut, 'Ctrl+X'))
        m.addAction(self._act('&Copy', self.copy, 'Ctrl+C'))
        m.addAction(self.a_copy_img)
        m.addAction(self._act('Copy as &SMILES', self.copy_smiles, 'Ctrl+Alt+C'))
        m.addAction(self._act('&Paste', self.paste, 'Ctrl+V'))
        self._paste_special_menu(m.addMenu('Paste S&pecial'))
        ca = m.addMenu('Copy &As')
        ca.addAction('&HELM', lambda: self.copy_sequence('helm'))
        ca.addAction('HELM (&natural analog, FASTA)', lambda: self.copy_sequence('fasta'))
        ca.addAction('&Three-letter text (H-Ala-Gly-OH)', lambda: self.copy_sequence('three'))
        m.addAction(self._act('&Delete', self.canvas.delete_selection))
        m.addAction(self._act('Select &all', self.canvas.select_all, 'Ctrl+A'))
        m.addAction(self._act('&Find…', self.find_dialog, 'Ctrl+F',
                              tip='Select captions, labels or a substructure (SMILES / SMARTS)'))

        m = mb.addMenu('&View')
        m.addAction(self._act('Zoom &in', lambda: self.canvas.set_zoom(self.canvas.zoom * 1.25), 'Ctrl+='))
        m.addAction(self._act('Zoom &out', lambda: self.canvas.set_zoom(self.canvas.zoom / 1.25), 'Ctrl+-'))
        m.addAction(self._act('&Actual size (100 %)', lambda: self.canvas.set_zoom(1.0), 'Ctrl+1'))
        m.addAction(self.a_fit)
        m.addAction(self._act('Whole &page', self.canvas.show_page, 'Ctrl+Shift+0'))
        m.addSeparator()
        self.a_stereo = self._toggle(m, 'Show &stereochemistry (R/S, E/Z)', 'show_stereo')
        self.a_numbers = self._toggle(m, 'Show atom &numbers', 'show_atom_numbers')
        self.a_reaction = QtGui.QAction('Show &reaction interpretation', self, checkable=True)
        self.a_reaction.setToolTip('Reactants red, products magenta, above the arrow cyan, below blue')
        self.a_reaction.toggled.connect(lambda v: (setattr(self.canvas, 'show_reaction', v), self.canvas.update()))
        m.addAction(self.a_reaction)
        self.a_partial = QtGui.QAction('Show &partial charges', self, checkable=True)
        self.a_partial.setToolTip('Gasteiger partial charges as discs behind the atoms: red negative (electron rich), '
                                  'blue positive (electron poor); point at an atom for its value')
        self.a_partial.toggled.connect(self._show_partial)
        m.addAction(self.a_partial)
        m.addAction(self.highlight_bar.toggleViewAction())
        m.addAction(self.fill_bar.toggleViewAction())
        a = QtGui.QAction('Display &warnings (atoms over their valence)', self, checkable=True)
        a.setChecked(True)
        a.toggled.connect(lambda v: (setattr(self.canvas, 'show_warnings', v), self.canvas.update()))
        m.addAction(a)
        self.a_warnings = a
        a = QtGui.QAction('Show r&ulers', self, checkable=True)
        a.toggled.connect(lambda v: (setattr(self.canvas, 'show_rulers', v), self.canvas.update()))
        m.addAction(a)
        m.addSeparator()
        theme.add_menu(self, m, self._retheme)

        m = mb.addMenu('&Object')
        self.a_group = self._act('&Group', self.canvas.group_selection, 'Ctrl+G',
                                 tip='The selected objects select, move and are analysed as one')
        self.a_ungroup = self._act('&Ungroup', self.canvas.ungroup_selection, 'Ctrl+Shift+G')
        m.addAction(self.a_group); m.addAction(self.a_ungroup)
        m.addSeparator()
        al = m.addMenu('&Align')
        for text, how, key in (('&Left edges', 'left', 'Ctrl+Alt+L'), ('Horizontal &centres', 'hcenter', 'Ctrl+Alt+H'),
                               ('&Right edges', 'right', 'Ctrl+Alt+R'), ('&Top edges', 'top', 'Ctrl+Alt+T'),
                               ('&Vertical centres', 'vcenter', 'Ctrl+Alt+V'), ('&Bottom edges', 'bottom', 'Ctrl+Alt+B')):
            al.addAction(self._act(text, lambda h=how: self.canvas.align(h), key))
        di = m.addMenu('&Distribute')
        di.addAction(self._act('&Horizontally', lambda: self.canvas.distribute('h')))
        di.addAction(self._act('&Vertically', lambda: self.canvas.distribute('v')))
        m.addAction(self._act('Arrange as reaction &scheme', self.canvas.arrange_scheme,
                              tip='Structures and arrows in one row, centred on one line, evenly spaced'))
        m.addSeparator()
        self._transform_actions(m)
        m.addSeparator()
        m.addMenu(self._color_menu(m))

        m = mb.addMenu('&Structure')
        m.addAction(self.a_clean)
        m.addAction(self._act('&Expand label', self.canvas.expand_labels, 'Ctrl+Shift+E'))
        m.addAction(self._act('Show &analysis (formula, mass, m/z …)', self.canvas.add_analysis, 'Ctrl+Shift+A'))
        m.addAction(self._act('Analyze &stoichiometry', self.canvas.analyze_stoichiometry))
        m.addAction(self._act('Check s&tructure', self.check_structure))
        m.addAction(self._act('Add &missing charges', self.canvas.fix_charges,
                              tip='The formal charges the bonds call for: N with four bonds is N+, O with one bond '
                                  'and no H is O−'))
        fm = m.addMenu('&Forms')
        fm.addAction(self._act('&Tautomers', lambda: self.canvas.add_forms('tautomers'),
                               tip='The other tautomers of the selected structure, to its right'))
        fm.addAction(self._act('&Resonance forms', lambda: self.canvas.add_forms('resonance'),
                               tip='The other resonance structures of the selected structure, to its right'))
        rm = m.addMenu('Reactio&n')
        rm.addAction(self._act('&Map atoms', self.canvas.map_reaction,
                               tip='Number matching atoms in reactants and products (:1, :2 ...)'))
        rm.addAction(self._act('&Clear atom map', self.canvas.clear_map))
        rm.addAction(self._act('Check &balance', self.check_balance))
        rm.addAction(self._act('&Run the curved arrows (next step)', self.run_arrows, 'Ctrl+Shift+R',
                               tip='Draws the structure the electron-pushing arrows make, to the right; or says what '
                                   'is wrong with the arrows'))
        m.addAction(self._act('Match m/&z…', self.match_mz, 'Ctrl+M',
                              tip='A measured m/z against the structures drawn and elemental compositions'))
        m.addAction(self._act('Annotate mass s&pectrum…', self.annotate_spectrum,
                              tip='Open a measured spectrum (CSV, JCAMP-DX, MSP, MGF, mzML) and label its peaks '
                                  'with the ions of the structures drawn'))
        m.addAction(self._act('MS f&ragments', self.show_fragments,
                              tip='The fragments the structure can logically give in EI, ESI+ or ESI− (the '
                                  'Fragments tab of the formula panel)'))
        if self.panel3d:
            m.addAction(self.a_send3d)
        m.addSeparator()
        if nmr.available():
            m.addAction(self._act('Predict ¹&H NMR…', lambda: self.open_nmr('1H')))
            m.addAction(self._act('Predict ¹³&C NMR…', lambda: self.open_nmr('13C')))
        m.addAction(self._act('Convert Structure to &Name', self.structure_to_iupac, 'Ctrl+Alt+N',
                              tip='Systematic IUPAC name, made offline (Openclatura), pasted under the structure'))
        m.addAction(self._act('&Name to structure (PubChem, online)…', self.name_to_structure, 'Ctrl+Shift+N'))
        m.addAction(self._act('Common na&me (PubChem, online)', self.structure_to_name))
        m.addSeparator()
        m.addAction(self.a_stereo); m.addAction(self.a_numbers); m.addAction(self.a_reaction)

        m = mb.addMenu('Bio&polymer')
        m.addAction(self._act('&Insert sequence…', self.insert_sequence, 'Ctrl+Shift+Q',
                              tip='A peptide, DNA or RNA sequence: type or paste FASTA, H-Ala-Gly-OH or HELM'))
        self._paste_special_menu(m.addMenu('Paste S&pecial'))
        m.addSeparator()
        self.bio_menu = m
        m.aboutToShow.connect(self._fill_bio_menu)
        self._bio_menu_fixed = len(m.actions())

        m = mb.addMenu('&Templates')
        for cat, items in templates.TEMPLATES.items():
            sub = m.addMenu(cat)
            for name, smiles in items:
                a = sub.addAction(name)
                a.triggered.connect(lambda _c=False, s=smiles: self.insert_text(s))

        m = mb.addMenu('F&ormat')
        m.addAction(self._act('&Document settings (bond length, lines, font, page)…', self.style_dialog))
        aps = m.addMenu('&Apply Style Sheet')
        for name in render.STYLE_SHEETS:
            aps.addAction(name, lambda n=name: self.apply_style_sheet(n))
        res = m.addMenu('Picture &resolution (PNG, clipboard)')
        grp = QtGui.QActionGroup(self)
        cur = int(self.settings.value('dpi', 600))
        for dpi in (300, 600, 1200):
            a = QtGui.QAction(f'{dpi} dpi', self, checkable=True); a.setChecked(dpi == cur)
            a.triggered.connect(lambda _c=False, d=dpi: self.settings.setValue('dpi', d))
            grp.addAction(a); res.addAction(a)
        self.a_transparent = QtGui.QAction('&Transparent background (PNG, clipboard)', self, checkable=True)
        self.a_transparent.setChecked(str(self.settings.value('transparent', 'false')).lower() == 'true')
        self.a_transparent.toggled.connect(lambda v: self.settings.setValue('transparent', v))
        m.addAction(self.a_transparent)

        m = mb.addMenu('&Help')
        m.addAction(self._act('&How to draw', self.show_help, 'F1'))
        m.addAction(self._act('&Keyboard shortcuts', self.show_shortcuts, 'Ctrl+/'))
        m.addAction(self._act('&About', self.about))

    def _toggle(self, menu, text, key):
        """A checkable menu item for a document display option (kept in the drawing, as ChemDraw does)."""
        a = QtGui.QAction(text, self, checkable=True)
        a.triggered.connect(lambda v: self._set_doc_option(key, v))
        menu.addAction(a)
        return a

    def _set_doc_option(self, key, value):
        self.canvas.snapshot()
        self.canvas.doc.style[key] = bool(value)
        self.canvas.edited()

    def _sync_toggles(self):
        s = render.style_of(self.canvas.doc)
        self.a_stereo.setChecked(bool(s['show_stereo']))
        self.a_numbers.setChecked(bool(s['show_atom_numbers']))

    def conditions_dialog(self, arrow):
        """Reagents and conditions over and under a reaction arrow."""
        c = self.canvas
        have = c.conditions_of(arrow)
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle('Reaction conditions')
        form = QtWidgets.QFormLayout(dlg)
        boxes = {}
        for where, label in (('above', 'Above the arrow (reagents)'), ('below', 'Below (solvent, temperature, time)')):
            e = QtWidgets.QPlainTextEdit()
            e.setFixedHeight(56)
            if have[where] is not None:
                e.setPlainText(c.doc.texts[have[where]].text)
            boxes[where] = e
            form.addRow(label, e)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                        QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        form.addRow(bb)
        if dlg.exec():
            c.set_conditions(arrow, boxes['above'].toPlainText(), boxes['below'].toPlainText())

    def check_balance(self):
        import reactions
        arrows = [i for k, i in self.canvas.sel if k == 'arrow']
        text, _diff = reactions.balance(self.canvas.doc, arrows[0] if arrows else None)
        QtWidgets.QMessageBox.information(self, 'Check balance', text)

    def match_mz(self):
        from ms_dialog import MatchDialog
        if getattr(self, '_match', None) is None:
            self._match = MatchDialog(self)
        self._match.show(); self._match.raise_(); self._match.activateWindow()
        self._match.mz.setFocus()

    def annotate_spectrum(self, path=None):
        from spectrum_dialog import SpectrumDialog
        if getattr(self, '_spectrum', None) is None:
            self._spectrum = SpectrumDialog(self)
        self._spectrum.show(); self._spectrum.raise_(); self._spectrum.activateWindow()
        if path:
            self._spectrum.load_path(path)
        elif self._spectrum.spectra:
            self._spectrum.update_view()                    # the drawing may have changed since
        return self._spectrum

    def find_dialog(self):
        text, ok = QtWidgets.QInputDialog.getText(
            self, 'Find', 'Text in a caption, an atom label or element (OMe, Br),\n'
            'or a substructure as SMILES / SMARTS (C(=O)O, c1ccccc1, [OX2H]):', text=getattr(self, '_last_find', ''))
        if ok and text.strip():
            self._last_find = text.strip()
            self.canvas.find(text.strip())

    def run_arrows(self):
        _items, problems = self.canvas.run_arrows()
        if problems:
            QtWidgets.QMessageBox.warning(self, 'Curved arrows', 'These arrows cannot be run:\n\n• ' +
                                          '\n• '.join(problems))

    def show_fragments(self):
        """Structure > MS fragments: the formula panel, on its Fragments tab."""
        self.info_dock.show()
        self.info_dock.raise_()
        self.info.tabs.setCurrentWidget(self.info.frag_tab)

    def _show_partial(self, on):
        self.canvas.show_partial = on
        self.canvas.update()
        if on:
            self.statusBar().showMessage('Partial charges (Gasteiger): red δ−, blue δ+. Point at an atom for its '
                                         'value.', 8000)

    def check_structure(self):
        problems = self.canvas.check_structure()
        if problems:
            QtWidgets.QMessageBox.warning(self, 'Check structure', '\n'.join(problems))
        else:
            self.statusBar().showMessage('Check structure: no problems found.', 5000)

    def name_to_structure(self):
        """A structure from its name, looked up at PubChem (the name is sent there)."""
        import pubchem
        name, ok = QtWidgets.QInputDialog.getText(self, 'Name to structure',
                                                  'Compound name (looked up online at PubChem):')
        if not ok or not name.strip():
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            smiles = pubchem.smiles_for_name(name)
        except pubchem.PubChemError as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, APP_NAME, str(e))
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        if not smiles:
            QtWidgets.QMessageBox.information(self, APP_NAME, f'PubChem has no compound called "{name}".')
            return
        import layout2d
        self.insert_text(layout2d.usual_form(smiles))

    def open_nmr(self, nucleus):
        """The NMR prediction window for the selected structure (or the one structure drawn)."""
        if not nmr.available():
            QtWidgets.QMessageBox.information(self, APP_NAME, 'NMR prediction needs the CSP5 module, which is not '
                                              'installed with this copy of the program.')
            return
        ids = self.info._target()
        if not ids:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'Draw (or select) a structure first.')
            return
        if len(self.canvas.doc.fragments(ids)) > 1:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'Select one structure (click one of its bonds twice, or '
                                              'drag around it) for its NMR prediction.')
            return
        from nmr_dialog import NMRDialog
        d = NMRDialog(self.canvas, ids, nucleus, self)
        d.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)   # kept alive by its parent until closed
        d.show()

    # ------------------------------------------------------------------------------------------- 3D tab
    def _tab_changed(self, index):
        """The 3D tab has no use for the drawing tools: the file bar stays, the rest comes back with the drawing."""
        if index == 1:
            self._drawing_bars = [w for w in self.findChildren(QtWidgets.QToolBar) + [self.info_dock]
                                  if w.isVisible() and w is not self.file_bar]
            for w in self._drawing_bars:
                w.hide()
            self.statusBar().showMessage('Left drag: rotate. Right drag: pan. Wheel: zoom. Click atoms to measure. '
                                         'Double-click: starting view.')
        else:
            for w in self._drawing_bars:
                w.show()
            self._drawing_bars = []
            self.statusBar().clearMessage()
            self.canvas.setFocus()

    def send_to_3d(self):
        """Builds the selected molecule (the whole structure, when only part of it is selected) in the 3D tab."""
        if not self.panel3d:
            return
        c = self.canvas
        seqs = c.selected('seq')
        atoms = c.selected_atoms()
        if seqs and not atoms:
            mol, error = biopolymer.to_mol(c.doc.seqs[seqs[0]])
            name = 'Sequence'
        else:
            frags = c.doc.fragments()
            chosen = [f for f in frags if atoms & set(f)] if atoms else frags
            if not chosen:
                only = list(c.doc.seqs)
                if len(only) != 1:
                    QtWidgets.QMessageBox.information(self, APP_NAME, 'Draw (or select) a molecule first.')
                    return
                mol, error = biopolymer.to_mol(c.doc.seqs[only[0]])
                name = 'Sequence'
            elif len(chosen) > 1:
                QtWidgets.QMessageBox.information(self, APP_NAME, 'Select one molecule (click one of its atoms or bonds '
                                                  'with the select tool, or drag around it) to send it to 3D.')
                return
            else:
                doc, grown = c.doc.expanded()
                mol, _order, error = chem.to_mol(doc, chem.expand_ids(grown, chosen[0]))
                name = ''
        if mol is None:
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'This structure has no 3D model: {error}')
            return
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        name = (name + ' ' if name else '') + CalcMolFormula(mol)
        if not self.panel3d.send(mol, name):
            self.statusBar().showMessage('Still building the last 3D model: try again in a moment.', 5000)
            return
        self.tabs.setCurrentIndex(1)

    def structure_to_iupac(self):
        """The systematic name of the selected structure (offline, Openclatura), pasted under it as a caption."""
        ids = self.info._target()
        if not ids:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'Draw (or select) a structure first.')
            return
        name = chem.iupac_name(self.canvas.doc, ids)
        if not name:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'No name could be made for this structure. Check it '
                                              '(Structure > Check structure), or try the PubChem lookup.')
            return
        self.canvas.add_caption_below(ids, name)
        self.statusBar().showMessage('Name (IUPAC 2013 rules, generated by Openclatura): ' + name, 8000)

    def structure_to_name(self):
        """The name PubChem has for the selected structure (its InChIKey is sent there); pasted as a caption."""
        import pubchem
        ids = self.info._target()
        info = chem.analyse(self.canvas.doc, ids)
        if not info.inchikey:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'Select (or draw) one complete structure first.')
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            found = pubchem.name_for_inchikey(info.inchikey)
        except pubchem.PubChemError as e:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.warning(self, APP_NAME, str(e))
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        if not found:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'PubChem does not have this structure (a new compound '
                                              'has no name there).')
            return
        title, iupac = found
        choice = title
        if iupac and iupac.lower() != title.lower():
            choice, ok = QtWidgets.QInputDialog.getItem(self, 'Structure to name', 'Paste which name under the structure?',
                                                        [title, iupac], 0, False)
            if not ok:
                return
        self.canvas.add_caption_below(ids, choice)

    def rotate_dialog(self):
        deg, ok = QtWidgets.QInputDialog.getDouble(self, 'Rotate', 'Degrees (clockwise):', 90, -360, 360, 1)
        if ok:
            self.canvas.rotate(deg)

    def scale_dialog(self):
        pct, ok = QtWidgets.QInputDialog.getDouble(self, 'Scale', 'Percent:', 100, 5, 2000, 1)
        if ok:
            self.canvas.scale_selection(pct / 100)

    def _retheme(self):
        for a in list(self.tool_actions.values()) + self.findChildren(QtGui.QAction):
            name = getattr(a, 'icon_name', None)
            if name:
                a.setIcon(icons.icon(name))
        if self.panel3d:
            self.panel3d.update()

    # ------------------------------------------------------------------------------------------- tools
    def set_tool(self, name, element=None):
        self.canvas.set_tool(name, element)
        key = f'atom:{element}' if name == 'atom' else name
        if key in self.tool_actions:
            self.tool_actions[key].setChecked(True)
        self.statusBar().showMessage(self.canvas.hint())

    def other_element(self):
        text, ok = QtWidgets.QInputDialog.getText(self, 'Element', 'Element symbol (e.g. Na, Si, B, Se):')
        if ok and text.strip():
            what = model.parse_label(text.strip())
            if what is None or what[0] != 'element':
                QtWidgets.QMessageBox.warning(self, APP_NAME, f'"{text}" is not an element symbol.')
                return
            self.set_tool('atom', what[1])
            for a in self.tool_actions.values():
                if a.isChecked():
                    a.setChecked(False)

    # ------------------------------------------------------------------------------------------- state
    def _on_change(self):
        self.dirty = True
        self._revision += 1
        self._title()
        self._sync_toggles()
        self.a_undo.setEnabled(bool(self.canvas.undo_stack))
        self.a_redo.setEnabled(bool(self.canvas.redo_stack))

    def autosave_now(self):
        """The autosave timer: a copy of the drawing while it has unsaved changes made since the last copy."""
        if not self.dirty or self.canvas.doc.is_empty():
            self.autosave.clear()
            return
        if self._revision != self._autosaved and self.autosave.write(self.canvas.doc, self.path):
            self._autosaved = self._revision

    def _saved(self):
        """The drawing on screen is what is on disk (or a new, empty one): no copy is needed."""
        self.dirty = False
        self.autosave.clear()
        self._autosaved = self._revision

    def recover(self):
        """At start-up: offers the drawings an earlier Kekulex left unsaved when it stopped (a crash, a power cut);
        the newest opens here, unsaved, so Save puts it back where it came from."""
        found = autosave.orphans()
        if not found:
            return False
        sed, meta = found[0]
        when = time.strftime('%d %b %Y, %H:%M', time.localtime(meta['saved']))
        name = os.path.basename(meta.get('path') or '') or 'an untitled drawing'
        more = f'\n\n({len(found) - 1} more will be offered the next time Kekulex starts.)' if len(found) > 1 else ''
        box = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Question, APP_NAME,
                                    f'Kekulex closed without saving {name} (last copy {when}).\n\n'
                                    f'Open the recovered drawing?{more}', parent=self)
        keep = box.addButton('Recover', QtWidgets.QMessageBox.ButtonRole.AcceptRole)
        box.addButton('Discard', QtWidgets.QMessageBox.ButtonRole.DestructiveRole)
        box.exec()
        if box.clickedButton() is not keep:
            autosave.discard(sed)
            return False
        try:
            with open(sed, encoding='utf-8') as f:
                doc = Doc.from_dict(json.load(f))
        except (OSError, ValueError) as e:
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'The recovered copy could not be read:\n{e}')
            return False
        self.canvas.set_doc(doc)
        self.path = meta.get('path') or None
        self._on_change()                                   # unsaved until the user saves it
        autosave.discard(sed)
        self.autosave_now()
        self.statusBar().showMessage('Recovered drawing - save it to keep it.', 10000)
        return True

    def _title(self):
        name = os.path.basename(self.path) if self.path else 'Untitled'
        self.setWindowTitle(f'{"*" if self.dirty else ""}{name} - {APP_NAME}')

    def _maybe_save(self):
        if not self.dirty or self.canvas.doc.is_empty():
            return True
        r = QtWidgets.QMessageBox.question(
            self, APP_NAME, 'Save the changes to the drawing?',
            QtWidgets.QMessageBox.StandardButton.Save | QtWidgets.QMessageBox.StandardButton.Discard |
            QtWidgets.QMessageBox.StandardButton.Cancel)
        if r == QtWidgets.QMessageBox.StandardButton.Save:
            return self.save()
        return r == QtWidgets.QMessageBox.StandardButton.Discard

    def closeEvent(self, e):
        if not self._maybe_save():
            e.ignore(); return
        theme.save_geometry(self, 'main')
        if self.panel3d:
            self.panel3d.wait()
        self.autosave_timer.stop()
        self.autosave.clear()                               # closed on purpose: nothing to recover
        e.accept()

    def _dir(self):
        return str(self.settings.value('last_dir', os.path.expanduser('~/Documents')))

    def _remember(self, path):
        self.settings.setValue('last_dir', os.path.dirname(path))
        rec = [p for p in (self.settings.value('recent') or []) if p != path]
        if isinstance(rec, str):
            rec = [rec]
        self.settings.setValue('recent', ([path] + rec)[:10])

    def _fill_recent(self):
        self.recent_menu.clear()
        rec = self.settings.value('recent') or []
        if isinstance(rec, str):
            rec = [rec]
        for p in rec:
            a = self.recent_menu.addAction(p)
            a.triggered.connect(lambda _c=False, p=p: self.open(p))
        if not rec:
            self.recent_menu.addAction('(none)').setEnabled(False)

    # ------------------------------------------------------------------------------------------- files
    def new(self):
        if self._maybe_save():
            self.path = None
            self.canvas.set_doc(Doc())
            self._saved()
            self._title()

    def open(self, path=None):
        if not self._maybe_save():
            return
        if not path:
            path, _f = QtWidgets.QFileDialog.getOpenFileName(self, 'Open', self._dir(), OPEN_FILTER)
            if not path:
                return
        notes = []
        try:
            if path.lower().endswith('.sed'):
                with open(path, encoding='utf-8') as f:
                    doc = Doc.from_dict(json.load(f))
                self.path = path
            elif path.lower().endswith('.cdxml'):
                doc, notes = cdxml.read(path)
                self.path = path
            else:
                doc = Doc()
                x = 0.0
                for mol in chem.mols_from_file(path):
                    items = doc.insert_rdkit(mol, (0, 0))
                    pts = doc.points_of(items)
                    dx = x - min(p[0] for p in pts) if pts else x
                    doc.transform(items, lambda px, py: (px + dx, py))
                    x = max(p[0] for p in doc.points_of(items)) + 2 if pts else x
                self.path = None
        except Exception as e:
            applog.log.exception('open %s', path)
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'Could not open {os.path.basename(path)}:\n{e}')
            return
        self.canvas.set_doc(doc)
        self._saved()
        self._remember(path)
        self._title()
        self._sync_toggles()
        if path.lower().endswith('.cdx'):
            notes.append('the structures only (a .cdx file is read through RDKit; save it as .cdxml in ChemDraw to '
                         'keep text, arrows and layout)')
        if notes:
            self.statusBar().showMessage('Not brought over from ' + os.path.basename(path) + ': ' + '; '.join(notes), 12000)

    def save(self):
        if not self.path:
            return self.save_as()
        try:
            if self.path.lower().endswith('.cdxml'):
                cdxml.write(self.canvas.doc, self.path, f'{APP_NAME} {__version__}')
                doc = self.canvas.doc
                if doc.seqs or any(sh.kind in model.BIO_KINDS for sh in doc.shapes.values()):
                    self.statusBar().showMessage('Sequences and bio objects are not written to .cdxml: save as .sed '
                                                 'to keep them (or expand sequences to structures first).', 10000)
            else:
                with open(self.path, 'w', encoding='utf-8') as f:
                    json.dump(self.canvas.doc.to_dict(), f, indent=1)
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'Could not save:\n{e}')
            return False
        self._saved()
        self._remember(self.path)
        self._title()
        return True

    def save_as(self):
        start = self.path or os.path.join(self._dir(), 'structure.sed')
        path, flt = QtWidgets.QFileDialog.getSaveFileName(self, 'Save drawing', start, SAVE_FILTER)
        if not path:
            return False
        if not path.lower().endswith(('.sed', '.cdxml')):
            path += '.cdxml' if 'cdxml' in flt else '.sed'
        self.path = path
        return self.save()

    def _subset(self):
        """The selection (or everything) as its own Doc, with the drawing's style."""
        if not self.canvas.sel:
            return self.canvas.doc
        d = Doc.from_dict(self.canvas.doc.to_dict(list(self.canvas.sel)))
        return d

    def export_picture(self):
        doc = self._subset()
        if doc.is_empty():
            QtWidgets.QMessageBox.information(self, APP_NAME, 'There is nothing to export.')
            return
        base = os.path.splitext(os.path.basename(self.path))[0] if self.path else 'structure'
        path, flt = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Export picture', os.path.join(self._dir(), base + '.png'),
            'PNG picture (*.png);;SVG vector picture (*.svg);;PDF (*.pdf)')
        if not path:
            return
        ext = os.path.splitext(path)[1].lower()
        if ext not in ('.png', '.svg', '.pdf'):
            ext = '.svg' if 'SVG' in flt else '.pdf' if 'PDF' in flt else '.png'
            path += ext
        try:
            if ext == '.svg':
                render.save_svg(doc, path)
            elif ext == '.pdf':
                render.save_pdf(doc, path)
            else:
                render.save_png(doc, path, int(self.settings.value('dpi', 600)), self.a_transparent.isChecked())
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'Could not export:\n{e}')
            return
        self.settings.setValue('last_dir', os.path.dirname(path))
        self.statusBar().showMessage(f'Exported {path}', 5000)

    def batch_export(self):
        from batch_dialog import BatchDialog
        BatchDialog(self).exec()

    def send_to_notebook(self):
        import notebook_client
        doc = self.canvas.doc
        if not doc.atoms:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'Draw something first.')
            return
        mols = _notebook_molecules(doc)
        notebook_client.send_drawing(self, doc.to_dict(), mols, '.'.join(m['smiles'] for m in mols), program=APP_NAME)

    def export_mol(self):
        doc = self._subset()
        if not doc.atoms:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'There is no structure to export.')
            return
        base = os.path.splitext(os.path.basename(self.path))[0] if self.path else 'structure'
        path, _f = QtWidgets.QFileDialog.getSaveFileName(self, 'Export MOL', os.path.join(self._dir(), base + '.mol'),
                                                         'MOL file (*.mol)')
        if not path:
            return
        block, _order = doc.expanded()[0].to_molblock(title=base)    # nicknames written out in full
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(block)
        self.statusBar().showMessage(f'Exported {path}', 5000)

    # ------------------------------------------------------------------------------------------- clipboard
    def _smiles_of(self, doc):
        parts = []
        for frag in doc.fragments():
            s = chem.analyse(doc, frag).smiles
            if s:
                parts.append(s)
        return '.'.join(parts)

    def copy(self, picture=True):
        doc = self._subset()
        if doc.is_empty():
            return
        md = QtCore.QMimeData()
        md.setData(MIME, json.dumps(self.canvas.doc.to_dict(list(self.canvas.sel) or None)).encode('utf-8'))
        smi = self._smiles_of(doc)
        if smi:
            md.setText(smi)
        if picture:
            dpi = int(self.settings.value('dpi', 600))
            img = render.to_image(doc, dpi, False)
            md.setImageData(img)
            buf = QtCore.QBuffer(); buf.open(QtCore.QIODevice.OpenModeFlag.WriteOnly)
            render.to_image(doc, dpi, self.a_transparent.isChecked()).save(buf, 'PNG')
            md.setData('image/png', buf.data())
            md.setData('image/svg+xml', render.svg_bytes(doc))
        QtWidgets.QApplication.clipboard().setMimeData(md)
        return True

    def copy_picture(self):
        if self.copy(picture=True):
            self.statusBar().showMessage('Copied as a picture: paste it into PowerPoint or Word.', 4000)

    def copy_smiles(self):
        smi = self._smiles_of(self._subset())
        if smi:
            QtWidgets.QApplication.clipboard().setText(smi)
            self.statusBar().showMessage('SMILES copied: ' + smi, 4000)

    def cut(self):
        if self.canvas.sel and self.copy():
            self.canvas.delete_selection()

    def _select_inserted(self):
        """What was just pasted or inserted stays selected, under the select tool, ready to drag or turn."""
        if self.canvas.tool not in ('select', 'lasso'):
            self.set_tool('select')

    def paste(self):
        import pictures
        md = QtWidgets.QApplication.clipboard().mimeData()
        if md.hasFormat(MIME):
            self.canvas.insert_items(json.loads(bytes(md.data(MIME)).decode('utf-8')))
            self._select_inserted()
            return
        img = pictures.from_mime(md)
        if img is not None and not (md.hasText() and self._is_structure_text(md.text())):
            self.canvas.insert_image(img)                   # a copied picture (or image file): placed as it is
            self.statusBar().showMessage('Picture pasted: drag it to move it, Object > Scale to resize.', 6000)
        elif md.hasText():
            text = md.text().strip()
            if text.upper().startswith(('PEPTIDE1{', 'RNA1{')):     # HELM: a sequence
                self.paste_sequence(None, text)
            else:
                self.insert_text(md.text())

    @staticmethod
    def _is_structure_text(text):
        """True when clipboard text reads as SMILES / MOL (then it wins over a picture copied with it)."""
        text = (text or '').strip()
        if not text or len(text) > 200000:
            return False
        if 'M  END' in text:
            return chem.mol_from_text(text) is not None
        lines = [x.strip() for x in text.replace('\r', '\n').split('\n') if x.strip()]
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        try:
            return bool(lines) and len(lines) < 200 and all(chem.mol_from_text(x) is not None for x in lines)
        finally:
            RDLogger.EnableLog('rdApp.*')

    def insert_picture(self):
        """File > Insert picture: an image file placed on the drawing."""
        import pictures
        path, _f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Insert picture', self._dir(), 'Pictures (' + ' '.join('*' + e for e in pictures.IMAGE_EXTS) + ')')
        if not path:
            return
        img = QtGui.QImage(path)
        if img.isNull():
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'Could not read {os.path.basename(path)} as a picture.')
            return
        self.canvas.insert_image(img)

    def insert_text(self, text):
        mols = []
        if 'M  END' in text:
            m = chem.mol_from_text(text)
            mols = [m] if m else []
        else:
            for part in text.replace('\r', '\n').split('\n'):
                part = part.strip()
                if not part:
                    continue
                m = chem.mol_from_text(part)
                if m is None:
                    mols = []
                    break
                mols.append(m)
        if not mols:
            QtWidgets.QMessageBox.information(self, APP_NAME, 'This is not a SMILES or MOL structure:\n' + text[:200])
            return False
        self.canvas.insert_items(mols=mols)
        self._select_inserted()
        return True

    def insert_smiles_box(self):
        if self.insert_text(self.smiles_box.text()):
            self.smiles_box.clear()

    def insert_dialog(self):
        text, ok = QtWidgets.QInputDialog.getMultiLineText(self, 'Insert structure',
                                                          'SMILES (one per line) or a MOL block:')
        if ok and text.strip():
            self.insert_text(text)

    # ------------------------------------------------------------------------------------------- drag and drop
    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls() or e.mimeData().hasImage():
            e.acceptProposedAction()

    def dropEvent(self, e):
        import pictures
        img = pictures.from_mime(e.mimeData())
        if img is not None:                                 # a picture (or image file) dropped: placed where it fell
            pos = self.canvas.mapFrom(self, e.position().toPoint())
            self.canvas.insert_image(img, self.canvas.to_model(QtCore.QPointF(pos)))
            e.acceptProposedAction()
            return
        for url in e.mimeData().urls():
            if url.isLocalFile():
                self.open(url.toLocalFile())
                break

    # ------------------------------------------------------------------------------------------- dialogs
    def style_dialog(self):
        d = StyleDialog(render.style_of(self.canvas.doc), self)
        if d.exec():
            self.canvas.snapshot()
            self.canvas.doc.style = {**self.canvas.doc.style, **d.style()}
            self.canvas.edited()

    def apply_style_sheet(self, name):
        """Format > Apply Style Sheet: the drawing takes a journal / document style (structures scale with the bond
        length, as all coordinates are in bond lengths); page and display options stay."""
        self.canvas.snapshot()
        self.canvas.doc.style = {**self.canvas.doc.style, **render.STYLE_SHEETS[name]}
        self.canvas.edited()
        self.statusBar().showMessage(f'Style sheet: {name}', 4000)

    def new_from_style_sheet(self, name):
        if self._maybe_save():
            doc = Doc()
            doc.style = dict(render.STYLE_SHEETS[name])
            self.path = None
            self.canvas.set_doc(doc)
            self._saved()
            self._title()

    # ------------------------------------------------------------------------------------------- printing
    def _print_on(self, printer):
        """Draw the page on a printer: the document's page at its real size, shrunk only if the paper is smaller."""
        doc = self.canvas.doc
        W, H = render.page_size(render.style_of(doc))
        p = QtGui.QPainter(printer)
        paper = printer.pageLayout().fullRectPoints()
        k = min(1.0, paper.width() / W, paper.height() / H)
        dev = printer.resolution() / 72.0
        p.scale(dev * k, dev * k)
        render.Renderer(doc).paint(p)
        p.end()

    def _printer(self):
        from PyQt6.QtPrintSupport import QPrinter
        printer = QPrinter(QPrinter.PrinterMode.HighResolution)
        printer.setFullPage(True)
        s = render.style_of(self.canvas.doc)
        size = {'A0': QtGui.QPageSize.PageSizeId.A0, 'A1': QtGui.QPageSize.PageSizeId.A1,
                'A2': QtGui.QPageSize.PageSizeId.A2,
                'A4': QtGui.QPageSize.PageSizeId.A4, 'A3': QtGui.QPageSize.PageSizeId.A3,
                'A5': QtGui.QPageSize.PageSizeId.A5, 'B5': QtGui.QPageSize.PageSizeId.B5,
                'Letter': QtGui.QPageSize.PageSizeId.Letter, 'Legal': QtGui.QPageSize.PageSizeId.Legal,
                'Tabloid': QtGui.QPageSize.PageSizeId.Tabloid}.get(s['page'], QtGui.QPageSize.PageSizeId.A4)
        printer.setPageSize(QtGui.QPageSize(size))
        printer.setPageOrientation(QtGui.QPageLayout.Orientation.Landscape if s['orientation'] == 'landscape'
                                   else QtGui.QPageLayout.Orientation.Portrait)
        printer.setDocName(os.path.basename(self.path) if self.path else 'Structure')
        return printer

    def print_document(self):
        from PyQt6.QtPrintSupport import QPrintDialog
        if self.canvas.doc.is_empty():
            QtWidgets.QMessageBox.information(self, APP_NAME, 'There is nothing to print.')
            return
        printer = self._printer()
        dlg = QPrintDialog(printer, self)
        if dlg.exec():
            self._print_on(printer)
            self.statusBar().showMessage('Sent to the printer.', 4000)

    def print_preview(self):
        from PyQt6.QtPrintSupport import QPrintPreviewDialog
        dlg = QPrintPreviewDialog(self._printer(), self)
        dlg.paintRequested.connect(self._print_on)
        dlg.resize(900, 1000)
        dlg.exec()

    def show_help(self):
        text = ('<h3>Drawing</h3><ul>'
                '<li><b>Bond tools</b>: click empty space for a bond, click an atom to grow a bond from it, '
                'drag from an atom to point it (15° steps; drop it on an atom to close a ring). Click a bond to make '
                'it double, then triple. With the double-bond tool, clicking a double bond moves its second line.</li>'
                '<li><b>Wedge / hash</b>: start at the stereocentre; click the bond again to reverse it.</li>'
                '<li><b>Chain</b>: drag; the number of bonds follows the mouse.</li>'
                '<li><b>Rings</b>: click empty space, a bond (fused ring) or an atom (attached ring). Benzene puts '
                'alternating double bonds where they fit.</li>'
                '<li><b>Hotkeys</b> as in ChemDraw: point at an atom or a bond and press a key (o = O, n = N, '
                'Shift+O = OMe, 2 = C=O, a = phenyl; over a bond 2 = double, a = fused benzene, w = wedge). '
                'Help &gt; Keyboard shortcuts lists them all.</li>'
                '<li><b>Labels</b>: double-click an atom (or click it with the Text tool, or press Enter over it) and '
                'type: an element (N, Cl), with a charge (O-, NH3+), or a nickname (OMe, CO2H, Boc, Ph, NHAc, OTs …). '
                'Nicknames stay short on the drawing and count in full for the formula; Structure &gt; Expand label '
                'draws them out.</li>'
                '<li><b>Captions</b>: Text tool, click and type; Esc or a click elsewhere ends. The style bar sets size, '
                'bold, italic, formula subscripts (X₂), alignment and colour.</li>'
                '<li><b>Arrows and shapes</b>: the drop-down buttons hold reaction, equilibrium, retrosynthetic, '
                'resonance, dashed, curved and fish-hook arrows, brackets, rectangles, ovals and lines; drag to draw.</li>'
                '<li><b>Object</b> menu: align, distribute, flip, rotate, scale, colour. <b>Structure</b> menu: clean up, '
                'expand label, analysis caption, R/S and E/Z labels. <b>Templates</b>: rings, heterocycles, amino '
                'acids, sugars, nucleobases, solvents.</li>'
                '<li><b>Biology</b>: the sequence tool (or Biopolymer &gt; Insert sequence) takes FASTA, '
                'H-Ala-Gly-OH or HELM text for peptides, DNA and RNA; double-click a sequence to edit it (cyclic, '
                'disulfide / lactam bridges, termini); Biopolymer &gt; Expand to structure, Copy As HELM / FASTA, '
                'Show analysis with the charge states. The bio drop-down draws membranes, DNA, receptors, '
                'organelles, cells and plasmid maps (click or drag).</li>'
                '<li><b>ChemDraw files</b>:.cdxml opens and saves (structures, captions, arrows, shapes, style); '
                '.cdx opens the structures.</li>'
                '<li><b>Select</b>: click, Shift+click or drag a box; drag to move; drop one atom on another to merge '
                'them. Arrow keys nudge. Drag the round handle at the top right of the selection to turn it freely '
                '(Shift: 15° steps).</li>'
                '<li><b>SMILES box</b> (top): paste a SMILES and press Enter. Ctrl+V pastes SMILES or MOL text too.</li>'
                '<li>Right or middle mouse button drag: pan. Ctrl+wheel: zoom. Ctrl+0: fit.</li></ul>'
                '<h3>Figures</h3><p>Ctrl+Shift+C copies the selection (or everything) as a picture for PowerPoint or '
                'Word. File &gt; Export picture writes PNG, SVG or PDF. The default style is ACS 1996 '
                '(Format &gt; Drawing style).</p>'
                '<h3>Masses</h3><p>The panel on the right shows the formula, the monoisotopic (exact) and average '
                'mass and the m/z of common positive (black) and negative (blue) ions of the selected structure. '
                'Atoms over their usual valence get a red box.</p>')
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle('How to draw'); box.setTextFormat(QtCore.Qt.TextFormat.RichText); box.setText(text)
        box.exec()

    def show_shortcuts(self):
        """The hotkey tables, built from canvas.py's own tables (so they cannot drift apart)."""
        import canvas as cv

        def key(k):
            return {' ': 'Space'}.get(k, f'Shift+{k}' if k.isalpha() and k.isupper() else k)

        def table(title, rows):
            cells = ''.join(f'<tr><td nowrap><b>{key(k)}</b>&nbsp;&nbsp;</td><td nowrap>{v}</td></tr>' for k, v in rows)
            return f'<h4>{title}</h4><table cellspacing="0" cellpadding="2">{cells}</table>'
        sprout = {'C': 'bond', 'O': 'C=O (2)', 'Ph': 'phenyl', 'Cy': 'cyclohexyl', 'cPent': 'cyclopentyl',
                  'ethynyl': 'ethynyl', 'cPr': 'cyclopropyl', 'cBu': 'cyclobutyl', 'tBu': 'tert-butyl', 'SO2': 'sulfonyl',
                  'Cp': 'cyclopentadienyl'}
        atoms = sorted(cv.ATOM_KEYS.items(), key=lambda kv: (kv[0].lower(), kv[0].isupper()))
        sprouts = []
        for k, (g, order, stereo) in cv.SPROUT_KEYS.items():
            what = 'wedged bond' if stereo == 'wedge' else 'hashed bond' if stereo == 'hash' else \
                '=CH2 (methylidene)' if order == 2 and g == 'C' else sprout.get(g, g)
            sprouts.append((k, what))
        sprouts.append(('9', 'two methyls'))
        fuse = [(k, ('benzene' if n == 6 else 'cyclopentadiene') if ar else f'{n}-membered ring')
                for k, (n, ar) in cv.FUSE_KEYS.items()]
        bonds = [('1 / 2 / 3', 'single / double / triple'), ('l / r / c', 'double bond line left / right / centred'),
                 ('w', 'wedge (again: reverse)'), ('h, Shift+W', 'hashed wedge'), ('b', 'bold'),
                 ('d', 'dashed (on a double bond: dashed second line)'), ('y', 'wavy'), ('Shift+H', 'hashed')] + fuse
        tools = [(k, v) for k, v in cv.TOOL_KEYS.items()]
        menu = [('Ctrl+Shift+K or Ctrl+Alt+K', 'clean up structure'), *([('Ctrl+3', 'send the selected molecule to 3D')] if self.panel3d else []), ('Ctrl+Alt+C', 'copy as SMILES'),
                ('Ctrl+Shift+C', 'copy as picture'),
                ('Ctrl+Shift+H / V', 'rotate 180° horizontal / vertical (same molecule)'),
                ('Ctrl+R', 'rotate…'), ('Ctrl+] / Ctrl+[', 'rotate 15°'), ('Ctrl+K', 'scale…'),
                ('Ctrl+E', 'export picture'), ('Ctrl+I', 'insert SMILES / MOL'),
                ('Ctrl+= / Ctrl+- / Ctrl+0', 'zoom in / out / fit'), ('+ / -', 'charge (over an atom)'),
                ('.', 'radical over an atom: •, then •• (carbene), then none'),
                ('c n o s … again', 'the same element key over its atom: one H fewer each time (NH2, NH, N, NH2 …)'),
                ('t', 'type the label of the atom under the mouse (CH2, CH3•, O•+ ...)'),
                ('Enter', 'type the label of the atom under the mouse'),
                ('Ctrl+Shift+E', 'expand label'), ('Ctrl+Shift+A', 'show analysis'),
                ('Ctrl+Shift+R', 'run the curved arrows (next mechanism step)'),
                ('Ctrl+Alt+L / H / R', 'align left / centres / right'), ('Ctrl+Alt+T / V / B', 'align top / middles / bottom'),
                ('Delete', 'delete the item under the mouse or the selection'), ('Esc', 'clear the selection')]
        html = ('<p>Point at an atom or bond with the mouse and press the key (with any tool). Capital = Shift.</p>'
                '<table><tr><td valign="top">' + table('Over an atom: label', [(k, v) for k, v in atoms]) +
                '</td><td width="24"></td><td valign="top">' + table('Over an atom: add', sprouts) +
                table('Over a bond', bonds) + '</td><td width="24"></td><td valign="top">' +
                table('Anywhere: tool', [(k, v.replace('bond1', 'bond').replace('benzene', 'benzene ring'))
                                         for k, v in tools]) + table('Menus and others', menu) + '</td></tr></table>')
        d = QtWidgets.QDialog(self)
        d.setWindowTitle('Keyboard shortcuts')
        theme.register_window(d)
        lay = QtWidgets.QVBoxLayout(d)
        view = QtWidgets.QTextBrowser(); view.setHtml(html)
        lay.addWidget(view)
        d.resize(900, 640)
        d.exec()

    def about(self):
        from rdkit import rdBase
        QtWidgets.QMessageBox.about(self, APP_NAME, f'<b>{APP_NAME}</b> {__version__}<br>Chemical structures for figures, '
                                    f'with formula, exact mass and adduct m/z.<br><br>Chemistry: RDKit {rdBase.rdkitVersion}; '
                                    f'interface: Qt {QtCore.QT_VERSION_STR}.<br>Log: {applog.LOG_FILE}')


def _notebook_molecules(doc):
    """The drawing's molecules for the lab notebook: SMILES, InChIKey, formula, MW and reactant / product around
    the reaction arrow (stoich.reaction_columns), left to right."""
    import chem
    import stoich
    roles = {}
    cols = stoich.reaction_columns(doc)
    if cols:
        for frag in cols[0]:
            roles[tuple(sorted(frag))] = 'reactant'
        for frag in cols[1]:
            roles[tuple(sorted(frag))] = 'product'
    out = []
    for frag in sorted(doc.fragments(), key=lambda f: sum(doc.atoms[a].x for a in f) / len(f)):
        info = chem.analyse(doc, frag)
        if info.smiles:
            out.append(dict(smiles=info.smiles, inchikey=info.inchikey, formula=info.formula, mw=round(info.average, 4),
                            role=roles.get(tuple(sorted(frag)))))
    return out


def main(argv=None):
    argv = sys.argv if argv is None else argv
    applog.setup()
    app = theme.make_app(argv, APP_NAME)
    w = MainWindow()
    w.show()
    QtCore.QTimer.singleShot(0, w.canvas.fit)
    if not os.environ.get('KEKULEX_NO_RECOVERY'):
        w.recover()                                         # drawings left unsaved by a crash
    if len(argv) > 1 and os.path.exists(argv[1]):
        w.open(argv[1])
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
