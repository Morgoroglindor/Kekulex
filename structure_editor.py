"""
structure_editor.py - Structure editor: draw chemical structures for figures and slides and read off their formula,
exact mass and the m/z of the usual LC-MS adducts.

    python structure_editor.py [drawing.sed | file.mol | file.sdf | file.smi]

Left: the drawing tools and elements. Middle: the drawing. Right: formula, masses, adduct m/z, SMILES, InChIKey.
Drawings are saved as .sed (JSON); structures go out as MOL / SMILES, pictures as PNG, SVG and PDF or to the
clipboard (Ctrl+Shift+C) for PowerPoint / Word.
"""
import json
import os
import sys

from PyQt6 import QtCore, QtGui, QtWidgets

import applog
import chem
import icons
import render
import theme
from canvas import Canvas
from info_panel import InfoPanel
from model import Doc
from version import __version__

APP_NAME = 'Structure editor'
MIME = 'application/x-structure-editor'
OPEN_FILTER = ('Structures (*.sed *.mol *.sdf *.sd *.smi *.smiles);;Structure editor drawings (*.sed);;'
               'MOL / SDF files (*.mol *.sdf *.sd);;SMILES files (*.smi *.smiles *.txt);;All files (*)')
TOOLS = (('select', 'Select / move  (Esc)'), None,
         ('bond1', 'Single bond  (1 over a bond)'), ('bond2', 'Double bond  (2 over a bond)'),
         ('bond3', 'Triple bond  (3 over a bond)'), ('wedge', 'Wedged bond (toward you)'),
         ('hash', 'Hashed bond (away from you)'), ('chain', 'Chain'), None,
         ('benzene', 'Benzene'), ('ring3', 'Cyclopropane'), ('ring4', 'Cyclobutane'), ('ring5', 'Cyclopentane'),
         ('ring6', 'Cyclohexane'), ('ring7', 'Cycloheptane'), ('ring8', 'Cyclooctane'), None,
         ('charge+', 'Positive charge  (+ over an atom)'), ('charge-', 'Negative charge  (- over an atom)'),
         ('erase', 'Eraser  (Delete over an item)'), None,
         ('arrow', 'Reaction arrow'), ('equilibrium', 'Equilibrium arrow'), ('text', 'Text'))
ELEMENTS = ('C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I')
STYLE_FIELDS = (('bond_length', 'Bond length (pt)', 4, 60, 0.1), ('line_width', 'Line width (pt)', 0.1, 5, 0.05),
                ('bold_width', 'Bold / wedge width (pt)', 0.3, 10, 0.1), ('hash_spacing', 'Hash spacing (pt)', 0.5, 10, 0.1),
                ('bond_spacing', 'Double bond spacing (× bond length)', 0.05, 0.5, 0.01),
                ('margin', 'Margin around labels (pt)', 0, 6, 0.1), ('font_size', 'Font size (pt)', 4, 48, 0.5))


class StyleDialog(QtWidgets.QDialog):
    """Format > Drawing style: bond length, line widths, font (ACS 1996 is the default)."""

    def __init__(self, style, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Drawing style')
        theme.register_window(self)
        lay = QtWidgets.QFormLayout(self)
        self.spins = {}
        for key, label, lo, hi, step in STYLE_FIELDS:
            s = QtWidgets.QDoubleSpinBox(); s.setRange(lo, hi); s.setSingleStep(step); s.setDecimals(2)
            lay.addRow(label + ':', s)
            self.spins[key] = s
        self.font = QtWidgets.QFontComboBox()
        lay.addRow('Font:', self.font)
        self.color_atoms = QtWidgets.QCheckBox('Colour heteroatoms (N blue, O red ...)')
        lay.addRow('', self.color_atoms)
        acs = QtWidgets.QPushButton('ACS 1996 (journals)')
        acs.clicked.connect(lambda: self.load(render.ACS_STYLE))
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                         QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.addButton(acs, QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        box.accepted.connect(self.accept); box.rejected.connect(self.reject)
        lay.addRow(box)
        self.load(style)

    def load(self, style):
        for key, s in self.spins.items():
            s.setValue(float(style[key]))
        self.font.setCurrentFont(QtGui.QFont(style['font']))
        self.color_atoms.setChecked(bool(style['color_atoms']))

    def style(self):
        s = {k: sp.value() for k, sp in self.spins.items()}
        s['font'] = self.font.currentFont().family()
        s['color_atoms'] = self.color_atoms.isChecked()
        return s


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.path = None
        self.dirty = False
        self.settings = theme.settings()
        self.canvas = Canvas(self)
        self.setCentralWidget(self.canvas)
        self.info = InfoPanel(self.canvas)
        dock = QtWidgets.QDockWidget('Formula and masses', self)
        dock.setObjectName('info')
        dock.setWidget(self.info)
        dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetMovable |
                         QtWidgets.QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
        self.canvas.changed.connect(self._on_change)
        self.canvas.zoom_changed.connect(lambda z: self.zoom_label.setText(f'{z * 100:.0f} %'))
        self._build_tools()
        self._build_menus()
        self.zoom_label = QtWidgets.QLabel('200 %')
        self.statusBar().addPermanentWidget(self.zoom_label)
        self.set_tool('bond1')
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
        for t in TOOLS:
            if t is None:
                gap(); continue
            name, tip = t
            a = QtGui.QAction(icons.icon(name), tip, self, checkable=True)
            a.icon_name = name
            a.triggered.connect(lambda _c=False, n=name: self.set_tool(n))
            group.addAction(a)
            self.tool_actions[name] = a
            b = QtWidgets.QToolButton(); b.setDefaultAction(a); b.setIconSize(QtCore.QSize(24, 24)); b.setAutoRaise(True)
            b.setFixedSize(34, 32)
            place(b)
        gap()
        for el in ELEMENTS + ('…',):
            b = QtWidgets.QToolButton()
            b.setText(el); b.setCheckable(el != '…'); b.setFixedSize(34, 26); b.setAutoRaise(True)
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
        self.a_fit = self._act('&Fit to window', self.canvas.fit, 'Ctrl+0', 'fit', 'Show the whole drawing')
        for a in (self.a_new, self.a_open, self.a_save, None, self.a_undo, self.a_redo, None, self.a_copy_img,
                  self.a_clean, self.a_fit):
            top.addSeparator() if a is None else top.addAction(a)
        top.addSeparator()
        top.addWidget(QtWidgets.QLabel(' SMILES: '))
        self.smiles_box = QtWidgets.QLineEdit()
        self.smiles_box.setPlaceholderText('type or paste SMILES and press Enter')
        self.smiles_box.setMaximumWidth(360)
        self.smiles_box.returnPressed.connect(self.insert_smiles_box)
        top.addWidget(self.smiles_box)

    def _build_menus(self):
        mb = self.menuBar()
        m = mb.addMenu('&File')
        m.addAction(self.a_new); m.addAction(self.a_open)
        self.recent_menu = m.addMenu('Open &recent')
        self.recent_menu.aboutToShow.connect(self._fill_recent)
        m.addAction(self.a_save)
        m.addAction(self._act('Save &as…', self.save_as, 'Ctrl+Shift+S'))
        m.addSeparator()
        m.addAction(self._act('&Insert SMILES / MOL…', self.insert_dialog, 'Ctrl+I'))
        m.addAction(self._act('Export &picture (PNG, SVG, PDF)…', self.export_picture, 'Ctrl+E'))
        m.addAction(self._act('Export &MOL…', self.export_mol))
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
        m.addAction(self._act('&Delete', self.canvas.delete_selection))
        m.addAction(self._act('Select &all', self.canvas.select_all, 'Ctrl+A'))
        m.addSeparator()
        m.addAction(self._act('Flip &horizontally', lambda: self.canvas.transform_selection('flip_h'), 'Ctrl+Shift+H'))
        m.addAction(self._act('Flip &vertically', lambda: self.canvas.transform_selection('flip_v'), 'Ctrl+Shift+V'))
        m.addAction(self._act('Rotate 15° clockwise', lambda: self.canvas.transform_selection('rot+'), 'Ctrl+R'))
        m.addAction(self._act('Rotate 15° counter-clockwise', lambda: self.canvas.transform_selection('rot-'), 'Ctrl+Shift+R'))
        m.addAction(self.a_clean)

        m = mb.addMenu('&View')
        m.addAction(self._act('Zoom &in', lambda: self.canvas.set_zoom(self.canvas.zoom * 1.25), 'Ctrl+='))
        m.addAction(self._act('Zoom &out', lambda: self.canvas.set_zoom(self.canvas.zoom / 1.25), 'Ctrl+-'))
        m.addAction(self._act('&Actual size (100 %)', lambda: self.canvas.set_zoom(1.0), 'Ctrl+1'))
        m.addAction(self.a_fit)
        m.addSeparator()
        theme.add_menu(self, m, self._retheme)

        m = mb.addMenu('F&ormat')
        m.addAction(self._act('Drawing &style…', self.style_dialog))
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
        m.addAction(self._act('&About', self.about))

    def _retheme(self):
        for a in list(self.tool_actions.values()) + self.findChildren(QtGui.QAction):
            name = getattr(a, 'icon_name', None)
            if name:
                a.setIcon(icons.icon(name))

    # ------------------------------------------------------------------------------------------- tools
    def set_tool(self, name, element=None):
        self.canvas.set_tool(name, element)
        key = f'atom:{element}' if name == 'atom' else name
        if key in self.tool_actions:
            self.tool_actions[key].setChecked(True)
        self.statusBar().showMessage(self.canvas.hint())

    def other_element(self):
        from canvas import parse_label
        text, ok = QtWidgets.QInputDialog.getText(self, 'Element', 'Element symbol (e.g. Na, Si, B, Se):')
        if ok and text.strip():
            el, _c = parse_label(text.strip())
            if el is None:
                QtWidgets.QMessageBox.warning(self, APP_NAME, f'"{text}" is not an element symbol.')
                return
            self.set_tool('atom', el)
            for a in self.tool_actions.values():
                if a.isChecked():
                    a.setChecked(False)

    # ------------------------------------------------------------------------------------------- state
    def _on_change(self):
        self.dirty = True
        self._title()
        self.a_undo.setEnabled(bool(self.canvas.undo_stack))
        self.a_redo.setEnabled(bool(self.canvas.redo_stack))

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
            self.dirty = False
            self._title()

    def open(self, path=None):
        if not self._maybe_save():
            return
        if not path:
            path, _f = QtWidgets.QFileDialog.getOpenFileName(self, 'Open', self._dir(), OPEN_FILTER)
            if not path:
                return
        try:
            if path.lower().endswith('.sed'):
                with open(path, encoding='utf-8') as f:
                    doc = Doc.from_dict(json.load(f))
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
        self.dirty = False
        self._remember(path)
        self._title()

    def save(self):
        if not self.path:
            return self.save_as()
        try:
            with open(self.path, 'w', encoding='utf-8') as f:
                json.dump(self.canvas.doc.to_dict(), f, indent=1)
        except OSError as e:
            QtWidgets.QMessageBox.warning(self, APP_NAME, f'Could not save:\n{e}')
            return False
        self.dirty = False
        self._remember(self.path)
        self._title()
        return True

    def save_as(self):
        start = self.path or os.path.join(self._dir(), 'structure.sed')
        path, _f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save drawing', start, 'Structure editor drawing (*.sed)')
        if not path:
            return False
        if not path.lower().endswith('.sed'):
            path += '.sed'
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
        block, _order = doc.to_molblock(title=base)
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

    def paste(self):
        md = QtWidgets.QApplication.clipboard().mimeData()
        if md.hasFormat(MIME):
            self.canvas.insert_items(json.loads(bytes(md.data(MIME)).decode('utf-8')))
        elif md.hasText():
            self.insert_text(md.text())

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
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e):
        for url in e.mimeData().urls():
            if url.isLocalFile():
                self.open(url.toLocalFile())
                break

    # ------------------------------------------------------------------------------------------- dialogs
    def style_dialog(self):
        d = StyleDialog(render.style_of(self.canvas.doc), self)
        if d.exec():
            self.canvas.snapshot()
            self.canvas.doc.style = d.style()
            self.canvas.edited()

    def show_help(self):
        text = ('<h3>Drawing</h3><ul>'
                '<li><b>Bond tools</b>: click empty space for a bond, click an atom to grow a bond from it, '
                'drag from an atom to point it (15° steps; drop it on an atom to close a ring). Click a bond to make '
                'it double, then triple. With the double-bond tool, clicking a double bond moves its second line.</li>'
                '<li><b>Wedge / hash</b>: start at the stereocentre; click the bond again to reverse it.</li>'
                '<li><b>Chain</b>: drag; the number of bonds follows the mouse.</li>'
                '<li><b>Rings</b>: click empty space, a bond (fused ring) or an atom (attached ring). Benzene puts '
                'alternating double bonds where they fit.</li>'
                '<li><b>Atoms</b>: point at an atom and press C, N, O, S, P, F, H, I, L (Cl), B (Br), + or -; '
                'or pick an element on the left and click atoms. Double-click an atom to type any element and charge.</li>'
                '<li><b>Select</b>: click, Shift+click or drag a box; drag to move; drop one atom on another to merge '
                'them. Arrow keys nudge.</li>'
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

    def about(self):
        from rdkit import rdBase
        QtWidgets.QMessageBox.about(self, APP_NAME, f'<b>{APP_NAME}</b> {__version__}<br>Chemical structures for figures, '
                                    f'with formula, exact mass and adduct m/z.<br><br>Chemistry: RDKit {rdBase.rdkitVersion}; '
                                    f'interface: Qt {QtCore.QT_VERSION_STR}.<br>Log: {applog.LOG_FILE}')


def main(argv=None):
    argv = sys.argv if argv is None else argv
    applog.setup()
    app = theme.make_app(argv, APP_NAME)
    w = MainWindow()
    w.show()
    QtCore.QTimer.singleShot(0, w.canvas.fit)
    if len(argv) > 1 and os.path.exists(argv[1]):
        w.open(argv[1])
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
