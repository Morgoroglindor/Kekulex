"""
info_panel.py - the panel beside the drawing: formula, exact and average mass, the m/z of the usual adducts, SMILES
and InChIKey of the selected structure (or of one picked from the list), updated as the drawing changes.
"""
from PyQt6 import QtCore, QtGui, QtWidgets

import chem


class InfoPanel(QtWidgets.QWidget):
    def __init__(self, canvas, parent=None):
        super().__init__(parent)
        self.canvas = canvas
        self.setMinimumWidth(260)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Structure:'))
        self.which = QtWidgets.QComboBox()
        self.which.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.which.setToolTip('Which structure the numbers are for. Selecting atoms picks "Selection".')
        self.which.currentIndexChanged.connect(lambda _i: self.refresh(keep_choice=True))
        row.addWidget(self.which, 1)
        lay.addLayout(row)
        self.tabs = QtWidgets.QTabWidget()
        lay.addWidget(self.tabs, 1)
        page = QtWidgets.QWidget()
        self.tabs.addTab(page, 'Analysis')
        outer, lay = lay, QtWidgets.QVBoxLayout(page)
        lay.setContentsMargins(4, 6, 4, 4)
        self.formula = QtWidgets.QLabel('–')
        f = self.formula.font(); f.setPointSizeF(f.pointSizeF() * 1.6); self.formula.setFont(f)
        self.formula.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        frow = QtWidgets.QHBoxLayout()
        frow.addWidget(self.formula, 1); frow.addWidget(self._paste_button('formula'))
        lay.addLayout(frow)
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        self.exact = self._value(); self.average = self._value(); self.charge = self._value()
        self.analysis = self._value(); self.analysis.setWordWrap(True)
        self.mz = self._value(); self.mz.setWordWrap(True)
        self.name = self._value(); self.name.setWordWrap(True)
        self.name.setToolTip('Systematic (IUPAC) name, generated offline by Openclatura from the IUPAC 2013 rules. '
                             'Check it before publishing.')
        for label, w, key in (('Name:', self.name, 'name'), ('Exact mass:', self.exact, 'exact'), ('m/z:', self.mz, 'mz'),
                              ('Mol. weight:', self.average, 'mw'), ('Elem. anal.:', self.analysis, 'analysis')):
            r = QtWidgets.QHBoxLayout(); r.addWidget(w, 1); r.addWidget(self._paste_button(key))
            form.addRow(label, r)
        form.addRow('Charge:', self.charge)
        lay.addLayout(form)
        self.table = QtWidgets.QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(['Ion', 'm/z'])
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setToolTip('Monoisotopic m/z (electron mass included). Ctrl+C copies the selected rows.')
        QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Copy, self.table, self.copy_rows,
                        context=QtCore.Qt.ShortcutContext.WidgetShortcut)
        lay.addWidget(self.table, 2)
        self.pattern_title = QtWidgets.QLabel('Isotope pattern')
        lay.addWidget(self.pattern_title)
        self.pattern = QtWidgets.QTableWidget(0, 2)
        self.pattern.setHorizontalHeaderLabels(['m/z', '%'])
        self.pattern.verticalHeader().hide()
        self.pattern.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.pattern.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.pattern.setToolTip('Isotope peaks of the ion selected above (relative to the highest)')
        lay.addWidget(self.pattern, 1)
        self.table.currentCellChanged.connect(lambda r, _c, _pr, _pc: self.show_pattern(r))
        for label, attr in (('SMILES', 'smiles'), ('InChIKey', 'inchikey')):
            lay.addWidget(QtWidgets.QLabel(label + ':'))
            r = QtWidgets.QHBoxLayout()
            e = QtWidgets.QLineEdit(); e.setReadOnly(True)
            b = QtWidgets.QToolButton(); b.setText('Copy')
            b.clicked.connect(lambda _c=False, e=e: QtWidgets.QApplication.clipboard().setText(e.text()))
            r.addWidget(e, 1); r.addWidget(b)
            lay.addLayout(r)
            setattr(self, attr, e)
        self.errors = QtWidgets.QLabel()
        self.errors.setWordWrap(True)
        self.errors.setStyleSheet('color: #C62828;')
        outer.addWidget(self.errors)
        self.props = QtWidgets.QTableWidget(0, 2)
        self.props.setHorizontalHeaderLabels(['Property', 'Value'])
        self.props.verticalHeader().hide()
        self.props.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.props.horizontalHeader().setStretchLastSection(True)
        self.props.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.props.setAlternatingRowColors(True)
        self.props.setToolTip('Calculated with RDKit. Double-click a row to paste it under the structure.')
        self.props.cellDoubleClicked.connect(self._paste_property_row)
        self.tabs.addTab(self.props, 'Properties')
        self.tabs.currentChanged.connect(lambda _i: self.refresh(keep_choice=True))
        self.timer = QtCore.QTimer(self, singleShot=True, interval=120)
        self.timer.timeout.connect(self.refresh)
        canvas.changed.connect(self.timer.start)
        canvas.selection_changed.connect(self.timer.start)
        self.refresh()

    def _paste_button(self, key):
        b = QtWidgets.QToolButton()
        b.setText('Paste'); b.setAutoRaise(True)
        b.setToolTip('Put this under the structure on the drawing, as a caption')
        b.clicked.connect(lambda _c=False, k=key: self.paste(k))
        return b

    def _target(self):
        choices = self._choices()
        k = self.which.currentIndex()
        ids = choices[k][1] if 0 <= k < len(choices) else None
        return set(ids) if ids else set(self.canvas.doc.atoms)

    def paste(self, key):
        """A caption with one analysis value under the structure the numbers are for."""
        info = getattr(self, 'info', None)
        if not info or not info.atoms:
            return
        if key == 'name':
            if not getattr(self, 'iupac', ''):
                return
            self.canvas.add_caption_below(self._target(), self.iupac)
            return
        if key == 'formula':
            text = f'Chemical Formula: {info.formula}{chem.charge_text(info.charge)}'
        elif key == 'exact':
            text = f'Exact Mass: {info.exact:.2f}'
        elif key == 'mw':
            text = f'Molecular Weight: {info.average:.2f}'
        elif key == 'mz':
            text = 'm/z: ' + ', '.join(f'{m:.2f} ({p:.1f}%)' for m, p in chem.isotope_pattern(info.counts, 0)[:4])
        else:
            text = f'Elemental Analysis: {info.analysis}'
        self.canvas.add_caption_below(self._target(), text, formula=True)

    def _paste_property_row(self, row, _col):
        name, value = self.props.item(row, 0).text(), self.props.item(row, 1).text()
        self.canvas.add_caption_below(self._target(), f'{name}: {value}')

    @staticmethod
    def _value():
        v = QtWidgets.QLabel('–')
        v.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        return v

    def _choices(self):
        """[(label, atom ids or None)]: the selection (if any), each structure, and all of them together."""
        doc = self.canvas.doc
        out = []
        sel = self.canvas.selected_atoms()
        if sel:
            out.append(('Selection', sel))
        frags = doc.fragments()
        if len(frags) == 1:
            out.append(('Structure', None))
        else:
            for k, f in enumerate(frags, 1):
                out.append((f'Structure {k}: {chem.analyse(doc, f).formula}', set(f)))
            if frags:
                out.append(('All structures together', None))
        return out

    def refresh(self, keep_choice=False):
        choices = self._choices()
        if not keep_choice:
            cur = self.which.currentText()
            self.which.blockSignals(True)
            self.which.clear()
            for label, _ids in choices:
                self.which.addItem(label)
            idx = 0
            if self.canvas.selected_atoms():
                idx = 0
            elif cur in [c[0] for c in choices]:
                idx = [c[0] for c in choices].index(cur)
            self.which.setCurrentIndex(idx)
            self.which.blockSignals(False)
        k = self.which.currentIndex()
        ids = choices[k][1] if 0 <= k < len(choices) else None
        info = chem.analyse(self.canvas.doc, ids)
        self.show_info(info)
        self.iupac = chem.iupac_name(self.canvas.doc, ids) if info.atoms and not info.errors else ''
        self.name.setText(self.iupac or ('–' if not info.atoms else '(no name found)'))
        if self.tabs.currentIndex() == 1:                   # properties only while their tab is shown
            rows = chem.properties(self.canvas.doc, ids) if info.atoms else []
            self.props.setRowCount(len(rows))
            for r, (name, value, tip) in enumerate(rows):
                a, b = QtWidgets.QTableWidgetItem(name), QtWidgets.QTableWidgetItem(value)
                a.setToolTip(tip); b.setToolTip(tip)
                self.props.setItem(r, 0, a); self.props.setItem(r, 1, b)

    def show_info(self, info):
        self.info = info
        if not info.atoms:
            self.formula.setText('–')
            for w in (self.exact, self.average, self.charge, self.analysis, self.mz):
                w.setText('–')
            self.table.setRowCount(0); self.pattern.setRowCount(0)
            self.smiles.clear(); self.inchikey.clear(); self.errors.clear()
            return
        self.formula.setText(info.formula_html)
        self.exact.setText(f'{info.exact:.5f}')
        self.average.setText(f'{info.average:.3f}')
        self.charge.setText(chem.charge_text(info.charge) or '0')
        self.analysis.setText(info.analysis)
        self.mz.setText(', '.join(f'{m:.2f} ({p:.1f}%)' for m, p in chem.isotope_pattern(info.counts, 0)[:4]))
        self.table.blockSignals(True)
        self.table.setRowCount(len(info.ions))
        for r, ion in enumerate(info.ions):
            a = QtWidgets.QTableWidgetItem(ion.label)
            b = QtWidgets.QTableWidgetItem(f'{ion.mz:.4f}')
            b.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            if ion.z < 0:
                a.setForeground(QtGui.QColor('#1565C0')); b.setForeground(QtGui.QColor('#1565C0'))
            self.table.setItem(r, 0, a); self.table.setItem(r, 1, b)
        self.table.blockSignals(False)
        row = min(max(self.table.currentRow(), 0), len(info.ions) - 1)
        self.show_pattern(row)
        self.smiles.setText(info.smiles); self.smiles.setCursorPosition(0)
        self.inchikey.setText(info.inchikey); self.inchikey.setCursorPosition(0)
        self.errors.setText('\n'.join(info.errors))

    def show_pattern(self, row):
        """The isotope peaks of ion row (M, M+1, M+2 ... with their relative heights)."""
        info = getattr(self, 'info', None)
        if info is None or not (0 <= row < len(info.ions)):
            self.pattern.setRowCount(0)
            return
        ion = info.ions[row]
        peaks = chem.isotope_pattern(ion.counts, ion.z)
        self.pattern_title.setText(f'Isotope pattern of {ion.label}')
        self.pattern.setRowCount(len(peaks))
        for r, (mz, pct) in enumerate(peaks):
            for c, text in enumerate((f'{mz:.4f}', f'{pct:.1f}')):
                it = QtWidgets.QTableWidgetItem(text)
                it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.pattern.setItem(r, c, it)

    def copy_rows(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}) or range(self.table.rowCount())
        text = '\n'.join(f'{self.table.item(r, 0).text()}\t{self.table.item(r, 1).text()}' for r in rows)
        QtWidgets.QApplication.clipboard().setText(text)
