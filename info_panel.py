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
        self.formula = QtWidgets.QLabel('–')
        f = self.formula.font(); f.setPointSizeF(f.pointSizeF() * 1.6); self.formula.setFont(f)
        self.formula.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.formula)
        form = QtWidgets.QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        self.exact = self._value(); self.average = self._value(); self.charge = self._value()
        form.addRow('Exact mass:', self.exact)
        form.addRow('Average MW:', self.average)
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
        lay.addWidget(self.table, 1)
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
        lay.addWidget(self.errors)
        self.timer = QtCore.QTimer(self, singleShot=True, interval=120)
        self.timer.timeout.connect(self.refresh)
        canvas.changed.connect(self.timer.start)
        canvas.selection_changed.connect(self.timer.start)
        self.refresh()

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

    def show_info(self, info):
        self.info = info
        if not info.atoms:
            self.formula.setText('–')
            for w in (self.exact, self.average, self.charge):
                w.setText('–')
            self.table.setRowCount(0)
            self.smiles.clear(); self.inchikey.clear(); self.errors.clear()
            return
        self.formula.setText(info.formula_html)
        self.exact.setText(f'{info.exact:.5f}')
        self.average.setText(f'{info.average:.3f}')
        self.charge.setText(chem.charge_text(info.charge) or '0')
        self.table.setRowCount(len(info.ions))
        for r, (label, z, mz) in enumerate(info.ions):
            a = QtWidgets.QTableWidgetItem(label)
            b = QtWidgets.QTableWidgetItem(f'{mz:.4f}')
            b.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
            if z < 0:
                a.setForeground(QtGui.QColor('#1565C0')); b.setForeground(QtGui.QColor('#1565C0'))
            self.table.setItem(r, 0, a); self.table.setItem(r, 1, b)
        self.smiles.setText(info.smiles); self.smiles.setCursorPosition(0)
        self.inchikey.setText(info.inchikey); self.inchikey.setCursorPosition(0)
        self.errors.setText('\n'.join(info.errors))

    def copy_rows(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}) or range(self.table.rowCount())
        text = '\n'.join(f'{self.table.item(r, 0).text()}\t{self.table.item(r, 1).text()}' for r in rows)
        QtWidgets.QApplication.clipboard().setText(text)
