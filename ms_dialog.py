"""Structure > Match m/z: a measured m/z against the structures drawn (their ions and adducts, fragments drawn as
ions) and against elemental compositions."""
from PyQt6 import QtCore, QtWidgets

import chem
import ms


class MatchDialog(QtWidgets.QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.setWindowTitle('Match m/z')
        self.canvas = window.canvas
        self.resize(760, 520)
        lay = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QHBoxLayout()
        self.mz = QtWidgets.QLineEdit()
        self.mz.setPlaceholderText('measured m/z, e.g. 121.0648')
        self.ppm = QtWidgets.QDoubleSpinBox(decimals=1, minimum=0.1, maximum=500.0, value=5.0, suffix=' ppm')
        self.z = QtWidgets.QComboBox()
        for text, z in (('+1', 1), ('−1', -1), ('+2', 2), ('−2', -2), ('+3', 3)):
            self.z.addItem(text, z)
        self.elements = QtWidgets.QLineEdit('C H N O S')
        self.elements.setToolTip('Elements for the compositions; a number caps one: C H N4 O S P Cl Br F Na K')
        go = QtWidgets.QPushButton('Search')
        go.setDefault(True)
        go.clicked.connect(self.search)
        self.mz.returnPressed.connect(self.search)
        for w in (QtWidgets.QLabel('m/z'), self.mz, self.ppm, QtWidgets.QLabel('charge'), self.z,
                  QtWidgets.QLabel('elements'), self.elements, go):
            form.addWidget(w)
        form.setStretch(1, 2); form.setStretch(6, 2)
        lay.addLayout(form)
        lay.addWidget(QtWidgets.QLabel('Drawn structures and fragments (double-click selects it):'))
        self.drawn = self._table(('Structure', 'Ion', 'Formula', 'm/z calc.', 'ppm'))
        self.drawn.cellDoubleClicked.connect(self._select)
        lay.addWidget(self.drawn, 1)
        lay.addWidget(QtWidgets.QLabel('Elemental compositions of the ion (RDBE: rings + double bonds; even / odd '
                                       'electron):'))
        self.comps = self._table(('Ion formula', 'm/z calc.', 'ppm', 'RDBE', 'Electrons'))
        lay.addWidget(self.comps, 2)
        self.note = QtWidgets.QLabel('')
        lay.addWidget(self.note)
        self._ids = []

    @staticmethod
    def _table(heads):
        t = QtWidgets.QTableWidget(0, len(heads))
        t.setHorizontalHeaderLabels(heads)
        t.verticalHeader().setVisible(False)
        t.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        t.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        t.horizontalHeader().setStretchLastSection(True)
        return t

    @staticmethod
    def _fill(table, rows):
        table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            for k, v in enumerate(row):
                it = QtWidgets.QTableWidgetItem(v)
                if k and v[:1] in '-+0123456789':
                    it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
                table.setItem(r, k, it)
        table.resizeColumnsToContents()

    def units(self):
        doc = self.canvas.doc
        out = []
        for k, u in enumerate(self.canvas.analysis_units(set(doc.atoms)), 1):
            out.append((f'{k}: {chem.analyse(doc, set(u)).formula}', set(u)))
        return out

    def search(self):
        try:
            mz = float(self.mz.text().replace(',', '.'))
        except ValueError:
            self.note.setText('Type the measured m/z (a number).')
            return
        tol, z = self.ppm.value(), self.z.currentData()
        drawn = [r for r in ms.drawn_matches(self.canvas.doc, self.units(), mz, tol)
                 if r[1].endswith('+') == (z > 0)]
        self._ids = [r[5] for r in drawn]
        self._fill(self.drawn, [(r[0], r[1], r[2], f'{r[3]:.4f}', f'{r[4]:+.1f}') for r in drawn])
        limits = ms.limits_from_text(self.elements.text()) or ms.DEFAULT_LIMITS
        comps = ms.compositions(mz, z, tol, limits)
        sign = chem.charge_text(z)
        self._fill(self.comps, [(f + sign, f'{calc:.4f}', f'{d:+.2f}', f'{rdbe:g}', e) for f, calc, d, rdbe, e in comps])
        self.note.setText(f'{len(drawn)} drawn match(es), {len(comps)} composition(s) within {tol:g} ppm.')

    def _select(self, row, _col):
        if 0 <= row < len(self._ids):
            self.canvas.set_selection([('atom', a) for a in self._ids[row]])
