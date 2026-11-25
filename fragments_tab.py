"""
fragments_tab.py - the Fragments tab of the formula panel: the MS fragments the structure can logically give
(fragments.predict), for EI (M•+), ESI+ ([M+H]+) or ESI− ([M−H]−), listed by m/z with the ion, what is lost, the rule
and how likely it is. Picking a row shows on the drawing which atoms stay (green), which leave (red) and the bonds that
break; Draw puts the fragment ion under the structure, after an arrow labelled with the loss.
"""
from PyQt6 import QtCore, QtGui, QtWidgets

import fragments

SHOW = (('Likely', fragments.LIKELY), ('Likely and possible', fragments.POSSIBLE), ('All', 0.0))
DOTS = {'likely': '●●●', 'possible': '●●○', 'minor': '●○○'}


class _Num(QtWidgets.QTableWidgetItem):
    """A cell that sorts by its number."""

    def __lt__(self, other):
        try:
            return float(self.data(QtCore.Qt.ItemDataRole.UserRole)) < float(other.data(QtCore.Qt.ItemDataRole.UserRole))
        except (TypeError, ValueError):
            return super().__lt__(other)


class FragmentsTab(QtWidgets.QWidget):
    def __init__(self, canvas, parent=None):
        super().__init__(parent)
        self.canvas = canvas
        self.frags = []
        self.pre = None
        self._key = None
        self._ids = None
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 6, 4, 4)
        row = QtWidgets.QHBoxLayout()
        self.mode = QtWidgets.QComboBox()
        for key, text in fragments.MODES:
            self.mode.addItem(text, key)
        self.mode.setToolTip('ESI+: collision-induced fragments of [M+H]+ (LC-MS/MS). ESI−: of [M−H]−. '
                             'EI: of the radical cation M•+ (GC-MS, 70 eV)')
        self.steps = QtWidgets.QSpinBox()
        self.steps.setRange(1, 3)
        self.steps.setValue(2)
        self.steps.setPrefix('Steps: ')
        self.steps.setToolTip('1: fragments of the precursor ion only. 2, 3: those fragments lose small molecules '
                              'again (H2O, CO, NH3 ...)')
        self.level = QtWidgets.QComboBox()
        for text, _v in SHOW:
            self.level.addItem(text)
        self.level.setCurrentIndex(1)
        self.level.setToolTip('Likelihood is a rough plausibility from the fragmentation rule and the ion it makes, '
                             'not a predicted intensity')
        row.addWidget(self.mode)
        row.addWidget(self.steps)
        row.addWidget(self.level, 1)
        lay.addLayout(row)
        self.head = QtWidgets.QLabel()
        self.head.setWordWrap(True)
        self.head.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.head)
        self.table = QtWidgets.QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(['m/z', 'Ion', 'Lost', 'How', ''])
        self.table.verticalHeader().hide()
        hh = self.table.horizontalHeader()
        for c in (0, 1, 2, 4):
            hh.setSectionResizeMode(c, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(True)
        self.table.setWordWrap(False)
        self.table.setToolTip('Pick a row to see it on the drawing; double-click to draw the fragment. '
                              'Ctrl+C copies the list.')
        self.table.itemSelectionChanged.connect(self._picked)
        self.table.cellDoubleClicked.connect(lambda _r, _c: self.draw())
        QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Copy, self.table, self.copy,
                        context=QtCore.Qt.ShortcutContext.WidgetShortcut)
        lay.addWidget(self.table, 1)
        brow = QtWidgets.QHBoxLayout()
        self.note = QtWidgets.QLabel()
        self.note.setWordWrap(True)
        self.note.setStyleSheet('color: palette(mid);')
        brow.addWidget(self.note, 1)
        copy = QtWidgets.QPushButton('Copy')
        copy.setToolTip('The list as tab-separated text (m/z, ion, lost, how, likelihood)')
        copy.clicked.connect(self.copy)
        self.draw_b = QtWidgets.QPushButton('Draw')
        self.draw_b.setToolTip('The picked fragment ion under the structure, after an arrow labelled with the loss')
        self.draw_b.clicked.connect(self.draw)
        brow.addWidget(copy)
        brow.addWidget(self.draw_b)
        lay.addLayout(brow)
        for w in (self.mode, self.level):
            w.currentIndexChanged.connect(lambda _i: self.recompute())
        self.steps.valueChanged.connect(lambda _v: self.recompute())

    # ------------------------------------------------------------------------------------------- filling
    def message(self, text):
        self._key = None
        self.frags, self.pre = [], None
        self.head.setText(text)
        self.table.setRowCount(0)
        self.note.clear()
        self.canvas.show_fragment(None)

    def refresh(self, ids, info):
        """The fragments of the structure with atoms ids (None: the whole drawing); info is its chem.analyse."""
        if not info.atoms:
            self.message('Draw a structure to see the fragments it can give.')
            return
        if info.fragments > 1:
            self.message('Pick one structure in the list above (or select it): fragments are worked out for one '
                         'molecule at a time.')
            return
        if info.errors:
            self.message('Fix the structure first: ' + info.errors[0])
            return
        self._ids = ids
        key = (info.smiles, tuple(sorted(ids)) if ids else tuple(sorted(self.canvas.doc.atoms)),
               self.mode.currentData(), self.steps.value(), self.level.currentIndex())
        if key == self._key:
            return
        self._key = key
        self._fill()

    def recompute(self):
        self._key = None
        if self.isVisible():
            self._fill()

    def _fill(self):
        doc = self.canvas.doc
        ids = self._ids
        pre, frags = fragments.predict(doc, ids, self.mode.currentData(), self.steps.value())
        if pre is None:
            self.message('This structure could not be read.')
            return
        self.pre = pre
        least = SHOW[self.level.currentIndex()][1]
        shown = [f for f in frags if f.score >= least]
        self.frags = shown
        if getattr(pre, 'drawn_ion', False):
            self.head.setText(f'The ion as drawn, {pre.ion}  m/z {pre.mz:.4f}: {len(shown)} fragments '
                              f'(of {len(frags)}). Draw the neutral molecule for EI / ESI ions.')
        else:
            self.head.setText(f'Precursor {pre.label}  {pre.ion}  m/z {pre.mz:.4f}: {len(shown)} fragments '
                              f'(of {len(frags)})')
        self.mode.setEnabled(not getattr(pre, 'drawn_ion', False))
        t = self.table
        t.setSortingEnabled(False)
        t.blockSignals(True)
        t.setRowCount(len(shown))
        right = QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter
        for r, f in enumerate(shown):
            mz = _Num(f'{f.mz:.4f}')
            mz.setData(QtCore.Qt.ItemDataRole.UserRole, f.mz)
            mz.setTextAlignment(right)
            mz.setData(QtCore.Qt.ItemDataRole.UserRole + 1, r)
            like = _Num(DOTS[f.likelihood])
            like.setData(QtCore.Qt.ItemDataRole.UserRole, f.score)
            like.setToolTip(f'{f.likelihood} ({f.score:.2f})')
            how = QtWidgets.QTableWidgetItem(f.how)
            tip = f.how + (('\nAlso: ' + '\n      '.join(f.also)) if f.also else '')
            if f.S is None:
                tip += '\n(worked out from the formula: no structure to draw)'
            how.setToolTip(tip)
            cells = (mz, QtWidgets.QTableWidgetItem(f.ion), QtWidgets.QTableWidgetItem('−' + f.lost), how, like)
            for c, it in enumerate(cells):
                if f.likelihood == 'minor':
                    it.setForeground(QtGui.QColor('#8A8A8A'))
                t.setItem(r, c, it)
        t.blockSignals(False)
        t.setSortingEnabled(True)
        if t.horizontalHeader().sortIndicatorSection() not in range(5) or not hasattr(self, '_sorted'):
            t.sortItems(0, QtCore.Qt.SortOrder.DescendingOrder)
            self._sorted = True
        self.note.setText('Green: atoms the ion keeps; red: atoms lost; red bars: bonds broken.' if shown else
                          'No fragment at this likelihood: choose "All".')
        self.canvas.show_fragment(None)

    # ------------------------------------------------------------------------------------------- using a row
    def current(self):
        rows = {i.row() for i in self.table.selectedItems()}
        if not rows:
            return None
        k = self.table.item(rows.pop(), 0).data(QtCore.Qt.ItemDataRole.UserRole + 1)
        return self.frags[k] if k is not None and 0 <= k < len(self.frags) else None

    def _picked(self):
        f = self.current()
        if f is None:
            self.canvas.show_fragment(None)
            return
        doc = self.canvas.doc
        keep, lost = set(f.atoms) & set(doc.atoms), set(f.lost_atoms) & set(doc.atoms)
        cut = {i for i, b in doc.bonds.items() if (b.a in keep and b.b in lost) or (b.b in keep and b.a in lost)}
        self.canvas.show_fragment((keep, lost, cut) if f.S is not None else (keep, set(), set()))

    def draw(self):
        f = self.current()
        if f is None:
            self.canvas.message.emit('Pick a fragment in the list first.')
            return
        piece, why = fragments.draw_doc(self.canvas.doc, f)
        if piece is None:
            self.canvas.message.emit(why)
            return
        parent = set(self.pre.atoms) & set(self.canvas.doc.atoms)
        self.canvas.draw_fragment(piece, parent, '−' + f.lost)

    def copy(self):
        t = self.table
        lines = ['m/z\tIon\tLost\tHow\tLikelihood']
        for r in range(t.rowCount()):
            k = t.item(r, 0).data(QtCore.Qt.ItemDataRole.UserRole + 1)
            f = self.frags[k]
            lines.append(f'{f.mz:.4f}\t{f.ion}\t−{f.lost}\t{f.how}\t{f.likelihood}')
        QtWidgets.QApplication.clipboard().setText('\n'.join(lines))
