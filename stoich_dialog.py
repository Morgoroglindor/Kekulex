"""
stoich_dialog.py - typing into a stoichiometry grid: two tables (reactants, products) with a column per compound.
White cells take a value with or without its unit ('250 mg', '1.5 mmol', '2 mL', '0.5 M', '95 %'); grey cells are
calculated. Right-click a row's name to hide or show it on the drawing. The drawing follows as you type.
"""
from PyQt6 import QtCore, QtGui, QtWidgets

import stoich
import theme

UNIT_HINT = {'mass': ' (g)', 'moles': ' (mol)', 'volume': ' (mL)', 'molarity': ' (M)', 'density': ' (g/mL)',
             'percent': ' (%)'}


class GridDialog(QtWidgets.QDialog):
    def __init__(self, canvas, grid_id, parent=None):
        super().__init__(parent)
        self.canvas = canvas
        self.gid = grid_id
        self.setWindowTitle('Stoichiometry')
        theme.register_window(self)
        canvas.snapshot()                                   # one undo step for the whole edit
        lay = QtWidgets.QVBoxLayout(self)
        tip = QtWidgets.QLabel('Type values with or without units (250 mg, 1.5 mmol, 2 mL, 0.5 M, 95 %). '
                               'Bold on the drawing = typed, normal = calculated. Right-click a row name to hide it.')
        tip.setWordWrap(True)
        lay.addWidget(tip)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Limiting reactant:'))
        self.limiting = QtWidgets.QComboBox()
        self.limiting.currentIndexChanged.connect(self._limiting_changed)
        row.addWidget(self.limiting, 1)
        lay.addLayout(row)
        self.tables = {}
        for side, title in (('reactants', 'Reactants'), ('products', 'Products')):
            lay.addWidget(QtWidgets.QLabel(f'<b>{title}</b>'))
            t = QtWidgets.QTableWidget()
            t.verticalHeader().setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
            t.verticalHeader().customContextMenuRequested.connect(lambda pos, t=t, s=side: self._row_menu(t, s, pos))
            t.itemChanged.connect(lambda item, s=side: self._edited(s, item))
            lay.addWidget(t)
            self.tables[side] = t
        self.error = QtWidgets.QLabel()
        self.error.setStyleSheet('color: #C62828;')
        lay.addWidget(self.error)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                         QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.accepted.connect(self.accept); box.rejected.connect(self.reject)
        lay.addWidget(box)
        self._filling = False
        self.fill()
        self.resize(820, 760)

    @property
    def grid(self):
        return self.canvas.doc.grids[self.gid]

    def fill(self):
        """Show what the grid holds and what follows from it."""
        self._filling = True
        try:
            from PyQt6.QtCore import QLocale
            dec = QLocale.system().decimalPoint() or ','
            table = stoich.compute(self.canvas.doc, self.grid)
            self.limiting.blockSignals(True)
            self.limiting.clear()
            for k, col in enumerate(table['reactants']):
                self.limiting.addItem(f'{k + 1}. {col.formula}')
                if col.cells['limiting'].value:
                    self.limiting.setCurrentIndex(k)
            self.limiting.blockSignals(False)
            for side, rows, inputs in (('reactants', stoich.REACTANT_ROWS, stoich.REACTANT_INPUTS),
                                       ('products', stoich.PRODUCT_ROWS, stoich.PRODUCT_INPUTS)):
                t = self.tables[side]
                cols = table[side]
                rows = [(k, label) for k, label in rows if k not in ('formula', 'limiting')]
                t.setColumnCount(len(cols)); t.setRowCount(len(rows))
                t.setHorizontalHeaderLabels([c.formula for c in cols])
                t.setVerticalHeaderLabels([label + UNIT_HINT.get(stoich.KIND.get(k, ''), '') +
                                           ('  (hidden)' if k in self.grid.hidden else '') for k, label in rows])
                t.rows = [k for k, _l in rows]
                for c, col in enumerate(cols):
                    for r, (key, _label) in enumerate(rows):
                        cell = col.cells.get(key)
                        num, unit = stoich.display(key, cell, dec) if cell else ('', '')
                        item = QtWidgets.QTableWidgetItem((num + (' ' + unit if unit else '')).strip())
                        if key in inputs:
                            f = item.font(); f.setBold(bool(cell and cell.typed)); item.setFont(f)
                            if not (cell and cell.typed):
                                item.setForeground(QtGui.QColor('#808080'))
                        else:
                            item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                            item.setBackground(QtGui.QColor(theme.T['alt']))
                        t.setItem(r, c, item)
                t.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
                t.verticalHeader().setDefaultSectionSize(24)
                t.setMinimumHeight(24 * (len(rows) + 1) + 6)
        finally:
            self._filling = False

    def _columns(self, side):
        doc = self.canvas.doc
        return [c for c in getattr(self.grid, side) if any(a in doc.atoms for a in c.get('atoms', []))]

    def _edited(self, side, item):
        if self._filling:
            return
        t = self.tables[side]
        key = t.rows[item.row()]
        col = self._columns(side)[item.column()]
        text = item.text().strip()
        try:
            value = stoich.parse(text, stoich.KIND.get(key, 'number'))
        except ValueError as e:
            self.error.setText(str(e))
            return
        self.error.clear()
        if value is None:                                   # emptied: back to calculated
            col.pop(key, None)
        else:
            col[key] = value
        self._apply()

    def _limiting_changed(self, k):
        for j, col in enumerate(self._columns('reactants')):
            col['limiting'] = j == k
        self._apply()

    def _row_menu(self, t, side, pos):
        r = t.verticalHeader().logicalIndexAt(pos)
        if r < 0:
            return
        key = t.rows[r]
        menu = QtWidgets.QMenu(self)
        hide = menu.addAction('Show this row' if key in self.grid.hidden else 'Hide this row')
        show_all = menu.addAction('Show all rows')
        clear = menu.addAction('Clear the typed values of this row')
        chosen = menu.exec(t.verticalHeader().mapToGlobal(pos))
        if chosen is hide:
            if key in self.grid.hidden:
                self.grid.hidden.remove(key)
            else:
                self.grid.hidden.append(key)
        elif chosen is show_all:
            self.grid.hidden.clear()
        elif chosen is clear:
            for col in getattr(self.grid, side):
                col.pop(key, None)
        else:
            return
        self._apply()

    def _apply(self):
        self.canvas._renderer = None
        self.canvas.update()
        self.fill()

    def accept(self):
        self.canvas.edited()
        super().accept()

    def reject(self):
        self.canvas.doc = self.canvas.undo_stack.pop()       # back to how it was before the dialog
        self.canvas.edited()
        super().reject()
