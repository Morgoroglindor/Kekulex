"""
nmr_dialog.py - the NMR prediction window: 1H or 13C, the solvent, a simulated spectrum, the signals (shift +/-
uncertainty, multiplicity, hydrogens, atoms) and the experimental-section line to copy. Clicking a signal selects its
atoms on the drawing; Label atoms writes the shifts beside them. The prediction runs in the background (the first one
loads the model, a few seconds).
"""
from PyQt6 import QtCore, QtGui, QtWidgets

import nmr
import render
import theme
from plotcanvas import PlotCanvas

EXCHANGEABLE_NOTE = 'OH / NH / SH: exchangeable, the predicted shift is unreliable'
EXCHANGED_NOTE = 'OH / NH / SH: exchanges with D in {}, not seen (left out of the spectrum and the text)'


class _Worker(QtCore.QThread):
    done = QtCore.pyqtSignal(object, str)                   # (Result or None, error text)

    def __init__(self, doc, atoms, nucleus, solvent):
        super().__init__()
        self.args = (doc, atoms, nucleus, solvent)

    def run(self):
        try:
            self.done.emit(nmr.predict(*self.args), '')
        except Exception as e:                              # shown in the window, also logged
            import logging
            logging.getLogger('kekulex').warning('NMR prediction failed', exc_info=True)
            self.done.emit(None, str(e))


class SpectrumView(PlotCanvas):
    """The simulated spectrum, ppm decreasing to the right as on a spectrometer; the selected signal in the accent.
    The mouse as in Spectra viewer (plotcanvas.py), but drag and right drag act on ppm only: drag a range to zoom, wheel
    zooms ppm, Ctrl+wheel the height, right drag pans, double-click or F shows the whole spectrum; a click picks the
    nearest signal."""
    x_only = True

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(160)
        self.result = None
        self.curve = None
        self.highlight = None
        self.mhz = 400
        self.on_reset = self.redraw

    def show_result(self, result, mhz=400, keep_zoom=False):
        self.result, self.mhz = result, mhz
        self.curve = nmr.spectrum(result, mhz) if result else None
        if not keep_zoom:
            self.view = None
        self.redraw()

    def redraw(self):
        ax = self.ax
        ax.clear()
        if self.result and self.curve is not None:
            x, y = self.curve
            ymax = float(y.max()) or 1.0
            ax.plot(x, y, lw=0.9, color='#1F3B73')
            ax.set_xlim(float(x.max()), float(x.min()))     # ppm decreasing to the right
            ax.set_ylim(-0.03 * ymax, 1.15 * ymax)
            h = self.result.nucleus == '1H'
            for k, s in enumerate(self.result.signals):     # each shift above the tallest line of its signal
                if nmr.exchanged(self.result, s):
                    continue
                pos =[p for p, _w in nmr.lines(s, self.mhz)] if h else [s.shift]
                span = (x >= min(pos) - 0.01) & (x <= max(pos) + 0.01)
                top = float(y[span].max()) if span.any() else 0.0
                ax.annotate(f'{s.shift:.2f}' if h else f'{s.shift:.1f}', (s.shift, top), xytext=(0, 3),
                            textcoords='offset points', ha='center', va='bottom', fontsize=7, rotation=90,
                            color=theme.current_accent if k == self.highlight else theme.T['fig_fg'])
            if self.highlight is not None and self.highlight < len(self.result.signals):
                ax.axvline(self.result.signals[self.highlight].shift, color=theme.current_accent, lw=0.8, ls='--')
            ax.set_xlabel(f'δ (ppm)   {self.mhz if h else self.mhz // 4} MHz, predicted')
        ax.set_yticks([])
        self.style()
        self.apply_view()
        self.draw_idle()


class NMRDialog(QtWidgets.QDialog):
    def __init__(self, canvas, atom_ids, nucleus='1H', parent=None):
        super().__init__(parent)
        self.canvas = canvas
        self.atoms = set(atom_ids) if atom_ids else None
        self.result = None
        self.worker = None
        self.setWindowTitle('NMR prediction')
        theme.register_window(self)
        lay = QtWidgets.QVBoxLayout(self)
        row = QtWidgets.QHBoxLayout()
        self.nucleus = QtWidgets.QComboBox(); self.nucleus.addItems(['¹H', '¹³C'])
        self.nucleus.setCurrentIndex(0 if nucleus == '1H' else 1)
        self.solvent = QtWidgets.QComboBox()
        for name, key in nmr.SOLVENTS:
            self.solvent.addItem(name, key)
        row.addWidget(QtWidgets.QLabel('Nucleus:')); row.addWidget(self.nucleus)
        row.addSpacing(12)
        row.addWidget(QtWidgets.QLabel('Solvent:')); row.addWidget(self.solvent)
        row.addSpacing(12)
        self.mhz = QtWidgets.QComboBox()
        for f in (300, 400, 500, 600, 800):
            self.mhz.addItem(f'{f} MHz', f)
        self.mhz.setCurrentIndex(1)
        self.mhz.setToolTip('Spectrometer frequency (¹H): spreads the multiplets in ppm and goes into the text')
        row.addWidget(QtWidgets.QLabel('Spectrometer:')); row.addWidget(self.mhz)
        row.addStretch(1)
        self.status = QtWidgets.QLabel()
        row.addWidget(self.status)
        lay.addLayout(row)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)  # drag the handle to give the spectrum more room
        top = QtWidgets.QWidget()
        tl = QtWidgets.QVBoxLayout(top); tl.setContentsMargins(0, 0, 0, 0); tl.setSpacing(2)
        self.view = SpectrumView()
        self.view.on_click = self._clicked_ppm
        tl.addWidget(self.view, 1)
        hint = QtWidgets.QLabel('Drag: zoom to the box · wheel: zoom ppm · Ctrl+wheel: zoom height · right drag: pan · '
                                'double-click or F: whole spectrum · click: pick the nearest signal')
        hint.setStyleSheet('color: #66707A;')
        tl.addWidget(hint)
        split.addWidget(top)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(['δ (ppm)', '±', 'Multiplicity', 'J (Hz)', 'Integral', 'Note'])
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.currentCellChanged.connect(lambda r, _c, _pr, _pc: self._signal_picked(r))
        split.addWidget(self.table)
        split.setStretchFactor(0, 3); split.setStretchFactor(1, 2)
        lay.addWidget(split, 1)
        self.text = QtWidgets.QLineEdit(); self.text.setReadOnly(True)
        lay.addWidget(self.text)
        buttons = QtWidgets.QHBoxLayout()
        for label, slot in (('Copy text', self.copy_text), ('Copy table', self.copy_table),
                            ('Paste text under structure', self.paste_text), ('Label atoms', self.label_atoms),
                            ('Save spectrum picture…', self.save_picture)):
            b = QtWidgets.QPushButton(label); b.clicked.connect(slot); buttons.addWidget(b)
        buttons.addStretch(1)
        close = QtWidgets.QPushButton('Close'); close.clicked.connect(self.close); buttons.addWidget(close)
        lay.addLayout(buttons)
        note = QtWidgets.QLabel('Shifts predicted by CSP5 (machine learning on measured spectra). Multiplicities and J '
                                'are first-order estimates from the structure (typical values; Karplus on a 3D '
                                'conformer for ring CH-CH). Check predictions against your data.')
        note.setWordWrap(True); note.setStyleSheet('color: #66707A;')
        lay.addWidget(note)
        self.nucleus.currentIndexChanged.connect(lambda _i: self.run())
        self.solvent.currentIndexChanged.connect(lambda _i: self.run())
        self.mhz.currentIndexChanged.connect(lambda _i: self._mhz_changed())
        self.resize(900, 780)
        self.run()

    def nucleus_key(self):
        return '1H' if self.nucleus.currentIndex() == 0 else '13C'

    def run(self):
        if self.worker is not None and self.worker.isRunning():
            self._pending = True
            return
        self._pending = False
        self.status.setText('Predicting… (the first time loads the model)')
        self.worker = _Worker(self.canvas.doc.copy(), self.atoms, self.nucleus_key(), self.solvent.currentData())
        self.worker.done.connect(self._done)
        self.worker.start()

    def _done(self, result, error):
        if getattr(self, '_pending', False):
            self.run()
            return
        if result is None:
            self.status.setText('')
            self.result = None
            self.view.show_result(None)
            self.table.setRowCount(0)
            self.text.setText('')
            QtWidgets.QMessageBox.warning(self, 'NMR prediction', error)
            return
        self.result = result
        self.status.setText(f'{len(result.signals)} signals')
        self.view.highlight = None
        self.view.show_result(result, self.mhz.currentData())
        self.table.setRowCount(len(result.signals))
        h = result.nucleus == '1H'
        for r, s in enumerate(result.signals):
            exch = h and s.multiplicity == 'br s'
            note = (EXCHANGED_NOTE.format(nmr.SOLVENT_NAMES.get(result.solvent, result.solvent))
                    if nmr.exchanged(result, s) else EXCHANGEABLE_NOTE if exch else '')
            cells = (f'{s.shift:.2f}' if h else f'{s.shift:.1f}', f'{s.std:.2f}' if s.std else '',
                     s.multiplicity if h else ('' if len(s.atoms) == 1 else f'{len(s.atoms)} equivalent C'),
                     s.j_text() if h else '', f'{s.n_h}H' if h else f'{len(s.atoms)}C', note)
            for c, text in enumerate(cells):
                it = QtWidgets.QTableWidgetItem(text)
                if c < 5:
                    it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
                if exch:
                    it.setForeground(QtGui.QColor('#9A6700'))
                self.table.setItem(r, c, it)
        self.table.resizeColumnsToContents()
        self.text.setText(nmr.journal_text(result, self.mhz.currentData())); self.text.setCursorPosition(0)

    def _mhz_changed(self):
        if self.result:
            self.view.show_result(self.result, self.mhz.currentData(), keep_zoom=True)
            self.text.setText(nmr.journal_text(self.result, self.mhz.currentData())); self.text.setCursorPosition(0)

    def _clicked_ppm(self, ppm):
        if self.result and self.result.signals and ppm is not None:
            row = min(range(len(self.result.signals)), key=lambda k: abs(self.result.signals[k].shift - ppm))
            self.table.selectRow(row)

    def _signal_picked(self, row):
        if not self.result or not (0 <= row < len(self.result.signals)):
            return
        self.view.highlight = row
        self.view.redraw()
        atoms = [a for a in self.result.signals[row].atoms if a in self.canvas.doc.atoms]
        self.canvas.set_selection([('atom', a) for a in atoms])

    # ------------------------------------------------------------------------------------------- output
    def copy_text(self):
        QtWidgets.QApplication.clipboard().setText(self.text.text())

    def copy_table(self):
        rows = []
        for r in range(self.table.rowCount()):
            rows.append('\t'.join(self.table.item(r, c).text() for c in range(self.table.columnCount())))
        QtWidgets.QApplication.clipboard().setText('\n'.join(rows))

    def paste_text(self):
        if self.result:
            atoms = self.atoms or set(self.canvas.doc.atoms)
            self.canvas.add_caption_below(atoms, self.text.text())

    def label_atoms(self):
        """Each predicted shift as a small caption beside its atom (one undo step)."""
        if not self.result:
            return
        c = self.canvas
        c.snapshot()
        BL = render.style_of(c.doc)['bond_length']
        size = render.style_of(c.doc)['caption_size'] * 0.75
        for s in self.result.signals:
            for a in s.atoms:
                if a not in c.doc.atoms:
                    continue
                at = c.doc.atoms[a]
                fx, fy = c.doc.free_direction(a)
                text = f'{s.shift:.2f}' if self.result.nucleus == '1H' else f'{s.shift:.1f}'
                c.doc.add_text(at.x + fx * 0.75, at.y + fy * 0.75 - size / BL / 2, text, align='center', size=size,
                               color='#1F3B73')
        c.edited()

    def save_picture(self):
        """The spectrum as shown (with the zoom), PNG at 300 dpi, SVG or PDF."""
        path, _f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save spectrum picture', 'nmr_prediction.png',
                                                         'PNG picture (*.png);;SVG (*.svg);;PDF (*.pdf)')
        if path:
            if not path.lower().endswith(('.png', '.svg', '.pdf')):
                path += '.png'
            self.view.fig.savefig(path, dpi=300)

    def closeEvent(self, e):
        if self.worker is not None and self.worker.isRunning():
            self.worker.wait(30000)
        super().closeEvent(e)
