"""
spectrum_dialog.py - Structure > Annotate mass spectrum: a measured spectrum (CSV / TXT / JCAMP-DX / MSP / MGF / mzML,
or pasted) as a stick plot, its peaks matched to the ions of the structures drawn (their adducts, M•+, the rule-based
fragments of fragments.py and the pieces of breaking one bond), labelled on the plot and listed. Insert puts the
annotated spectrum on the drawing as a picture.
"""
import io
import re

from PyQt6 import QtCore, QtGui, QtWidgets

import chem
import ms
import msfile
import theme
from plotcanvas import PlotCanvas

FILE_FILTER = ('Mass spectra (*.csv *.txt *.tsv *.xy *.jdx *.dx *.msp *.mgf *.mzML *.mzml);;All files (*)')
MATCH_COLOR = '#D81B60'


def ion_text(match):
    """The ion's formula with its charge and radical dot: C2H3O+, C4H8O•+, C4H9O+ for [M+H]+."""
    m = re.search(r'(•?)(\d*[+−-])$', match[1])
    return match[2] + (m.group(1) + m.group(2) if m else '')


def draw_spectrum(ax, rows, fg, stick, label_color, most_labels=14):
    """Sticks of the peaks (relative intensity, %), the matched ones coloured and labelled with their ion formula;
    labels that would overlap a taller labelled peak's are left out."""
    if not rows:
        return
    top = max(r[1] for r in rows) or 1.0
    xs = [r[0] for r in rows]
    span = (max(xs) - min(xs)) or 1.0
    for mz, y, match, _err in rows:
        ax.vlines(mz, 0, 100 * y / top, color=label_color if match else stick, lw=1.1 if match else 0.8)
    placed = []
    for mz, y, match, _err in sorted(rows, key=lambda r: -r[1]):
        if not match or len(placed) >= most_labels:
            continue
        if any(abs(mz - p) < span * 0.035 for p in placed):
            continue
        placed.append(mz)
        ax.annotate(f'{mz:.4f}\n{ion_text(match)}' if mz != round(mz) else f'{mz:g}\n{ion_text(match)}',
                    (mz, 100 * y / top), xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=7.5,
                    color=label_color)
    ax.set_xlabel('m/z', color=fg)
    ax.set_ylabel('Relative intensity (%)', color=fg)
    ax.set_ylim(0, 118)
    pad = span * 0.05
    ax.set_xlim(min(xs) - pad, max(xs) + pad)


class SpectrumDialog(QtWidgets.QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.setWindowTitle('Annotate mass spectrum')
        self.window_ = window
        self.canvas = window.canvas
        self.spectra = []
        self.rows = []
        self.resize(980, 720)
        lay = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QHBoxLayout()
        open_b = QtWidgets.QPushButton('Open spectrum…')
        open_b.clicked.connect(self.open_file)
        paste_b = QtWidgets.QPushButton('Paste')
        paste_b.setToolTip('m/z and intensity pairs from the clipboard (two columns, or JCAMP-DX / MSP text)')
        paste_b.clicked.connect(self.paste)
        self.which = QtWidgets.QComboBox()
        self.which.setMinimumWidth(260)
        self.which.currentIndexChanged.connect(self.update_view)
        self.mode = QtWidgets.QComboBox()
        self.mode.addItem('Positive ions', True)
        self.mode.addItem('Negative ions', False)
        self.tol = QtWidgets.QDoubleSpinBox(decimals=3, minimum=0.001, maximum=500.0, value=10.0)
        self.unit = QtWidgets.QComboBox()
        self.unit.addItems(['ppm', 'Da'])
        self.unit.currentTextChanged.connect(self._unit_changed)
        self.min_rel = QtWidgets.QDoubleSpinBox(decimals=1, minimum=0.0, maximum=100.0, value=2.0, suffix=' %')
        self.min_rel.setToolTip('Peaks below this share of the tallest are left out')
        self.cleave = QtWidgets.QCheckBox('Bond cleavages')
        self.cleave.setChecked(True)
        self.cleave.setToolTip('Also the pieces of breaking one bond (outside rings) of each structure: [PH+H]+, '
                               '[PH]•+ and P+ (negative: [PH−H]−)')
        for w in (open_b, paste_b, self.which):
            top.addWidget(w)
        top.addStretch(1)
        lay.addLayout(top)
        opts = QtWidgets.QHBoxLayout()
        for w in (self.mode, QtWidgets.QLabel('within'), self.tol, self.unit, QtWidgets.QLabel('peaks over'),
                  self.min_rel, self.cleave):
            opts.addWidget(w)
        opts.addStretch(1)
        lay.addLayout(opts)
        for w in (self.mode, self.unit):
            w.currentIndexChanged.connect(self.update_view)
        for w in (self.tol, self.min_rel):
            w.valueChanged.connect(self.update_view)
        self.cleave.toggled.connect(self.update_view)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.plot = PlotCanvas(self)
        self.plot.x_only = True
        split.addWidget(self.plot)
        self.table = QtWidgets.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(['m/z', '%', 'Ion', 'Formula', 'm/z calc.', 'Error', 'Structure'])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.cellDoubleClicked.connect(self._select)
        split.addWidget(self.table)
        split.setSizes([420, 260])
        lay.addWidget(split, 1)
        bottom = QtWidgets.QHBoxLayout()
        self.note = QtWidgets.QLabel('Open a spectrum, or paste m/z and intensity pairs.')
        self.note.setWordWrap(True)
        bottom.addWidget(self.note, 1)
        copy_b = QtWidgets.QPushButton('Copy table')
        copy_b.clicked.connect(self.copy_table)
        insert_b = QtWidgets.QPushButton('Insert into drawing')
        insert_b.setToolTip('The annotated spectrum as a picture beside the structures')
        insert_b.clicked.connect(self.insert)
        close_b = QtWidgets.QPushButton('Close')
        close_b.clicked.connect(self.close)
        for b in (copy_b, insert_b, close_b):
            bottom.addWidget(b)
        lay.addLayout(bottom)

    # ------------------------------------------------------------------------------------------------- loading
    def open_file(self):
        path, _f = QtWidgets.QFileDialog.getOpenFileName(self, 'Open mass spectrum', '', FILE_FILTER)
        if path:
            self.load_path(path)

    def load_path(self, path):
        try:
            spectra = msfile.read(path)
        except Exception as e:                              # a damaged or unknown file: say so
            QtWidgets.QMessageBox.warning(self, 'Annotate mass spectrum', f'Could not read the file:\n{e}')
            return False
        return self.set_spectra(spectra, path)

    def paste(self):
        text = QtWidgets.QApplication.clipboard().text()
        self.set_spectra(msfile.parse_text(text), 'the clipboard')

    def set_spectra(self, spectra, source=''):
        spectra = [s for s in spectra if len(s.mz)]
        if not spectra:
            QtWidgets.QMessageBox.information(self, 'Annotate mass spectrum', f'No spectrum found in {source}.')
            return False
        self.spectra = spectra
        self.which.blockSignals(True)
        self.which.clear()
        for s in spectra:
            bits = [s.title]
            if s.rt is not None:
                bits.append(f'{s.rt:.2f} min')
            if s.level > 1:
                bits.append(f'MS{s.level}' + (f' of {s.precursor:.4f}' if s.precursor else ''))
            self.which.addItem(', '.join(bits))
        best = max(range(len(spectra)), key=lambda k: spectra[k].tic) if len(spectra) > 1 else 0
        self.which.setCurrentIndex(best)
        self.which.blockSignals(False)
        self.which.setVisible(len(spectra) > 1)
        peaks = msfile.centroid(spectra[best])
        unit_res = msfile.unit_resolution(peaks)
        self.unit.blockSignals(True)
        self.unit.setCurrentText('Da' if unit_res else 'ppm')
        self.tol.setValue(0.4 if unit_res else 10.0)
        self.unit.blockSignals(False)
        self.update_view()
        return True

    def _unit_changed(self, unit):
        self.tol.blockSignals(True)
        self.tol.setValue(0.4 if unit == 'Da' else 10.0)
        self.tol.blockSignals(False)

    # ------------------------------------------------------------------------------------------------- matching
    def units(self):
        doc = self.canvas.doc
        return [(f'{k}: {chem.analyse(doc, set(u)).formula}', set(u))
                for k, u in enumerate(self.canvas.analysis_units(set(doc.atoms)), 1)]

    def update_view(self, *_):
        if not self.spectra:
            return
        sp = self.spectra[max(0, self.which.currentIndex())]
        peaks = msfile.centroid(sp, self.min_rel.value() / 100.0)
        if len(peaks) > 400:                                # the tallest; a full profile scan would drown the plot
            peaks = sorted(sorted(peaks, key=lambda p: -p[1])[:400])
        cands = ms.candidate_ions(self.canvas.doc, self.units(), self.mode.currentData(), self.cleave.isChecked())
        self.rows = ms.annotate(peaks, cands, self.tol.value(), self.unit.currentText())
        self._fill()
        self._plot()
        n = sum(1 for r in self.rows if r[2])
        self.note.setText(f'{n} of {len(self.rows)} peaks matched to {len(cands)} possible ions of the structures '
                          'drawn.' + ('' if self.canvas.doc.atoms else ' Draw the structure to match its ions.'))

    def _fill(self):
        top = max((r[1] for r in self.rows), default=1.0) or 1.0
        self.table.setRowCount(len(self.rows))
        unit = self.unit.currentText()
        for k, (mz, y, m, err) in enumerate(sorted(self.rows, key=lambda r: -r[1])):
            vals = [f'{mz:.4f}', f'{100 * y / top:.1f}']
            vals += [m[1], m[2], f'{m[3]:.4f}', f'{err:+.1f} ppm' if unit == 'ppm' else f'{err:+.3f} Da', m[0]] if m \
                else [''] * 5
            for c, v in enumerate(vals):
                it = QtWidgets.QTableWidgetItem(v)
                if c in (0, 1, 4, 5):
                    it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
                if m:
                    it.setData(QtCore.Qt.ItemDataRole.UserRole, sorted(m[4]))
                self.table.setItem(k, c, it)
        self.table.resizeColumnsToContents()

    def _plot(self):
        T = theme.T
        ax = self.plot.layout(1)[0]
        draw_spectrum(ax, self.rows, T['fig_fg'], T['fig_fg'], MATCH_COLOR)
        self.plot.style()
        self.plot.apply_view()
        self.plot.draw_idle()

    def _select(self, row, _col):
        it = self.table.item(row, 0)
        ids = it.data(QtCore.Qt.ItemDataRole.UserRole) if it else None
        if ids:
            have = [('atom', a) for a in ids if a in self.canvas.doc.atoms]
            if have:
                self.canvas.set_selection(have)

    # ------------------------------------------------------------------------------------------------- output
    def table_text(self):
        heads = [self.table.horizontalHeaderItem(c).text() for c in range(self.table.columnCount())]
        lines = ['\t'.join(heads)]
        for r in range(self.table.rowCount()):
            lines.append('\t'.join(self.table.item(r, c).text() if self.table.item(r, c) else ''
                                   for c in range(self.table.columnCount())))
        return '\n'.join(lines)

    def copy_table(self):
        QtWidgets.QApplication.clipboard().setText(self.table_text())

    def image(self, dpi=200):
        """The annotated spectrum as a QImage (black on white, the zoom of the plot kept)."""
        from matplotlib.figure import Figure
        fig = Figure(figsize=(6.4, 3.0), layout='constrained')
        ax = fig.add_subplot(111)
        draw_spectrum(ax, self.rows, 'black', '#333333', MATCH_COLOR)
        if self.plot.view and self.plot.view.get('x'):
            ax.set_xlim(self.plot.view['x'])
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=dpi, facecolor='white')
        return QtGui.QImage.fromData(buf.getvalue(), 'PNG')

    def insert(self):
        if not self.rows:
            return
        self.canvas.insert_image(self.image())
