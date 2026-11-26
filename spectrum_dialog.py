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
ISOTOPE_COLOR = '#F48FB1'                                  # isotope peaks of a matched ion


def ion_text(match):
    """The ion's formula with its charge and radical dot: C2H3O+, C4H8O•+, C4H9O+ for [M+H]+."""
    if not match[1].startswith('['):                        # already the ion's formula (a rule-based fragment)
        return match[1]
    m = re.search(r'(•?)(\d*[+−-])$', match[1])
    return match[2] + (m.group(1) + m.group(2) if m else '')


def draw_spectrum(ax, rows, fg, stick, label_color, most_labels=None):
    """Sticks of the peaks (relative intensity, %), every matched one coloured and labelled with its m/z and ion
    formula. Each label is measured on the figure and lifted until it covers no other label and no peak (a dotted
    line joins a lifted label to its peak); the plot grows upward to hold them all. most_labels: the tallest so many
    only (None: all)."""
    if not rows:
        return
    top = max(r[1] for r in rows) or 1.0
    xs = [r[0] for r in rows]
    span = (max(xs) - min(xs)) or 1.0
    pad = span * 0.05
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    for mz, y, match, _err in rows:
        iso = bool(match) and match[0].startswith('isotope of')
        ax.vlines(mz, 0, 100 * y / top, color=ISOTOPE_COLOR if iso else label_color if match else stick,
                  lw=1.1 if match else 0.8)
    ax.set_xlabel('m/z', color=fg)
    ax.set_ylabel('Relative intensity (%)', color=fg)
    matched = sorted((r for r in rows if r[2] and not r[2][0].startswith('isotope of')), key=lambda r: -r[1])
    if most_labels:
        matched = matched[:most_labels]
    fig = ax.figure
    if not hasattr(fig.canvas, 'get_renderer'):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        FigureCanvasAgg(fig)
    sticks = [(mz, 100 * y / top) for mz, y, _m, _e in rows]
    ymax = 118.0
    for _attempt in range(6):                               # grow the plot until every label fits
        ax.set_ylim(0, ymax)
        renderer = fig.canvas.get_renderer()
        to_disp = ax.transData.transform
        placed, texts, highest = [], [], 0.0
        for mz, y, match, _err in matched:
            h = 100 * y / top
            label = f'{mz:.4f}' if mz != round(mz) else f'{mz:g}'
            t = ax.text(mz, h + 1.5, label + '\n' + ion_text(match), ha='center', va='bottom', fontsize=7.5,
                        color=label_color)
            base = h + 1.5
            for _k in range(200):
                t.set_y(base)
                bb = t.get_window_extent(renderer).expanded(1.06, 1.05)
                hit = any(bb.overlaps(o) for o in placed)
                if not hit:
                    for sx, sh in sticks:                   # a taller peak running through the label
                        px, py = to_disp((sx, sh))
                        if abs(sx - mz) > 1e-9 and bb.x0 <= px <= bb.x1 and py > bb.y0:
                            hit = True
                            break
                if not hit:
                    break
                base += 1.5
            placed.append(bb)
            texts.append((t, mz, h, base))
            highest = max(highest, ax.transData.inverted().transform((0, bb.y1))[1])
        if highest <= ymax - 2 or _attempt == 5:
            break
        for t, *_r in texts:
            t.remove()
        ymax = highest * 1.08 + 4
    for t, mz, h, base in texts:
        if base > h + 2.5:
            ax.plot([mz, mz], [h + 0.5, base], color=label_color, lw=0.5, ls=':')


class SpectrumDialog(QtWidgets.QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.setWindowTitle('Annotate mass spectrum')
        self.window_ = window
        self.canvas = window.canvas
        self.spectra = []
        self.rows = []
        scr = window.screen().availableGeometry() if window is not None and window.screen() else None
        if scr is not None:                                 # big enough to show every label and column
            self.resize(int(scr.width() * 0.88), int(scr.height() * 0.88))
        else:
            self.resize(1400, 900)
        self.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
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
        self.find_b = QtWidgets.QPushButton('Find the structure in this run…')
        self.find_b.setToolTip('The MS/MS scans of this LC-MS run whose precursor is an ion of a drawn structure, '
                               'scored by how much of their intensity its predicted fragments explain, and the MS1 '
                               'scan where each of its ions is tallest')
        self.find_b.clicked.connect(self.find_in_run)
        self.find_b.setEnabled(False)
        for w in (self.mode, QtWidgets.QLabel('within'), self.tol, self.unit, QtWidgets.QLabel('peaks over'),
                  self.min_rel, self.cleave):
            opts.addWidget(w)
        opts.addStretch(1)
        opts.addWidget(self.find_b)
        lay.addLayout(opts)
        self.rank = QtWidgets.QLabel()
        self.rank.setWordWrap(True)
        self.rank.setToolTip('Explains: the share of the measured intensity matched by the structure\'s ions, '
                             'fragments and their isotope peaks. Similarity: cosine of the measured and predicted '
                             'sticks (a predicted stick as tall as the fragment is likely). Draw several candidate '
                             'structures to compare them.')
        lay.addWidget(self.rank)
        for w in (self.mode, self.unit):
            w.currentIndexChanged.connect(self.update_view)
        for w in (self.tol, self.min_rel):
            w.valueChanged.connect(self.update_view)
        self.cleave.toggled.connect(self.update_view)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.plot = PlotCanvas(self)
        self.plot.x_only = True
        split.addWidget(self.plot)
        self.table = QtWidgets.QTableWidget(0, 9)
        self.table.setHorizontalHeaderLabels(['m/z', '%', 'Ion', 'Formula', 'm/z calc.', 'Error', 'Structure / rule',
                                              'Isotope peaks (measured vs calc.)', 'Other possible ions'])
        self.table.setWordWrap(True)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.cellDoubleClicked.connect(self._select)
        split.addWidget(self.table)
        split.setSizes([560, 340])
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
        self.rows, self.fits = ms.isotopes(self.rows, self.tol.value(), self.unit.currentText())
        self.others = ms.all_matches(peaks, cands, self.tol.value(), self.unit.currentText())
        self._fill()
        self._plot()
        n = sum(1 for r in self.rows if r[2])
        iso = sum(1 for r in self.rows if r[2] and r[2][0].startswith('isotope of'))
        self.note.setText(f'{n} of {len(self.rows)} peaks matched ({iso} of them isotope peaks) to {len(cands)} possible '
                          'ions of the structures drawn.' +
                          ('' if self.canvas.doc.atoms else ' Draw the structure to match its ions.'))
        self.find_b.setEnabled(len(self.spectra) > 1 and bool(self.canvas.doc.atoms))
        self._rank(peaks)

    def _rank(self, peaks):
        """The drawn structures ranked by how well they explain this spectrum."""
        units = self.units()
        if not units or not peaks:
            self.rank.setText('')
            return
        ranked = ms.rank_structures(self.canvas.doc, units, peaks, self.mode.currentData(), self.tol.value(),
                                    self.unit.currentText(), self.cleave.isChecked())
        self.ranked = ranked
        lines = []
        for k, r in enumerate(ranked):
            line = (f'{r["name"]}: explains {100 * r["explained"]:.0f} % of the intensity ({r["matched"]} of '
                    f'{r["peaks"]} peaks), similarity {r["similarity"]:.2f}')
            lines.append(f'<b>{line}</b>' if k == 0 and len(ranked) > 1 else line)
        head = 'Best fit first: ' if len(ranked) > 1 else 'Fit: '
        self.rank.setText(head + '<br>'.join(lines))

    def _fill(self):
        top = max((r[1] for r in self.rows), default=1.0) or 1.0
        self.table.setRowCount(len(self.rows))
        unit = self.unit.currentText()
        others = getattr(self, 'others', {})
        for k, (mz, y, m, err) in enumerate(sorted(self.rows, key=lambda r: -r[1])):
            vals = [f'{mz:.4f}', f'{100 * y / top:.1f}']
            vals += [m[1], m[2], f'{m[3]:.4f}', f'{err:+.1f} ppm' if unit == 'ppm' else f'{err:+.3f} Da', m[0]] if m \
                else [''] * 5
            alt = [c for c, _e in others.get(mz, []) if not m or (c[1], c[0]) != (m[1], m[0])]
            seen, names = set(), []
            for c in alt:
                if c[1] not in seen:
                    seen.add(c[1]); names.append(c[1])
            vals.append(getattr(self, 'fits', {}).get(mz, ''))
            vals.append(', '.join(names))
            for c, v in enumerate(vals):
                it = QtWidgets.QTableWidgetItem(v)
                if c in (0, 1, 4, 5):
                    it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter)
                if m:
                    it.setData(QtCore.Qt.ItemDataRole.UserRole, sorted(m[4]))
                if c == 8 and alt:
                    it.setToolTip('\n'.join(f'{a[1]}  {a[3]:.4f}  ({a[0]})' for a in alt))
                elif c == 6 and m:
                    it.setToolTip(m[0])
                if m and m[0].startswith('isotope of'):
                    it.setForeground(QtGui.QColor(ISOTOPE_COLOR))
                self.table.setItem(k, c, it)
        self.table.resizeColumnsToContents()
        hh = self.table.horizontalHeader()
        for c in (6, 7, 8):
            hh.setSectionResizeMode(c, QtWidgets.QHeaderView.ResizeMode.Interactive)
            self.table.setColumnWidth(c, min(self.table.columnWidth(c), 520))
        self.table.resizeRowsToContents()

    def find_in_run(self):
        """Search the opened LC-MS run for a drawn structure (the best-fitting one, or the one picked)."""
        units = self.units()
        if not units or len(self.spectra) < 2:
            return
        names = [u[0] for u in units]
        k = 0
        if len(units) > 1:
            name, ok = QtWidgets.QInputDialog.getItem(self, 'Find the structure in this run', 'Structure:', names, 0,
                                                      False)
            if not ok:
                return
            k = names.index(name)
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            hits = ms.search_run(self.canvas.doc, units[k][1], self.spectra, self.mode.currentData(),
                                 max(self.tol.value(), 5.0) if self.unit.currentText() == 'ppm' else self.tol.value(),
                                 self.unit.currentText())
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle(f'{names[k]} in this run')
        dlg.resize(820, 520)
        lay = QtWidgets.QVBoxLayout(dlg)
        n2 = sum(1 for h in hits if h['kind'] == 'MS2')
        lay.addWidget(QtWidgets.QLabel(f'{n2} MS/MS scans of its ions; MS1: where each ion is tallest. Double-click a '
                                       'row to open that scan.' if hits else
                                       'None of its ions were found in this run.'))
        t = QtWidgets.QTableWidget(len(hits), 7)
        t.setHorizontalHeaderLabels(['Scan', 'RT (min)', 'Ion', 'Precursor', 'Explained', 'Fragments found',
                                     'Intensity'])
        t.verticalHeader().hide()
        t.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        t.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        for r, h in enumerate(hits):
            sp = self.spectra[h['index']]
            vals = [f'{sp.title} ({h["kind"]})', f'{h["rt"]:.2f}' if h['rt'] is not None else '', h['ion'],
                    f'{h["precursor"]:.4f}' if h['precursor'] else '',
                    f'{100 * h["explained"]:.0f} %' if h['explained'] is not None else '',
                    f'{h["matched"]} of {h["peaks"]}' if h['matched'] is not None else '', f'{h["intensity"]:.3g}']
            for c, v in enumerate(vals):
                it = QtWidgets.QTableWidgetItem(v)
                it.setData(QtCore.Qt.ItemDataRole.UserRole, h['index'])
                t.setItem(r, c, it)
        t.resizeColumnsToContents()

        def open_scan(row, _col):
            self.which.setCurrentIndex(t.item(row, 0).data(QtCore.Qt.ItemDataRole.UserRole))
        t.cellDoubleClicked.connect(open_scan)
        lay.addWidget(t, 1)
        close = QtWidgets.QPushButton('Close')
        close.clicked.connect(dlg.accept)
        lay.addWidget(close, 0, QtCore.Qt.AlignmentFlag.AlignRight)
        self.run_dialog = dlg
        dlg.show()

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
        labels = sum(1 for r in self.rows if r[2])
        fig = Figure(figsize=(max(6.4, min(16.0, 0.42 * labels)), 3.0 + min(2.5, 0.08 * labels)),
                     layout='constrained')
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
