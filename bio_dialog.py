"""
bio_dialog.py - the editors for biopolymer sequences (type or paste FASTA, three-letter text or HELM; termini,
cyclic, bridges, residues per line, style) and plasmid maps (size in bp, name, regions and markers).
"""
import re

from PyQt6 import QtCore, QtWidgets

import biopolymer
import theme


class SequenceDialog(QtWidgets.QDialog):
    """Insert or edit a sequence. values() gives the model.Sequence fields to set."""

    def __init__(self, seq=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Biopolymer sequence' if seq is None else 'Edit sequence')
        theme.register_window(self)
        lay = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        self.polymer = QtWidgets.QComboBox()
        for p, name in (('PEPTIDE', 'Peptide / protein'), ('DNA', 'DNA'), ('RNA', 'RNA')):
            self.polymer.addItem(name, p)
        self.polymer.addItem('Detect from the text', None)
        form.addRow('Type:', self.polymer)
        self.text = QtWidgets.QPlainTextEdit()
        self.text.setPlaceholderText('One-letter (FASTA) codes: GIVEQCCTSICSLYQLENYCN\n'
                                     'Three-letter text: H-Ala-Gly-Pro-NH2\n'
                                     'HELM: PEPTIDE1{A.C.D}$$$$  or  RNA1{R(A)P.R(C)P.R(G)}$$$$')
        self.text.setMinimumHeight(90)
        form.addRow('Sequence:', self.text)
        self.n_term = QtWidgets.QComboBox(); self.n_term.addItem('H (free amine)', 'H'); self.n_term.addItem('Ac (acetyl)', 'Ac')
        self.c_term = QtWidgets.QComboBox(); self.c_term.addItem('OH (free acid)', 'OH'); self.c_term.addItem('NH₂ (amide)', 'NH2')
        row = QtWidgets.QHBoxLayout(); row.addWidget(QtWidgets.QLabel('N-terminus:')); row.addWidget(self.n_term)
        row.addSpacing(10); row.addWidget(QtWidgets.QLabel('C-terminus:')); row.addWidget(self.c_term); row.addStretch(1)
        self.term_row = QtWidgets.QWidget(); self.term_row.setLayout(row)
        form.addRow('', self.term_row)
        self.cyclic = QtWidgets.QCheckBox('Cyclic (head-to-tail amide)')
        form.addRow('', self.cyclic)
        self.bridges = QtWidgets.QLineEdit()
        self.bridges.setPlaceholderText('e.g. 6-11, 7-20  (Cys-Cys disulfide, Lys-Asp / Lys-Glu lactam)')
        form.addRow('Bridges:', self.bridges)
        self.p5 = QtWidgets.QCheckBox("5′-phosphate"); self.p3 = QtWidgets.QCheckBox("3′-phosphate")
        prow = QtWidgets.QHBoxLayout(); prow.addWidget(self.p5); prow.addWidget(self.p3); prow.addStretch(1)
        self.phos_row = QtWidgets.QWidget(); self.phos_row.setLayout(prow)
        form.addRow('', self.phos_row)
        self.per_line = QtWidgets.QSpinBox(); self.per_line.setRange(1, 100); self.per_line.setValue(10)
        form.addRow('Residues per line:', self.per_line)
        self.three = QtWidgets.QCheckBox('Three-letter labels (Ala, Gly ...)')
        self.ids = QtWidgets.QCheckBox('Residue numbers')
        srow = QtWidgets.QHBoxLayout(); srow.addWidget(self.three); srow.addWidget(self.ids); srow.addStretch(1)
        form.addRow('Display:', srow)
        lay.addLayout(form)
        self.info = QtWidgets.QLabel(); self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.info)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                        QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        if seq is not None:
            self.polymer.setCurrentIndex(biopolymer.POLYMERS.index(seq.polymer))
            self.text.setPlainText(''.join(seq.residues))
            self.n_term.setCurrentIndex(0 if seq.n_term == 'H' else 1)
            self.c_term.setCurrentIndex(0 if seq.c_term == 'OH' else 1)
            self.cyclic.setChecked(seq.cyclic)
            self.bridges.setText(', '.join(f'{i + 1}-{j + 1}' for i, j in seq.bridges))
            self.p5.setChecked(seq.p5); self.p3.setChecked(seq.p3)
            self.per_line.setValue(seq.per_line); self.three.setChecked(seq.three_letter); self.ids.setChecked(seq.show_ids)
        else:
            self.polymer.setCurrentIndex(3)
        for w in (self.polymer,):
            w.currentIndexChanged.connect(self._update)
        for w in (self.n_term, self.c_term):
            w.currentIndexChanged.connect(self._update)
        for w in (self.cyclic, self.p5, self.p3):
            w.toggled.connect(self._update)
        self.text.textChanged.connect(self._update)
        self.bridges.textChanged.connect(self._update)
        self.resize(560, 480)
        self._update()

    def _values(self):
        """The fields, or raises ValueError."""
        polymer, residues, extras = biopolymer.parse(self.text.toPlainText(), self.polymer.currentData())
        v = dict(polymer=polymer, residues=residues, per_line=self.per_line.value(),
                 three_letter=self.three.isChecked(), show_ids=self.ids.isChecked(),
                 cyclic=self.cyclic.isChecked() and polymer == 'PEPTIDE',
                 n_term=self.n_term.currentData(), c_term=self.c_term.currentData(),
                 p5=self.p5.isChecked(), p3=self.p3.isChecked(), bridges=[])
        for k in ('n_term', 'c_term', 'p5', 'p3', 'cyclic'):  # what the text itself says (HELM, H-...-NH2) wins
            if k in extras:
                v[k] = extras[k]
        text = self.bridges.text().strip()
        if text:
            for part in re.split(r'[,;\s]+', text):
                m = re.match(r'^(\d+)\s*[-–:]\s*(\d+)$', part)
                if not m:
                    raise ValueError(f'Bridge "{part}": write two residue numbers, e.g. 6-11.')
                i, j = int(m.group(1)) - 1, int(m.group(2)) - 1
                v['bridges'].append([min(i, j), max(i, j)])
        elif 'bridges' in extras:
            v['bridges'] = extras['bridges']
        if polymer != 'PEPTIDE':
            v['bridges'] = []
        return v

    def _update(self):
        pep = self.polymer.currentData() in ('PEPTIDE', None)
        try:
            v = self._values()
            pep = v['polymer'] == 'PEPTIDE'
            from model import Sequence
            seq = Sequence(0, 0, **v)
            info = biopolymer.analyse(seq)
            if info.errors:
                self.info.setText(f'{biopolymer.summary(seq)}. <span style="color:#B3261E">{info.errors[0]}</span>')
            else:
                self.info.setText(f'{biopolymer.summary(seq)} · {info.formula_html} · exact mass {info.exact:.4f} · '
                                  f'MW {info.average:.2f}')
        except ValueError as e:
            self.info.setText(f'<span style="color:#B3261E">{e}</span>' if self.text.toPlainText().strip() else '')
        self.term_row.setEnabled(pep and not self.cyclic.isChecked())
        self.cyclic.setEnabled(pep); self.bridges.setEnabled(pep); self.three.setEnabled(pep)
        self.phos_row.setEnabled(not pep)

    def _accept(self):
        try:
            self.result_values = self._values()
        except ValueError as e:
            QtWidgets.QMessageBox.warning(self, self.windowTitle(), str(e))
            return
        self.accept()

    def values(self):
        return self.result_values


class PlasmidDialog(QtWidgets.QDialog):
    """Size, name, regions (start, end, label) and markers (position, label) of a plasmid map."""

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Plasmid map')
        theme.register_window(self)
        lay = QtWidgets.QVBoxLayout(self)
        form = QtWidgets.QFormLayout()
        self.name = QtWidgets.QLineEdit(params.get('name', ''))
        self.name.setPlaceholderText('e.g. pUC19')
        form.addRow('Name:', self.name)
        self.bp = QtWidgets.QSpinBox(); self.bp.setRange(10, 10_000_000); self.bp.setValue(int(params.get('bp', 5000)))
        self.bp.setSuffix(' bp')
        form.addRow('Size:', self.bp)
        lay.addLayout(form)
        lay.addWidget(QtWidgets.QLabel('Regions (drawn as arrows from start to end, clockwise):'))
        self.regions = self._table(['Start', 'End', 'Label'], [(r['start'], r['end'], r.get('label', ''))
                                                            for r in params.get('regions', [])])
        lay.addWidget(self.regions)
        lay.addWidget(QtWidgets.QLabel('Markers (a tick and a label, e.g. restriction sites):'))
        self.markers = self._table(['Position', 'Label'], [(m['pos'], m.get('label', ''))
                                                           for m in params.get('markers', [])])
        lay.addWidget(self.markers)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                        QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self.resize(460, 520)

    def _table(self, headers, rows):
        box = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(box); v.setContentsMargins(0, 0, 0, 0)
        t = QtWidgets.QTableWidget(0, len(headers))
        t.setHorizontalHeaderLabels(headers)
        t.horizontalHeader().setStretchLastSection(True)
        t.verticalHeader().hide()
        for row in rows:
            r = t.rowCount(); t.insertRow(r)
            for c, val in enumerate(row):
                t.setItem(r, c, QtWidgets.QTableWidgetItem(str(val)))
        h = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton('Add'); rem = QtWidgets.QPushButton('Remove')
        add.clicked.connect(lambda: (t.insertRow(t.rowCount()), t.setCurrentCell(t.rowCount() - 1, 0)))
        rem.clicked.connect(lambda: t.removeRow(t.currentRow()) if t.currentRow() >= 0 else None)
        h.addWidget(add); h.addWidget(rem); h.addStretch(1)
        v.addWidget(t); v.addLayout(h)
        box.table = t
        return box

    @staticmethod
    def _rows(box):
        t = box.table
        out = []
        for r in range(t.rowCount()):
            cells = [(t.item(r, c).text().strip() if t.item(r, c) else '') for c in range(t.columnCount())]
            if any(cells):
                out.append((r + 1, cells))
        return out

    def _accept(self):
        bp = self.bp.value()
        regions, markers = [], []
        try:
            for n, (start, end, label) in self._rows(self.regions):
                s, e = int(start), int(end)
                if not (1 <= s <= bp and 1 <= e <= bp):
                    raise ValueError(f'Region {n}: start and end must lie between 1 and {bp}.')
                regions.append({'start': s, 'end': e, 'label': label})
            for n, (pos, label) in self._rows(self.markers):
                p = int(pos)
                if not 1 <= p <= bp:
                    raise ValueError(f'Marker {n}: the position must lie between 1 and {bp}.')
                markers.append({'pos': p, 'label': label})
        except ValueError as e:
            msg = str(e) if 'must' in str(e) else 'Positions are whole numbers of base pairs.'
            QtWidgets.QMessageBox.warning(self, 'Plasmid map', msg)
            return
        self.params = {'bp': bp, 'name': self.name.text().strip(), 'regions': regions, 'markers': markers}
        self.accept()


class TLCDialog(QtWidgets.QDialog):
    """The Rf values of a TLC plate, one line per lane (e.g. 0.25, 0.6); a blank line leaves an empty lane out."""

    def __init__(self, plate, parent=None):
        super().__init__(parent)
        self.setWindowTitle('TLC plate')
        theme.register_window(self)
        self.plate = plate
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(QtWidgets.QLabel('Rf values of the spots, one lane per line (0 = origin, 1 = solvent front):'))
        self.text = QtWidgets.QPlainTextEdit()
        self.text.setPlainText('\n'.join(', '.join(f"{s['rf']:.2f}" for s in lane) for lane in plate.lanes))
        lay.addWidget(self.text)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok |
                                        QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(self._accept); bb.rejected.connect(self.reject)
        lay.addWidget(bb)
        self.resize(380, 300)

    def _accept(self):
        import tlc
        lanes = []
        for n, line in enumerate(self.text.toPlainText().splitlines(), 1):
            vals = [v for v in re.split(r'[,;\s]+', line.strip()) if v]
            if not vals:
                continue
            old = self.plate.lanes[n - 1] if n - 1 < len(self.plate.lanes) else []
            lane = []
            for k, v in enumerate(vals):
                try:
                    rf = float(v.replace(',', '.'))
                except ValueError:
                    QtWidgets.QMessageBox.warning(self, 'TLC plate', f'Lane {n}: "{v}" is not a number.')
                    return
                if not 0 <= rf <= 1:
                    QtWidgets.QMessageBox.warning(self, 'TLC plate', f'Lane {n}: Rf {v} is outside 0-1.')
                    return
                lane.append(dict(old[k], rf=round(rf, 3)) if k < len(old) else tlc.new_spot(rf))
            lanes.append(lane)
        if not lanes:
            QtWidgets.QMessageBox.warning(self, 'TLC plate', 'A plate needs at least one lane with a spot.')
            return
        self.lanes = lanes
        self.accept()
