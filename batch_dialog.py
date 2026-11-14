"""File > Batch export: the dialog."""
import os

from PyQt6 import QtCore, QtWidgets

import batch


class BatchDialog(QtWidgets.QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window_ = window
        self.setWindowTitle('Batch export')
        self.resize(560, 260)
        form = QtWidgets.QFormLayout(self)
        row = QtWidgets.QHBoxLayout()
        self.src = QtWidgets.QLineEdit()
        self.src.setPlaceholderText('an .sdf, or a .smi / .txt / .csv list (SMILES, then a name)')
        pick = QtWidgets.QPushButton('Choose…'); pick.clicked.connect(self._pick)
        row.addWidget(self.src, 1); row.addWidget(pick)
        form.addRow('Structures', row)
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems(['One grid sheet', 'One picture per structure'])
        form.addRow('Make', self.mode)
        self.fmt = QtWidgets.QComboBox()
        self.fmt.addItems(['PNG', 'SVG', 'PDF'])
        form.addRow('Format', self.fmt)
        self.cols = QtWidgets.QSpinBox(minimum=1, maximum=12, value=4)
        form.addRow('Columns (grid)', self.cols)
        self.names = QtWidgets.QCheckBox('Name under each structure'); self.names.setChecked(True)
        form.addRow('', self.names)
        self.style_ = QtWidgets.QCheckBox("This drawing's style (bond length, font, lines)"); self.style_.setChecked(True)
        form.addRow('', self.style_)
        self.note = QtWidgets.QLabel('')
        form.addRow(self.note)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        go = buttons.addButton('Export…', QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        go.clicked.connect(self.run)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def _pick(self):
        path, _f = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Structures to export', self.window_._dir(),
            'Structure lists (*.sdf *.sd *.smi *.smiles *.txt *.csv *.mol *.cdx)')
        if path:
            self.src.setText(path)

    def run(self):
        path = self.src.text().strip()
        try:
            items = batch.read(path)
        except (OSError, ValueError) as e:
            self.note.setText(str(e) or 'The file could not be read.')
            return
        ext = '.' + self.fmt.currentText().lower()
        style = dict(self.window_.canvas.doc.style) if self.style_.isChecked() else None
        dpi = int(self.window_.settings.value('dpi', 600))
        transparent = self.window_.a_transparent.isChecked()
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            if self.mode.currentIndex() == 0:
                base = os.path.splitext(path)[0] + '_grid' + ext
                out, _f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save the grid sheet', base,
                                                                f'{self.fmt.currentText()} (*{ext})')
                if not out:
                    return
                if not out.lower().endswith(ext):
                    out += ext
                written = batch.export_grid(items, out, self.cols.value(), style, self.names.isChecked(), dpi,
                                            transparent)
            else:
                folder = QtWidgets.QFileDialog.getExistingDirectory(self, 'Folder for the pictures',
                                                                    os.path.dirname(path))
                if not folder:
                    return
                written = batch.export_each(items, folder, ext, style, self.names.isChecked(), dpi, transparent)
        except OSError as e:
            self.note.setText(f'Could not export: {e}')
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.window_.statusBar().showMessage(f'Batch export: {len(items)} structures, {len(written)} file(s) written.',
                                             8000)
        self.accept()
