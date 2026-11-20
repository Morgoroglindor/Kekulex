"""
notebook_client.py - the lab notebook, from LC viewer, Xpectra and Kekulex: send files and drawings to it, and open
its files ("Open in LC viewer" / "Open in Xpectra" on a file's page: lcviewer://open?server=…&dataset=…).

The canonical copy is LabNotebook_project/tools/desktop/notebook_client.py; each program has an identical copy.

The notebook's address and the personal access key (made on the notebook's "Access keys" page) are asked once and
kept for this Windows user, shared by the three programs (QSettings LabNotebook/Client).

    send_files(parent, paths, program='LC viewer', hint='lc')     the dialog: sample code (found in the file name)
    send_drawing(parent, doc_dict, molecules, smiles, program='Kekulex')
    open_link(url, parent) -> local path of the downloaded raw file, or None
"""
import json
import mimetypes
import os
import re
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import uuid
import webbrowser
import zipfile

from PyQt6 import QtCore, QtWidgets

SCHEMES = ('lcviewer', 'xpectra')
SAMPLE_CODE = re.compile(r'(?<![A-Za-z0-9-])([A-Z][A-Z0-9]{0,11}-(?:\d{3,}-[1-9]\d*|S\d{3,}))(?![A-Za-z0-9-])')
ENTRY_CODE = re.compile(r'(?<![A-Za-z0-9-])([A-Z][A-Z0-9]{0,11}-\d{3,})(?![A-Za-z0-9-])')


def _settings():
    return QtCore.QSettings('LabNotebook', 'Client')


def config():
    s = _settings()
    return (s.value('server', '') or '').rstrip('/'), s.value('key', '') or ''


def save_config(server, key):
    s = _settings()
    s.setValue('server', server.strip().rstrip('/'))
    s.setValue('key', key.strip())


def find_code(text, pattern=SAMPLE_CODE):
    m = pattern.search(os.path.basename(str(text)).upper())
    return m.group(1) if m else ''


class NotebookError(Exception):
    pass


def request(method, path, data=None, files=None, server=None, key=None, raw=False):
    """Call the notebook's API with the access key; returns the JSON (or bytes and headers when raw)."""
    server, key = (server, key) if server else config()
    if not server or not key:
        raise NotebookError('The notebook address and access key are not set.')
    headers = {'Authorization': f'Key {key}', 'Accept': 'application/json'}
    body = None
    if files:
        b = uuid.uuid4().hex
        parts = []
        for k, v in (data or {}).items():
            parts.append(f'--{b}\r\nContent-Disposition: form-data; name="{k}"\r\n\r\n{v}\r\n'.encode())
        for k, (name, content) in files.items():
            ctype = mimetypes.guess_type(name)[0] or 'application/octet-stream'
            parts.append(f'--{b}\r\nContent-Disposition: form-data; name="{k}"; filename="{name}"\r\nContent-Type: {ctype}\r\n\r\n'.encode())
            parts.append(content)
            parts.append(b'\r\n')
        parts.append(f'--{b}--\r\n'.encode())
        body = b''.join(parts)
        headers['Content-Type'] = f'multipart/form-data; boundary={b}'
    elif data is not None:
        body = json.dumps(data).encode()
        headers['Content-Type'] = 'application/json'
    req = urllib.request.Request(f'{server}/api/{path.lstrip("/")}', data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            content = r.read()
            if raw:
                return content, dict(r.headers)
            return json.loads(content or b'null')
    except urllib.error.HTTPError as e:
        try:
            detail = json.loads(e.read() or b'{}')
            msg = detail.get('detail') if isinstance(detail, dict) and 'detail' in detail else '; '.join(
                f'{k}: {" ".join(map(str, v)) if isinstance(v, list) else v}' for k, v in (detail.items() if isinstance(detail, dict) else []))
        except ValueError:
            msg = ''
        if e.code in (401, 403) and not msg:
            msg = 'The access key was refused (revoked, or for another notebook).'
        raise NotebookError(msg or f'The notebook answered {e.code}.') from e
    except urllib.error.URLError as e:
        raise NotebookError(f'Cannot reach {server}: {e.reason}') from e


def _file_bytes(path):
    """A file as it is; a folder (Agilent .D, Waters .raw, Bruker .d) zipped, as the notebook expects it."""
    if not os.path.isdir(path):
        with open(path, 'rb') as f:
            return os.path.basename(path), f.read()
    root = os.path.normpath(path)
    tmp = tempfile.SpooledTemporaryFile(max_size=64 * 1024 * 1024)
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
        for dirpath, _dirs, names in os.walk(root):
            for n in names:
                full = os.path.join(dirpath, n)
                z.write(full, os.path.join(os.path.basename(root), os.path.relpath(full, root)))
    tmp.seek(0)
    return os.path.basename(root) + '.zip', tmp.read()


# ------------------------------------------------------------------------------------------------------- dialogs
class SetupDialog(QtWidgets.QDialog):
    """The notebook's address and the access key (once per Windows user, for all three programs)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Lab notebook: connection')
        server, key = config()
        lay = QtWidgets.QFormLayout(self)
        self.server = QtWidgets.QLineEdit(server or 'https://')
        self.key = QtWidgets.QLineEdit(key)
        self.key.setEchoMode(QtWidgets.QLineEdit.EchoMode.Password)
        lay.addRow('Notebook address', self.server)
        lay.addRow('Access key', self.key)
        hint = QtWidgets.QLabel('Make a key on the notebook: click your name at the top right → “Make a key”. It is kept for this '
                                'Windows user and used by LC viewer, Xpectra and Kekulex.')
        hint.setWordWrap(True)
        lay.addRow(hint)
        self.status = QtWidgets.QLabel('')
        lay.addRow(self.status)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        test = bb.addButton('Test', QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        test.clicked.connect(self.test)
        bb.accepted.connect(self.accept_if_ok)
        bb.rejected.connect(self.reject)
        lay.addRow(bb)
        self.resize(520, self.sizeHint().height())

    def _try(self):
        me = request('GET', 'auth/me/', server=self.server.text().strip().rstrip('/'), key=self.key.text().strip())
        if not me or not me.get('user'):
            raise NotebookError('The key was not accepted.')
        return me['user']

    def test(self):
        try:
            u = self._try()
            self.status.setText(f'✔ Connected as {u["name"]}')
        except NotebookError as e:
            self.status.setText(f'✖ {e}')

    def accept_if_ok(self):
        try:
            self._try()
        except NotebookError as e:
            self.status.setText(f'✖ {e}')
            return
        save_config(self.server.text(), self.key.text())
        self.accept()


def ensure_config(parent):
    server, key = config()
    if server and key:
        return True
    return SetupDialog(parent).exec() == QtWidgets.QDialog.DialogCode.Accepted


class SendDialog(QtWidgets.QDialog):
    def __init__(self, parent, paths, program, hint):
        super().__init__(parent)
        self.program, self.hint, self.results = program, hint, []
        self.setWindowTitle(f'{program}: send to the lab notebook')
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(QtWidgets.QLabel('Files to send (the raw files, unchanged):'))
        self.list = QtWidgets.QListWidget()
        for p in paths:
            it = QtWidgets.QListWidgetItem(os.path.basename(os.path.normpath(p)))
            it.setData(QtCore.Qt.ItemDataRole.UserRole, p)
            it.setFlags(it.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(QtCore.Qt.CheckState.Checked)
            it.setToolTip(p)
            self.list.addItem(it)
        lay.addWidget(self.list)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel('Sample code'))
        self.code = QtWidgets.QLineEdit(next((find_code(p) for p in paths if find_code(p)), ''))
        self.code.setPlaceholderText('BT-001-1 (empty: the code in each file name)')
        row.addWidget(self.code, 1)
        lay.addLayout(row)
        note = QtWidgets.QLabel('A file named like BT-001-1_run2.lcd goes to sample BT-001-1 by itself.')
        note.setStyleSheet('color: gray')
        lay.addWidget(note)
        self.status = QtWidgets.QLabel('')
        self.status.setWordWrap(True)
        self.status.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.status.setOpenExternalLinks(True)
        lay.addWidget(self.status)
        bb = QtWidgets.QDialogButtonBox()
        self.send_btn = bb.addButton('Send', QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        setup = bb.addButton('Connection…', QtWidgets.QDialogButtonBox.ButtonRole.ActionRole)
        close = bb.addButton(QtWidgets.QDialogButtonBox.StandardButton.Close)
        self.send_btn.clicked.connect(self.send)
        setup.clicked.connect(lambda: SetupDialog(self).exec())
        close.clicked.connect(self.reject)
        lay.addWidget(bb)
        self.resize(560, 380)

    def send(self):
        server, _key = config()
        lines = []
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            for i in range(self.list.count()):
                it = self.list.item(i)
                if it.checkState() != QtCore.Qt.CheckState.Checked:
                    continue
                path = it.data(QtCore.Qt.ItemDataRole.UserRole)
                code = self.code.text().strip().upper() or find_code(path)
                try:
                    name, content = _file_bytes(path)
                    r = request('POST', 'upload/', data={'code': code or name, 'hint': self.hint, 'program': self.program},
                                files={'file': (name, content)})
                    url = f'{server}{r["url"]}'
                    state = {'parsed': 'read', 'reading': 'being read', 'failed': 'stored, not read'}.get(r['status'], r['status'])
                    lines.append(f'✔ {name} → <a href="{url}">{r["sample_code"]}</a> ({state})')
                    self.results.append(url)
                except NotebookError as e:
                    lines.append(f'✖ {os.path.basename(path)}: {e}')
                QtWidgets.QApplication.processEvents()
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.status.setText('<br>'.join(lines) or 'Nothing ticked.')


def send_files(parent, paths, program, hint='auto'):
    paths = [p for p in paths if p and os.path.exists(p)]
    if not paths:
        QtWidgets.QMessageBox.information(parent, 'Lab notebook', 'Open a file first.')
        return
    if not ensure_config(parent):
        return
    SendDialog(parent, paths, program, hint).exec()


def send_drawing(parent, doc_dict, molecules, smiles, program='Kekulex'):
    """Put a drawing into a notebook entry (its Structures and reaction). molecules: [{inchikey, smiles, formula, role, mw}]."""
    if not ensure_config(parent):
        return
    code, ok = QtWidgets.QInputDialog.getText(parent, 'Send the drawing to the lab notebook',
                                              'Entry code (e.g. BT-001) — its drawing is replaced; the old one stays in its history:')
    code = (code or '').strip().upper()
    if not ok or not code:
        return
    try:
        target = request('GET', f'resolve/?code={urllib.parse.quote(code)}')
        if target.get('kind') != 'experiment':
            raise NotebookError(f'{code} is a sample: give the code of an entry (like BT-001).')
        exp = request('GET', f'experiments/{target["id"]}/')
        if exp.get('locked'):
            raise NotebookError(f'{code} is signed: its drawing can no longer change.')
        request('PATCH', f'experiments/{exp["id"]}/', data=dict(
            drawing=doc_dict, smiles=smiles, molecules=molecules,
            inchikeys=' '.join(m['inchikey'] for m in molecules if m.get('inchikey')), base_version=exp['version']))
    except NotebookError as e:
        QtWidgets.QMessageBox.warning(parent, 'Lab notebook', str(e))
        return
    server, _ = config()
    if QtWidgets.QMessageBox.question(parent, 'Lab notebook', f'The drawing is in {code}. Open the entry in the browser?') == QtWidgets.QMessageBox.StandardButton.Yes:
        webbrowser.open(f'{server}/e/{exp["id"]}')


def is_link(arg):
    return isinstance(arg, str) and arg.split(':', 1)[0].lower() in SCHEMES and '://' in arg


def open_link(url, parent=None):
    """lcviewer://open?server=…&dataset=… → download the raw file (with the access key) and return its local path."""
    q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
    server, dataset = (q.get('server') or [''])[0].rstrip('/'), (q.get('dataset') or [''])[0]
    saved_server, _key = config()
    if not dataset.isdigit():
        QtWidgets.QMessageBox.warning(parent, 'Lab notebook', f'Not a notebook link: {url}')
        return None
    if server and saved_server != server:
        s = _settings()
        s.setValue('server', server)                    # the link says which notebook; the key may still need setting
    if not ensure_config(parent):
        return None
    try:
        content, headers = request('GET', f'datasets/{dataset}/download/', raw=True)
    except NotebookError as e:
        QtWidgets.QMessageBox.warning(parent, 'Lab notebook', str(e))
        return None
    m = re.search(r'filename="?([^";]+)"?', headers.get('Content-Disposition', ''))
    name = os.path.basename(m.group(1)) if m else f'dataset-{dataset}'
    folder = os.path.join(tempfile.gettempdir(), 'LabNotebook', dataset)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)
    with open(path, 'wb') as f:
        f.write(content)
    if name.lower().endswith('.zip'):                     # a vendor folder: open the folder inside
        with zipfile.ZipFile(path) as z:
            root = os.path.realpath(folder)
            if any(not os.path.realpath(os.path.join(folder, n)).startswith(root + os.sep) for n in z.namelist()):
                QtWidgets.QMessageBox.warning(parent, 'Lab notebook', 'The zip holds names outside its folder: not opened.')
                return None
            z.extractall(folder)
        inner = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
        if inner:
            return os.path.join(folder, inner[0])
    return path
