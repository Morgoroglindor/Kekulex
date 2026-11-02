"""
applog.py - log file and crash handler for Kekulex.

    applog.setup()          # once, at start-up

Messages go to %APPDATA%\\Kekulex\\kekulex.log (rotated at 1 MB, 3 old files kept). An exception nobody caught is written
there with its traceback and shown in a dialog, instead of vanishing (the installed .exe has no console).
"""
import logging
import logging.handlers
import os
import sys
import threading
import traceback

from version import __version__

LOG_DIR = os.path.join(os.environ.get('APPDATA') or os.path.expanduser('~/.config'), 'Kekulex')
LOG_FILE = os.path.join(LOG_DIR, 'kekulex.log')
log = logging.getLogger('kekulex')
_showing = False


def setup(level=logging.INFO):
    """Log to LOG_FILE (and to stderr when there is a console) and install the crash handler. Safe to call twice."""
    if getattr(log, '_lcv_setup', False):
        return
    log.setLevel(level)
    fmt = logging.Formatter('%(asctime)s %(levelname)-7s %(name)s: %(message)s')
    try:
        os.makedirs(LOG_DIR, exist_ok=True)
        fh = logging.handlers.RotatingFileHandler(LOG_FILE, maxBytes=1_000_000, backupCount=3, encoding='utf-8')
        fh.setFormatter(fmt); log.addHandler(fh)
    except OSError:
        pass                                                # read-only profile: keep the console handler only
    if sys.stderr is not None and getattr(sys.stderr, 'isatty', lambda: False)():
        sh = logging.StreamHandler(); sh.setFormatter(fmt); log.addHandler(sh)
    log._lcv_setup = True
    sys.excepthook = _excepthook
    threading.excepthook = lambda a: _excepthook(a.exc_type, a.exc_value, a.exc_traceback)
    log.info('Kekulex %s, Python %s', __version__, sys.version.split()[0])


def _excepthook(etype, value, tb):
    global _showing
    if issubclass(etype, KeyboardInterrupt):
        sys.__excepthook__(etype, value, tb); return
    log.critical('Unhandled exception', exc_info=(etype, value, tb))
    try:
        from PyQt6 import QtWidgets
        if QtWidgets.QApplication.instance() is None or _showing or threading.current_thread() is not threading.main_thread():
            return
        _showing = True                                     # one dialog at a time, even if the error repeats
        try:
            box = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Icon.Critical, 'Kekulex',
                                        f'An unexpected error occurred:\n\n{etype.__name__}: {value}\n\n'
                                        f'The details are in\n{LOG_FILE}')
            box.setDetailedText(''.join(traceback.format_exception(etype, value, tb)))
            box.exec()
        finally:
            _showing = False
    except Exception:
        pass                                                # never let the handler itself crash the program
