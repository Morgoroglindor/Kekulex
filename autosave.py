"""Autosave and crash recovery.

While a drawing has unsaved changes, a copy goes to the autosave folder every INTERVAL_S seconds (next to the log,
%APPDATA%/Kekulex/autosave). Saving, starting a new drawing, opening a file or closing the window normally removes
the copy, so a copy that is still there when Kekulex starts - and whose program is no longer running - is a drawing
lost to a crash or a power cut, and Kekulex offers to bring it back.

Each copy is <pid>-<n>.sed (the drawing, as a .sed file) with <pid>-<n>.json beside it ({path, saved, pid}).
"""
import json
import os
import sys
import time

import applog

DIR = os.environ.get('KEKULEX_AUTOSAVE_DIR') or os.path.join(applog.LOG_DIR, 'autosave')
INTERVAL_S = 60
_count = 0


def _alive(pid):
    """True while process pid is running."""
    if pid == os.getpid():
        return True
    if sys.platform == 'win32':
        import ctypes
        k32 = ctypes.WinDLL('kernel32', use_last_error=True)
        h = k32.OpenProcess(0x1000, False, pid)              # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return ctypes.get_last_error() == 5             # access denied: it exists
        try:
            code = ctypes.c_ulong()
            return bool(k32.GetExitCodeProcess(h, ctypes.byref(code))) and code.value == 259   # STILL_ACTIVE
        finally:
            k32.CloseHandle(h)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class Slot:
    """One window's autosave copy."""

    def __init__(self, folder=None):
        global _count
        _count += 1
        self.folder = folder or DIR
        base = os.path.join(self.folder, f'{os.getpid()}-{_count}')
        self.sed, self.meta = base + '.sed', base + '.json'

    def write(self, doc, path=None):
        """Writes the copy (the drawing doc, which was last saved to path or never); False when it could not."""
        try:
            os.makedirs(self.folder, exist_ok=True)
            tmp = self.sed + '.tmp'
            with open(tmp, 'w', encoding='utf-8') as f:
                json.dump(doc.to_dict(), f)
            os.replace(tmp, self.sed)                       # never a half-written copy
            with open(self.meta, 'w', encoding='utf-8') as f:
                json.dump(dict(path=path or '', saved=time.time(), pid=os.getpid()), f)
            return True
        except OSError:
            applog.log.warning('autosave failed', exc_info=True)
            return False

    def clear(self):
        for p in (self.sed, self.meta, self.sed + '.tmp'):
            try:
                os.remove(p)
            except OSError:
                pass


def orphans(folder=None):
    """[(sed path, meta dict)] of copies left by programs that are no longer running, newest first."""
    folder = folder or DIR
    out = []
    try:
        names = os.listdir(folder)
    except OSError:
        return out
    for name in names:
        if not name.endswith('.sed'):
            continue
        sed = os.path.join(folder, name)
        meta = {}
        try:
            with open(sed[:-4] + '.json', encoding='utf-8') as f:
                meta = json.load(f)
        except (OSError, ValueError):
            pass
        pid = meta.get('pid')
        if pid is None:
            try:
                pid = int(name.split('-')[0])
            except ValueError:
                pid = -1
        if pid > 0 and _alive(pid):
            continue
        meta.setdefault('saved', os.path.getmtime(sed))
        out.append((sed, meta))
    out.sort(key=lambda o: -o[1]['saved'])
    return out


def discard(sed):
    for p in (sed, sed[:-4] + '.json', sed + '.tmp'):
        try:
            os.remove(p)
        except OSError:
            pass
