"""
winkeys.py - keys that other programs take for themselves (Windows only). A program can register a key system-wide
(Notion takes Ctrl+Shift+K, ChemDraw's Clean Up); Windows then hands it to that program even while Kekulex is in front.
While a Kekulex window is active, a low-level keyboard hook sees such keys first, keeps them from the other program and
runs Kekulex's own action instead. The hook is removed whenever Kekulex goes to the background.

    keys = winkeys.KeyClaim(app)
    keys.claim('Ctrl+Shift+K', action.trigger)
"""
import ctypes
import sys
from ctypes import wintypes

from PyQt6 import QtCore, QtGui

WH_KEYBOARD_LL = 13
WM_KEYDOWN, WM_SYSKEYDOWN = 0x0100, 0x0104
VK_SHIFT, VK_CONTROL, VK_MENU, VK_LWIN, VK_RWIN = 0x10, 0x11, 0x12, 0x5B, 0x5C
LLKHF_INJECTED = 0x10


class KBDLLHOOKSTRUCT(ctypes.Structure):
    _fields_ = [('vkCode', wintypes.DWORD), ('scanCode', wintypes.DWORD), ('flags', wintypes.DWORD),
                ('time', wintypes.DWORD), ('dwExtraInfo', ctypes.c_size_t)]


def parse(text):
    """'Ctrl+Shift+K' -> (virtual key, ctrl, shift, alt), for letters and digits."""
    seq = QtGui.QKeySequence(text)[0]
    key, mods = seq.key(), seq.keyboardModifiers()
    M = QtCore.Qt.KeyboardModifier
    vk = int(key.value)
    if not (0x30 <= vk <= 0x39 or 0x41 <= vk <= 0x5A):
        raise ValueError(f'{text}: only letters and digits can be claimed')
    return vk, bool(mods & M.ControlModifier), bool(mods & M.ShiftModifier), bool(mods & M.AltModifier)


def matches(claims, vk, ctrl, shift, alt, win=False):
    """The callback claimed for this key with exactly these modifiers, or None."""
    if win:
        return None
    return claims.get((vk, ctrl, shift, alt))


class KeyClaim(QtCore.QObject):
    def __init__(self, app):
        super().__init__(app)
        self.claims = {}
        self._hook = None
        self._proc = None
        if sys.platform != 'win32':
            return
        self.user32 = ctypes.WinDLL('user32', use_last_error=True)
        self.kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        self.HOOKPROC = ctypes.WINFUNCTYPE(ctypes.c_ssize_t, ctypes.c_int, wintypes.WPARAM, wintypes.LPARAM)
        self.user32.SetWindowsHookExW.argtypes = (ctypes.c_int, self.HOOKPROC, wintypes.HINSTANCE, wintypes.DWORD)
        self.user32.SetWindowsHookExW.restype = wintypes.HHOOK
        self.user32.CallNextHookEx.argtypes = (wintypes.HHOOK, ctypes.c_int, wintypes.WPARAM, wintypes.LPARAM)
        self.user32.CallNextHookEx.restype = ctypes.c_ssize_t
        self.user32.UnhookWindowsHookEx.argtypes = (wintypes.HHOOK,)
        self.user32.GetAsyncKeyState.argtypes = (ctypes.c_int,)
        self.user32.GetAsyncKeyState.restype = ctypes.c_short
        self.kernel32.GetModuleHandleW.argtypes = (wintypes.LPCWSTR,)
        self.kernel32.GetModuleHandleW.restype = wintypes.HMODULE
        self._proc = self.HOOKPROC(self._on_key)
        app.applicationStateChanged.connect(self._state)
        app.aboutToQuit.connect(self.release)
        self._state(app.applicationState())

    def claim(self, text, callback):
        self.claims[parse(text)] = callback

    @property
    def active(self):
        return self._hook is not None

    def _state(self, state):
        if state == QtCore.Qt.ApplicationState.ApplicationActive:
            self.install()
        else:
            self.release()

    def install(self):
        if self._hook is None and self._proc is not None:
            self._hook = self.user32.SetWindowsHookExW(WH_KEYBOARD_LL, self._proc, self.kernel32.GetModuleHandleW(None), 0) or None

    def release(self):
        if self._hook is not None:
            self.user32.UnhookWindowsHookEx(self._hook)
            self._hook = None

    def _down(self, vk):
        return bool(self.user32.GetAsyncKeyState(vk) & 0x8000)

    def _on_key(self, code, wparam, lparam):
        # must answer at once: Windows drops slow hooks. The action runs afterwards, from Qt's event loop.
        if code >= 0 and wparam in (WM_KEYDOWN, WM_SYSKEYDOWN):
            k = ctypes.cast(lparam, ctypes.POINTER(KBDLLHOOKSTRUCT)).contents
            cb = matches(self.claims, k.vkCode, self._down(VK_CONTROL), self._down(VK_SHIFT), self._down(VK_MENU),
                         self._down(VK_LWIN) or self._down(VK_RWIN))
            if cb is not None:
                QtCore.QTimer.singleShot(0, cb)
                return 1                                    # kept from the program that registered it
        return self.user32.CallNextHookEx(self._hook, code, wparam, lparam)
