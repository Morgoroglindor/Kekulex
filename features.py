"""
features.py - the parts of the program chosen in the installer. Each optional part installs a small marker file in
_internal\\features (see installer.iss); running from the sources, everything counts as installed.
"""
import os
import sys

NMR = 'nmr'
VIEW_3D = '3d'


def installed(name):
    base = getattr(sys, '_MEIPASS', None)
    return base is None or os.path.exists(os.path.join(base, 'features', name + '.txt'))
