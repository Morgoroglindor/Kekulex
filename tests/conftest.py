import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
if sys.platform == 'win32':
    os.environ.setdefault('QT_QPA_FONTDIR', os.path.join(os.environ.get('WINDIR', r'C:\Windows'), 'Fonts'))


@pytest.fixture(scope='session')
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(['test'])


@pytest.fixture(autouse=True)
def _settings_file(tmp_path, monkeypatch):
    monkeypatch.setenv('SED_SETTINGS_FILE', str(tmp_path / 'settings.ini'))
