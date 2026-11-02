"""Parts chosen in the installer: without them the program still runs and their commands are gone."""
import os
import sys

import pytest
from rdkit import Chem

import features
import nmr
from model import Doc


def test_from_the_sources_everything_counts():
    assert features.installed(features.NMR) and features.installed(features.VIEW_3D)


def test_markers_in_a_frozen_program(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, '_MEIPASS', str(tmp_path), raising=False)
    assert not features.installed(features.NMR)
    (tmp_path / 'features').mkdir()
    (tmp_path / 'features' / 'nmr.txt').write_text('x')
    assert features.installed(features.NMR) and not features.installed(features.VIEW_3D)
    assert not nmr.available() or (tmp_path / 'features' / 'nmr.txt').exists()


def _menu_texts(win):
    out = []
    for top in win.menuBar().actions():
        out += [a.text() for a in top.menu().actions()]
    return out


def test_program_without_the_optional_parts(qapp, monkeypatch):
    import kekulex as se
    monkeypatch.setattr(features, 'installed', lambda name: False)
    w = se.MainWindow()
    try:
        w.show()
        qapp.processEvents()
        assert w.panel3d is None and w.tabs.count() == 1 and not w.tabs.tabBar().isVisible()
        assert not nmr.available()
        texts = ' '.join(_menu_texts(w))
        assert 'NMR' not in texts and '3D' not in texts
        d = Doc()
        d.insert_rdkit(Chem.MolFromSmiles('CCO'))
        w.canvas.set_doc(d)
        w.send_to_3d()
        w.canvas.select_all()
        w._retheme()
    finally:
        w.dirty = False
        w.close()


def test_program_with_everything(qapp):
    import kekulex as se
    w = se.MainWindow()
    try:
        assert w.tabs.count() == 2
        assert '3D' in ' '.join(_menu_texts(w))
    finally:
        w.dirty = False
        w.close()
