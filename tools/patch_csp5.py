"""tools/patch_csp5.py - patch the CSP5 0.3.1 source (MIT) so it installs on Windows (see requirements-nmr.txt):
skip the optional native peak-matching helpers (gcc-only) and let the package import without them. Prediction itself
does not use them.

    python tools\\patch_csp5.py path\\to\\csp5-0.3.1
"""
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
setup = root / 'setup.py'
s = setup.read_text()
old = '''        if os.name == "nt":
            raise RuntimeError("Building CSP5 native matching backends on Windows is unsupported")'''
new = '''        if os.name == "nt":                    # patched for Structure editor: prediction only, no peak matching
            self.announce("Skipping native matching backends on Windows", level=2)
            return'''
assert old in s or new in s
setup.write_text(s.replace(old, new))
(root / 'src' / 'csp5' / '_native' / '__init__.py').write_text('''"""Native matching backends bundled with CSP5 (patched: optional when not built, e.g. on Windows)."""

try:
    from .dp_backend import match_indices_dp
    from .murty_backend import murty_k_best
except RuntimeError:                           # not built on this platform: peak matching unavailable
    def match_indices_dp(*_a, **_k):
        raise RuntimeError("CSP5 peak matching is not available on this platform")

    def murty_k_best(*_a, **_k):
        raise RuntimeError("CSP5 peak matching is not available on this platform")

__all__ = ["match_indices_dp", "murty_k_best"]
''')
print('patched')
