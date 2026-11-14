"""
msfile.py - mass spectra from files and pasted text, for annotating them with the drawn structures:

    spectra = msfile.read(path)            # [Spectrum]: one for CSV / TXT / JCAMP-DX, one per scan in mzML / MGF / MSP
    spectra = msfile.parse_text(text)      # pasted 'm/z  intensity' lines
    peaks = msfile.centroid(spectrum)      # [(m/z, intensity)]: profile data reduced to its peak tops

Readers: two-column text (comma, tab, semicolon or spaces; header lines skipped), JCAMP-DX (##PEAK TABLE / ##XYDATA
as X,Y pairs), NIST MSP ("Num Peaks:" then pairs), MGF (BEGIN IONS ... END IONS) and mzML (32 / 64-bit arrays, zlib
or none). An mzML file's scans come with their retention time, MS level and total ion current.
"""
import base64
import os
import re
import zlib
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Spectrum:
    title: str
    mz: np.ndarray
    intensity: np.ndarray
    rt: float | None = None                                 # minutes
    level: int = 1
    precursor: float | None = None
    meta: dict = field(default_factory=dict)

    @property
    def tic(self):
        return float(self.intensity.sum()) if len(self.intensity) else 0.0


_NUM = re.compile(r'[-+]?(?:\d+[.,]?\d*|[.,]\d+)(?:[eE][-+]?\d+)?')


def _pairs(lines):
    """(m/z, intensity) arrays from lines holding two numbers (more are ignored; a decimal comma is read when the
    separator is not a comma)."""
    mz, inten = [], []
    for line in lines:
        line = line.strip()
        if not line or line[0] in '#;%!' or line.startswith('##'):
            continue
        if ';' in line or '\t' in line or ' ' in line.strip():
            parts = re.split(r'[;\t ]+', line.strip())
            nums = [p.replace(',', '.') for p in parts if _NUM.fullmatch(p)]
        else:
            nums = [p for p in line.split(',') if _NUM.fullmatch(p.strip())]
        if len(nums) < 2:
            continue
        try:
            x, y = float(nums[0]), float(nums[1])
        except ValueError:
            continue
        mz.append(x); inten.append(y)
    return _sorted(mz, inten)


def _sorted(mz, inten):
    mz, inten = np.asarray(mz, float), np.asarray(inten, float)
    order = np.argsort(mz, kind='stable')
    return mz[order], inten[order]


def parse_text(text, title='Pasted spectrum'):
    """Spectra in text: JCAMP-DX, MSP or MGF when it looks like one, else 'm/z intensity' lines."""
    head = text.lstrip()[:4000].upper()
    if head.startswith('##') or '##TITLE' in head:
        return _jcamp(text, title)
    if 'BEGIN IONS' in head:
        return _mgf(text)
    if 'NUM PEAKS' in head:
        return _msp(text)
    mz, inten = _pairs(text.splitlines())
    return [Spectrum(title, mz, inten)] if len(mz) else []


def _jcamp(text, title):
    out, cur, name = [], [], title
    reading = False
    for line in text.splitlines():
        s = line.strip()
        if s.upper().startswith('##TITLE='):
            name = s[8:].strip() or title
        elif s.upper().startswith(('##PEAK TABLE', '##XYPOINTS', '##XYDATA', '##DATA TABLE')):
            reading = True
            cur = []
        elif s.startswith('##'):
            if reading and cur:
                mz, inten = _sorted(*zip(*cur)) if cur else ([], [])
                out.append(Spectrum(name, mz, inten))
            reading = False
        elif reading:
            nums = [float(v) for v in re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', s)]   # "x,y x,y" / "x y; x y"
            cur += [(nums[k], nums[k + 1]) for k in range(0, len(nums) - 1, 2)]
    if reading and cur:
        out.append(Spectrum(name, *_sorted(*zip(*cur))))
    return out


def _msp(text):
    out, name, pairs, want = [], 'Spectrum', [], 0
    for line in text.splitlines() + ['']:
        s = line.strip()
        low = s.lower()
        if low.startswith('name:'):
            name = s[5:].strip()
        elif low.startswith('num peaks:'):
            want, pairs = int(re.findall(r'\d+', s)[0]), []
        elif want and s:
            nums = re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', s)
            pairs += [(float(nums[k]), float(nums[k + 1])) for k in range(0, len(nums) - 1, 2)]
            if len(pairs) >= want:
                out.append(Spectrum(name, *_sorted(*zip(*pairs[:want]))))
                want = 0
        elif not s and want and pairs:
            out.append(Spectrum(name, *_sorted(*zip(*pairs))))
            want = 0
    return out


def _mgf(text):
    out, meta, pairs, inside = [], {}, [], False
    for line in text.splitlines():
        s = line.strip()
        if s.upper() == 'BEGIN IONS':
            meta, pairs, inside = {}, [], True
        elif s.upper() == 'END IONS':
            if pairs:
                rt = meta.get('RTINSECONDS')
                pep = meta.get('PEPMASS', '').split()
                out.append(Spectrum(meta.get('TITLE', f'Spectrum {len(out) + 1}'), *_sorted(*zip(*pairs)),
                                    rt=float(rt) / 60 if rt and _NUM.fullmatch(rt) else None, level=2,
                                    precursor=float(pep[0]) if pep and _NUM.fullmatch(pep[0]) else None, meta=meta))
            inside = False
        elif inside and '=' in s and not s[0].isdigit():
            k, v = s.split('=', 1)
            meta[k.upper()] = v
        elif inside and s:
            nums = s.split()
            if len(nums) >= 2 and _NUM.fullmatch(nums[0]) and _NUM.fullmatch(nums[1]):
                pairs.append((float(nums[0]), float(nums[1])))
    return out


# ------------------------------------------------------------------------------------------------------------- mzML
_MS = '{http://psi.hupo.org/ms/mzml}'


def _decode(binary, params):
    data = base64.b64decode(binary or '')
    if 'MS:1000574' in params:                               # zlib compression
        data = zlib.decompress(data)
    if 'MS:1000521' in params:                               # 32-bit float
        return np.frombuffer(data, '<f4').astype(float)
    return np.frombuffer(data, '<f8').astype(float)


def _mzml(path, most=20000):
    import xml.etree.ElementTree as ET
    out = []
    for _ev, el in ET.iterparse(path, events=('end',)):
        tag = el.tag.rsplit('}', 1)[-1]
        if tag != 'spectrum':
            continue
        cv = {c.get('accession'): c.get('value') for c in el.iter() if c.tag.endswith('cvParam')}
        level = int(cv.get('MS:1000511') or 1)
        rt = None
        for c in el.iter():
            if c.tag.endswith('cvParam') and c.get('accession') == 'MS:1000016':
                rt = float(c.get('value'))
                if c.get('unitName', '') in ('second', '') and c.get('unitAccession') != 'UO:0000031':
                    rt /= 60.0
        arrays = {}
        for arr in el.iter():
            if not arr.tag.endswith('binaryDataArray'):
                continue
            params = {c.get('accession') for c in arr.iter() if c.tag.endswith('cvParam')}
            b = next((x.text for x in arr.iter() if x.tag.endswith('binary')), '')
            kind = 'mz' if 'MS:1000514' in params else 'int' if 'MS:1000515' in params else None
            if kind:
                arrays[kind] = _decode(b, params)
        prec = None
        for c in el.iter():
            if c.tag.endswith('cvParam') and c.get('accession') == 'MS:1000744':
                prec = float(c.get('value'))
        if 'mz' in arrays and 'int' in arrays and len(arrays['mz']) == len(arrays['int']):
            sid = el.get('id', '')
            scan = re.search(r'scan=(\d+)', sid)
            title = f'scan {scan.group(1) if scan else el.get("index", len(out))}'
            out.append(Spectrum(title, *_sorted(arrays['mz'], arrays['int']), rt=rt, level=level, precursor=prec,
                                meta={'profile': 'MS:1000128' in cv}))
        el.clear()
        if len(out) >= most:
            break
    return out


def read(path):
    """The spectra in a file ([] when it holds none that can be read)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.mzml':
        return _mzml(path)
    with open(path, 'rb') as f:
        raw = f.read()
    for enc in ('utf-8-sig', 'utf-16', 'latin-1'):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    if ext in ('.jdx', '.dx', '.jcamp'):
        return _jcamp(text, os.path.basename(path))
    if ext == '.mgf':
        return _mgf(text)
    if ext == '.msp':
        return _msp(text)
    return parse_text(text, os.path.basename(path))


# ---------------------------------------------------------------------------------------------------------- peaks
def is_profile(sp):
    """True for profile (continuum) data: many points closely and evenly spaced."""
    if sp.meta.get('profile'):
        return True
    if len(sp.mz) < 200:
        return False
    return float(np.median(np.diff(sp.mz))) < 0.05


def centroid(sp, min_rel=0.0):
    """[(m/z, intensity)] of the spectrum's peaks: its points as they are (centroided data), or the top of each
    profile peak (a Gaussian through its three highest points)."""
    mz, y = sp.mz, sp.intensity
    if not len(mz):
        return []
    if not is_profile(sp):
        pts = list(zip(mz.tolist(), y.tolist()))
    else:
        pts = []
        for k in range(1, len(y) - 1):
            if y[k] > 0 and y[k] >= y[k - 1] and y[k] > y[k + 1]:
                a, b, c = y[k - 1], y[k], y[k + 1]
                if a > 0 and c > 0:                         # a Gaussian through the top three points
                    la, lb, lc = np.log(a), np.log(b), np.log(c)
                    den = la - 2 * lb + lc
                    shift = 0.5 * (la - lc) / den if den < 0 else 0.0
                    step = (mz[k + 1] - mz[k - 1]) / 2
                    pts.append((float(mz[k] + max(-1.0, min(1.0, shift)) * step), float(b)))
                else:
                    w = y[k - 1:k + 2]
                    pts.append((float((mz[k - 1:k + 2] * w).sum() / w.sum()), float(b)))
    top = max((p[1] for p in pts), default=0.0)
    return [p for p in pts if top and p[1] >= top * min_rel]


def unit_resolution(peaks):
    """True when the m/z values are whole numbers (a low-resolution EI library spectrum): match in Da, not ppm."""
    if not peaks:
        return False
    return sum(1 for m, _y in peaks if abs(m - round(m)) < 1e-6) >= 0.9 * len(peaks)

