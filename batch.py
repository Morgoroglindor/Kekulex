"""File > Batch export: many structures from an SDF or SMILES list to pictures - one file each, or one grid sheet
with a name under each structure (for theses and supporting information)."""
import os
import re

from rdkit import Chem

import chem
import render
from model import Doc


def read(path):
    """[(name, molecule)] from .sdf (the title line, else a numbered name), .smi / .txt / .csv (SMILES, then an
    optional name), .mol or .cdx."""
    low = path.lower()
    out = []
    if low.endswith(('.smi', '.smiles', '.txt', '.csv')):
        with open(path, encoding='utf-8', errors='replace') as f:
            for k, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = re.split(r'[\t,; ]+', line, maxsplit=1)
                mol = Chem.MolFromSmiles(parts[0])
                if mol is None:
                    continue
                out.append(((parts[1].strip().strip('"') if len(parts) > 1 else '') or f'{len(out) + 1}', mol))
        if not out:
            raise ValueError('No SMILES could be read from this file.')
        return out
    for k, mol in enumerate(chem.mols_from_file(path), 1):
        name = mol.GetProp('_Name').strip() if mol.HasProp('_Name') else ''
        out.append((name or str(k), mol))
    return out


def _doc_of(mol, style, at=(0.0, 0.0)):
    d = Doc()
    d.style = dict(style or {})
    m = Chem.Mol(mol)
    m.RemoveAllConformers()                                  # always the usual 2D drawing
    d.insert_rdkit(m, at)                                   # laid out the usual way (layout2d)
    return d


def _safe(name):
    name = re.sub(r'[\\/:*?"<>|]+', '_', name).strip() or 'structure'
    return name[:80]


def _save(doc, path, dpi, transparent):
    ext = os.path.splitext(path)[1].lower()
    if ext == '.svg':
        render.save_svg(doc, path)
    elif ext == '.pdf':
        render.save_pdf(doc, path)
    else:
        render.save_png(doc, path, dpi, transparent)


def export_each(items, folder, ext='.png', style=None, names=True, dpi=600, transparent=False):
    """One picture per structure in folder (named after the structure); returns the paths written."""
    os.makedirs(folder, exist_ok=True)
    written, used = [], set()
    for name, mol in items:
        d = _doc_of(mol, style)
        if names:
            _caption(d, name)
        base = _safe(name)
        stem, k = base, 2
        while stem.lower() in used:
            stem = f'{base}_{k}'; k += 1
        used.add(stem.lower())
        path = os.path.join(folder, stem + ext)
        _save(d, path, dpi, transparent)
        written.append(path)
    return written


def _caption(d, name, below=None):
    """The name centred under the structure in d (below: the y to use)."""
    xs = [a.x for a in d.atoms.values()]
    ys = [a.y for a in d.atoms.values()]
    if not xs:
        return
    rend = render.Renderer(d)
    bottom = max((rend.item_rect('atom', a).bottom() for a in d.atoms), default=max(ys) * rend.BL) / rend.BL
    d.add_text((min(xs) + max(xs)) / 2, (below if below is not None else bottom + 0.7), name, align='center')


def grid_doc(items, columns=4, style=None, names=True):
    """One drawing with the structures in a grid (rows of columns), each with its name under it."""
    docs = [(name, _doc_of(mol, style)) for name, mol in items]
    sizes = []
    for _n, d in docs:
        xs = [a.x for a in d.atoms.values()] or [0.0]
        ys = [a.y for a in d.atoms.values()] or [0.0]
        sizes.append((max(xs) - min(xs), max(ys) - min(ys)))
    cw = max((w for w, _h in sizes), default=1.0) + 2.5
    ch = max((h for _w, h in sizes), default=1.0) + (3.2 if names else 2.0)
    out = Doc()
    out.style = dict(style or {})
    for k, (name, d) in enumerate(docs):
        row, col = divmod(k, max(1, columns))
        cx, cy = col * cw + cw / 2, row * ch + (ch - (1.2 if names else 0)) / 2
        xs = [a.x for a in d.atoms.values()] or [0.0]
        ys = [a.y for a in d.atoms.values()] or [0.0]
        dx = cx - (min(xs) + max(xs)) / 2
        dy = cy - (min(ys) + max(ys)) / 2
        out.insert_dict(d.to_dict(), dx, dy)
        if names:
            out.add_text(cx, row * ch + ch - 1.6, name, align='center')
    return out


def export_grid(items, path, columns=4, style=None, names=True, dpi=600, transparent=False):
    d = grid_doc(items, columns, style, names)
    _save(d, path, dpi, transparent)
    return [path]
