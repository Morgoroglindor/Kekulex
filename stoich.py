"""
stoich.py - the stoichiometry grid: which compounds a reaction has, what is typed into the grid, and what follows.

    cols = stoich.reaction_columns(doc)            # (reactant atom sets, product atom sets) from the arrow
    table = stoich.compute(doc, grid)              # {'reactants': [Column], 'products': [Column]}
    stoich.parse(text, 'mass')                     # '250 mg' -> 0.25 (g)

Typed values (shown bold) and calculated ones (normal), as in a lab notebook:
  reactants  Limiting?, Equivalents, Sample Mass (g), %Weight (purity, 100 % if empty), Molarity (M), Density (g/mL),
             Volume (mL; with a molarity or a density), Reactant Moles (mol); Reactant Mass is calculated.
  products   Equivalents, %Completion, Measured Mass (g), Purity (%); Expected Mass / Moles, Product Mass / Moles and
             %Yield are calculated.
The limiting reactant's moles come from what was typed for it (moles; else volume x molarity; else volume x density
x purity / MW; else mass x purity / MW); every other reactant follows from its equivalents, unless an amount was typed
for it too (then its equivalents follow). Expected moles of a product = limiting moles x equivalents x completion.
"""
import re
from dataclasses import dataclass, field

import chem

REACTANT_ROWS = (('formula', 'Formula'), ('mw', 'MW'), ('limiting', 'Limiting?'), ('equiv', 'Equivalents'),
                 ('mass', 'Sample Mass'), ('pct', '%Weight'), ('molarity', 'Molarity'), ('density', 'Density'),
                 ('volume', 'Volume'), ('moles', 'Reactant Moles'), ('rmass', 'Reactant Mass'))
PRODUCT_ROWS = (('formula', 'Formula'), ('mw', 'MW'), ('equiv', 'Equivalents'), ('completion', '%Completion'),
                ('exp_mass', 'Expected Mass'), ('exp_moles', 'Expected Moles'), ('measured', 'Measured Mass'),
                ('purity', 'Purity'), ('pmass', 'Product Mass'), ('pmoles', 'Product Moles'), ('yield', '%Yield'))
REACTANT_INPUTS = ('equiv', 'mass', 'pct', 'molarity', 'density', 'volume', 'moles')
PRODUCT_INPUTS = ('equiv', 'completion', 'measured', 'purity')
# what kind of quantity each row is (for units and parsing)
KIND = {'mass': 'mass', 'rmass': 'mass', 'exp_mass': 'mass', 'measured': 'mass', 'pmass': 'mass',
        'moles': 'moles', 'exp_moles': 'moles', 'pmoles': 'moles', 'volume': 'volume', 'molarity': 'molarity',
        'density': 'density', 'pct': 'percent', 'completion': 'percent', 'purity': 'percent', 'yield': 'percent',
        'equiv': 'number', 'mw': 'number'}
UNITS = {'mass': {'g': 1.0, 'mg': 1e-3, 'µg': 1e-6, 'ug': 1e-6, 'kg': 1e3},
         'moles': {'mol': 1.0, 'mmol': 1e-3, 'µmol': 1e-6, 'umol': 1e-6, 'kmol': 1e3},
         'volume': {'ml': 1.0, 'l': 1e3, 'µl': 1e-3, 'ul': 1e-3},
         'molarity': {'m': 1.0, 'mm': 1e-3, 'µm': 1e-6, 'um': 1e-6, 'mol/l': 1.0},
         'density': {'g/ml': 1.0, 'g/cm3': 1.0, 'kg/l': 1.0, 'g/l': 1e-3},
         'percent': {'%': 1.0, '': 1.0}, 'number': {'': 1.0, 'eq': 1.0, 'equiv': 1.0}}


@dataclass
class Cell:
    value: object = None                                    # float, str (formula) or bool (limiting); None: blank
    typed: bool = False                                     # typed by the user (drawn bold)


@dataclass
class Column:
    atoms: list
    formula: str = ''
    formula_html: str = ''
    mw: float = 0.0
    cells: dict = field(default_factory=dict)               # row key -> Cell


def parse(text, kind):
    """A typed value in the row's base unit (g, mol, mL, M, g/mL, %), or None. Accepts '250 mg', '1,5mmol', '2 mL',
    '0.5 M', '95%'. Raises ValueError for text that is not a number with a known unit."""
    text = (text or '').strip()
    if not text:
        return None
    m = re.fullmatch(r'([-+]?\d+(?:[.,]\d*)?(?:[eE][-+]?\d+)?|[.,]\d+)\s*([^\d\s].*)?', text)
    if not m:
        raise ValueError(f'"{text}" is not a number')
    number = float(m.group(1).replace(',', '.'))
    unit = (m.group(2) or '').strip()
    units = UNITS[kind]
    if kind == 'molarity':                                  # M and mM differ only by case
        key = {'M': 'm', 'mM': 'mm', 'µM': 'µm', 'uM': 'um'}.get(unit, unit.lower())
    else:
        key = unit.lower().replace('μ', 'µ')
    if key not in units:
        raise ValueError(f'unknown unit "{unit}" (use {", ".join(u for u in units if u)})')
    return number * units[key]


def reaction_columns(doc, arrow=None):
    """The reaction around an arrow (the given one, else the first straight arrow): ([reactant atom lists],
    [product atom lists]) ordered left to right. Structures above and below the arrow (reagents, catalysts) count
    as reactants, structures right of the head as products."""
    arrows = [i for i, a in doc.arrows.items() if not a.kind.startswith('curved')]
    if arrow is None:
        if not arrows:
            return None
        arrow = arrows[0]
    ar = doc.arrows[arrow]
    lo, hi = sorted((ar.x1, ar.x2))
    forward = ar.x2 >= ar.x1
    reactants, products = [], []
    for frag in doc.fragments():
        cx = sum(doc.atoms[a].x for a in frag) / len(frag)
        before = cx < lo if forward else cx > hi
        after = cx > hi if forward else cx < lo
        if after:
            products.append((cx, frag))
        elif before or lo <= cx <= hi:
            reactants.append((cx, frag))
    reactants.sort(); products.sort()
    return [f for _x, f in reactants], [f for _x, f in products]


def _num(col, key):
    v = col.get(key)
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


def compute(doc, grid):
    """{'reactants': [Column], 'products': [Column], 'limiting_moles': float or None} for a Grid."""
    out = {'reactants': [], 'products': [], 'limiting_moles': None}
    reac = [c for c in grid.reactants if any(a in doc.atoms for a in c.get('atoms', []))]
    prod = [c for c in grid.products if any(a in doc.atoms for a in c.get('atoms', []))]

    def column(c):
        atoms = [a for a in c['atoms'] if a in doc.atoms]
        info = chem.analyse(doc, atoms)
        col = Column(atoms, info.formula, info.formula_html, info.average)
        col.cells['formula'] = Cell(info.formula)
        col.cells['mw'] = Cell(info.average)
        return col

    def own_moles(c, mw):
        w = (_num(c, 'pct') or 100.0) / 100
        if _num(c, 'moles'):
            return _num(c, 'moles')
        vol = _num(c, 'volume')
        if vol and _num(c, 'molarity'):
            return vol / 1000 * _num(c, 'molarity')
        if vol and _num(c, 'density') and mw:
            return vol * _num(c, 'density') * w / mw
        if _num(c, 'mass') and mw:
            return _num(c, 'mass') * w / mw
        return None
    lim = next((k for k, c in enumerate(reac) if c.get('limiting')), 0)
    cols = [column(c) for c in reac]
    lim_moles = lim_eq = None
    if reac:
        lim_eq = _num(reac[lim], 'equiv') or 1.0
        lim_moles = own_moles(reac[lim], cols[lim].mw)
    out['limiting_moles'] = lim_moles
    for k, (c, col) in enumerate(zip(reac, cols)):
        mw = col.mw
        w = (_num(c, 'pct') or 100.0) / 100
        col.cells['limiting'] = Cell(k == lim, typed=bool(c.get('limiting')))
        for key in REACTANT_INPUTS:                          # what was typed
            if _num(c, key) is not None:
                col.cells[key] = Cell(_num(c, key), typed=True)
        moles = own_moles(c, mw)
        eq = _num(c, 'equiv')
        if k != lim and moles is None and lim_moles is not None:
            moles = lim_moles * (eq or 1.0) / lim_eq
        if k != lim and moles is not None and eq is None and lim_moles:
            col.cells['equiv'] = Cell(moles / lim_moles * lim_eq)
        elif eq is None and lim_moles is not None:
            col.cells['equiv'] = Cell(1.0 if k == lim else lim_eq)
        if lim_moles is not None and 'pct' not in col.cells:
            col.cells['pct'] = Cell(100.0)
        if moles is not None:
            col.cells.setdefault('moles', Cell(moles))
            pure = moles * mw
            col.cells['rmass'] = Cell(pure)
            col.cells.setdefault('mass', Cell(pure / w))
            if _num(c, 'volume') is None:
                if _num(c, 'molarity'):
                    col.cells['volume'] = Cell(moles / _num(c, 'molarity') * 1000)
                elif _num(c, 'density'):
                    col.cells['volume'] = Cell(pure / w / _num(c, 'density'))
        out['reactants'].append(col)
    for c in prod:
        col = column(c)
        mw = col.mw
        for key in PRODUCT_INPUTS:
            if _num(c, key) is not None:
                col.cells[key] = Cell(_num(c, key), typed=True)
        eq = _num(c, 'equiv') or 1.0
        completion = _num(c, 'completion')
        if lim_moles is not None:
            col.cells.setdefault('equiv', Cell(1.0))
            col.cells.setdefault('completion', Cell(100.0))
            exp_moles = lim_moles * eq / (lim_eq or 1.0) * (completion if completion is not None else 100.0) / 100
            col.cells['exp_moles'] = Cell(exp_moles)
            col.cells['exp_mass'] = Cell(exp_moles * mw)
        measured = _num(c, 'measured')
        if measured is not None:
            purity = _num(c, 'purity')
            col.cells.setdefault('purity', Cell(100.0))
            pmass = measured * (purity if purity is not None else 100.0) / 100
            col.cells['pmass'] = Cell(pmass)
            if mw:
                col.cells['pmoles'] = Cell(pmass / mw)
                if lim_moles:
                    col.cells['yield'] = Cell(pmass / mw / col.cells['exp_moles'].value * 100
                                              if col.cells['exp_moles'].value else None)
        out['products'].append(col)
    return out


def display(key, cell, decimal=','):
    """(number text, unit text) of a cell as the grid shows it: masses in g / mg / µg, moles in mol / mmol / µmol,
    volumes in mL / µL, two decimals, the locale's decimal separator."""
    v = cell.value
    if v is None:
        return '', ''
    if key == 'limiting':
        return ('Yes' if v else 'No'), ''
    if key == 'formula':
        return str(v), ''
    kind = KIND.get(key, 'number')
    unit = ''
    if kind == 'mass':
        v, unit = (v, 'g') if abs(v) >= 1 or v == 0 else (v * 1e3, 'mg') if abs(v) >= 1e-3 else (v * 1e6, 'µg')
    elif kind == 'moles':
        v, unit = (v, 'mol') if abs(v) >= 1 or v == 0 else (v * 1e3, 'mmol') if abs(v) >= 1e-3 else (v * 1e6, 'µmol')
    elif kind == 'volume':
        v, unit = (v, 'mL') if abs(v) >= 1 or v == 0 else (v * 1e3, 'µL')
    elif kind == 'molarity':
        unit = 'M'
    elif kind == 'density':
        unit = 'g/mL'
    elif kind == 'percent':
        unit = '%'
    return f'{v:.2f}'.replace('.', decimal), unit


def new_grid_columns(doc, arrow=None):
    """Grid columns for the reaction: ([{'atoms': ..., 'limiting': first}], [{'atoms': ...}]) or None."""
    found = reaction_columns(doc, arrow)
    if not found or not found[0]:
        return None
    reac, prod = found
    return ([dict(atoms=list(f), limiting=(k == 0)) for k, f in enumerate(reac)], [dict(atoms=list(f)) for f in prod])
