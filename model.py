"""
model.py - the drawing: atoms and bonds (any number of separate structures), reaction arrows and text, with the
geometry the editor needs (where the next bond goes, rings, chains, the side of a double bond, implicit hydrogens) and
the MOL conversion RDKit reads.

    doc = Doc()
    a = doc.add_atom(0, 0); b = doc.add_atom(1, 0, 'O'); doc.add_bond(a, b, 2)
    doc.to_molblock()            # V2000 MOL text (y up, 1.5 A bonds)
    Doc.from_dict(doc.to_dict()) # the .sed file is this dict as JSON

Coordinates are in bond lengths (one bond = 1.0) with y pointing down, as on the screen. Every item (atom, bond, arrow,
text, shape) has an integer id from one counter, so a selection can hold any of them as (kind, id). An atom with a
label is a nickname (OMe, Boc ...): drawn as that text, written out in full by expanded() for the chemistry.
"""
import copy
import math
import re
from dataclasses import asdict, dataclass, field

MERGE = 0.3                                                 # a new atom this close to an existing one is that atom
MOL_BOND = 1.5                                              # bond length written into MOL files (A, RDKit's own)

# Usual valences (the lowest that fits the bonds is used) and what a charge does to them
VALENCE = {'H': (1,), 'B': (3,), 'C': (4,), 'N': (3,), 'O': (2,), 'F': (1,), 'Si': (4,), 'P': (3, 5), 'S': (2, 4, 6),
           'Cl': (1, 3, 5, 7), 'Se': (2, 4, 6), 'Br': (1, 3, 5), 'I': (1, 3, 5, 7), 'As': (3, 5), 'Te': (2, 4, 6),
           'Li': (1,), 'Na': (1,), 'K': (1,), 'Mg': (2,)}
# Groups the hotkeys put on an atom, drawn out in full; the first atom of each SMILES is the one that attaches
GROUPS = {'Me': 'C', 'Et': 'CC', 'Ac': 'C(C)=O', 'CO2Me': 'C(=O)OC', 'CF3': 'C(F)(F)F', 'Cbz': 'C(=O)OCc1ccccc1',
          'MgBr': '[Mg]Br', 'NO2': '[N+](=O)[O-]', 'OMe': 'OC', 'Fmoc': 'C(=O)OCC1c2ccccc2-c2ccccc21', 'Ph': 'c1ccccc1',
          'Boc': 'C(=O)OC(C)(C)C', 'N3': 'N=[N+]=[N-]', 'tBu': 'C(C)(C)C', 'SO2': 'S(=O)=O', 'Cy': 'C1CCCCC1',
          'cPent': 'C1CCCC1', 'cPr': 'C1CC1', 'cBu': 'C1CCC1', 'Cp': 'C1C=CC=C1', 'ethynyl': 'C#C'}
# Valence electrons, for charged atoms: an ion makes min(e, 8 - e) bonds with e = electrons - charge (the octet), as
# ChemDraw counts them: NH2 -> NH3+ -> NH2 2+, C+ and C- take three bonds, O+ three, O- one, B- four.
VALENCE_ELECTRONS = {'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Se': 6, 'Br': 7,
                     'I': 7, 'As': 5, 'Te': 6, 'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2}
# Nicknames: an atom shows the short label (OMe, Boc ...) and counts as the whole group for the formula and SMILES.
FIXED_NICKNAMES = {k: v for k, v in GROUPS.items() if k not in ('SO2', 'Cy', 'cPent', 'cPr', 'cBu', 'Cp', 'ethynyl')}
FIXED_NICKNAMES.update({
    'CO2H': 'C(=O)O', 'COOH': 'C(=O)O', 'CO2Et': 'C(=O)OCC', 'CHO': 'C=O', 'CN': 'C#N', 'OAc': 'OC(C)=O',
    'NHAc': 'NC(C)=O', 'SO3H': 'S(=O)(=O)O', 'SO2Me': 'S(=O)(=O)C', 'TMS': '[Si](C)(C)C', 'TBS': '[Si](C)(C)C(C)(C)C',
    'TBDMS': '[Si](C)(C)C(C)(C)C', 'TIPS': '[Si](C(C)C)(C(C)C)C(C)C', 'Ts': 'S(=O)(=O)c1ccc(C)cc1',
    'Ms': 'S(=O)(=O)C', 'Tf': 'S(=O)(=O)C(F)(F)F', 'OTs': 'OS(=O)(=O)c1ccc(C)cc1', 'OMs': 'OS(=O)(=O)C',
    'OTf': 'OS(=O)(=O)C(F)(F)F', 'Bn': 'Cc1ccccc1', 'Bz': 'C(=O)c1ccccc1', 'OBn': 'OCc1ccccc1', 'OEt': 'OCC',
    'OPh': 'Oc1ccccc1', 'iPr': 'C(C)C', 'nPr': 'CCC', 'Pr': 'CCC', 'nBu': 'CCCC', 'Bu': 'CCCC', 'sBu': 'C(C)CC',
    'iBu': 'CC(C)C', 'NMe2': 'N(C)C', 'NEt2': 'N(CC)CC', 'NHMe': 'NC', 'NHBoc': 'NC(=O)OC(C)(C)C',
    'NHCbz': 'NC(=O)OCc1ccccc1', 'NHFmoc': 'NC(=O)OCC1c2ccccc2-c2ccccc21', 'PMB': 'Cc1ccc(OC)cc1',
    'THP': 'C1CCCCO1', 'CONH2': 'C(N)=O', 'SMe': 'SC', 'CH2OH': 'CO', 'CCl3': 'C(Cl)(Cl)Cl', 'Tr': 'C(c1ccccc1)(c1ccccc1)c1ccccc1',
    'Ph3P': 'P(c1ccccc1)(c1ccccc1)c1ccccc1', 'PPh3': 'P(c1ccccc1)(c1ccccc1)c1ccccc1', 'Piv': 'C(=O)C(C)(C)C',
    'OPiv': 'OC(=O)C(C)(C)C', 'Mes': 'c1c(C)cc(C)cc1C', 'Np': 'c1cccc2ccccc12', 'Py': 'c1ccccn1'})


class _Nicknames(dict):
    """The nicknames above, and any label that reads as a condensed formula (SO3K, CH2CH3, N(CH3)2 ...:
    condensed.py), so that typed groups count in full like the named ones. Iterating gives the named ones only."""
    def _formula(self, key):
        if not isinstance(key, str) or dict.__contains__(self, key):
            return None
        import condensed
        return condensed.smiles(key)

    def __contains__(self, key):
        return dict.__contains__(self, key) or self._formula(key) is not None

    def __missing__(self, key):
        smiles = self._formula(key)
        if smiles is None:
            raise KeyError(key)
        return smiles

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default


NICKNAMES = _Nicknames(FIXED_NICKNAMES)
LABEL_TOKEN = r'[A-Z][a-z]*\d*'


def _token_ok(tok):
    letters = tok.rstrip('0123456789')
    return letters in VALENCE or letters in NICKNAMES or letters in ('Mg', 'Li', 'Na', 'K', 'Zn', 'Cu')


def label_order(label, reverse):
    """The label's tokens in reading order for an atom bonded on its right (reverse: CO2Me -> MeO2C, NO2 -> O2N), and
    the index of the token that sits on the atom. Labels that do not split into chemical tokens (TMS, iPr) never turn."""
    import re
    sign = re.search(r'[+\-−]+$', label)                     # the count before it stays with its group (NMe3)
    if reverse and sign and sign.start() > 0:               # NMe3+ -> Me3N+: the charge stays at the end
        stem, anchor = label_order(label[:sign.start()], True)
        if len(stem) > 1:
            return stem + [sign.group(0)], anchor
    toks = re.findall(LABEL_TOKEN + r'|[^A-Z]+', label)
    if reverse and len(toks) > 1 and ''.join(toks) == label and all(re.fullmatch(LABEL_TOKEN, t) and _token_ok(t) for t in toks):
        return toks[::-1], len(toks) - 1
    whole =_nickname_tokens(label)                         # OTBS -> TBSO, NHBoc -> BocHN: a nickname turns as one
    if reverse and whole and len(whole) > 1:
        return whole[::-1], len(whole) - 1
    if reverse and '(' in label:                            # N(CH3)2 -> (H3C)2N: a bracket turns as one piece
        units = _label_units(re.findall(LABEL_TOKEN + r'|\(|\)\d*', label))
        if units and ''.join(_flat(units)) == label:
            out = _flat(_turned(units))
            heavy = [k for k, t in enumerate(out) if re.fullmatch(LABEL_TOKEN, t) and t.rstrip('0123456789') != 'H']
            return out, heavy[-1]
    return toks, 0


def _nickname_tokens(label):
    """label split into element tokens and whole nicknames (TBS, Boc, Ph ...), or None when it does not split so."""
    import re
    names = sorted((k for k in FIXED_NICKNAMES if k.isalpha() and len(k) > 1 and k[0].isupper()), key=len, reverse=True)
    out, i = [], 0
    while i < len(label):
        nick = next((n for n in names if label.startswith(n, i) and not label[i + len(n):i + len(n) + 1].islower()), None)
        if nick:
            m = re.match(r'\d*', label[i + len(nick):])
            out.append(nick + m.group(0)); i += len(nick) + len(m.group(0))
            continue
        m = re.match(LABEL_TOKEN, label[i:])
        if not m or not _token_ok(m.group(0)):
            return None
        out.append(m.group(0)); i += len(m.group(0))
    return out if any(t.rstrip('0123456789') in FIXED_NICKNAMES and len(t.rstrip('0123456789')) > 1 for t in out) else None


def _label_units(toks):
    """Tokens grouped by brackets: ['N', ['(', 'C', 'H3', ')2']]; None if they do not all read as chemistry."""
    import re
    stack = [[]]
    for t in toks:
        if t == '(':
            stack.append(['('])
        elif t.startswith(')'):
            if len(stack) < 2:
                return None
            group = stack.pop() + [t]
            stack[-1].append(group)
        elif re.fullmatch(LABEL_TOKEN, t) and _token_ok(t):
            stack[-1].append(t)
        else:
            return None
    return stack[0] if len(stack) == 1 and len(stack[0]) > 1 else None


def _turned(units):
    out = []
    for u in reversed(units):
        out.append(['('] + _turned(u[1:-1]) + [u[-1]] if isinstance(u, list) else u)
    return out


def _flat(units):
    return [t for u in units for t in (_flat(u) if isinstance(u, list) else [u])]


def parse_label(text):
    """What a typed atom label means: ('element', symbol, charge) - also 'OH', 'NH2', 'Na+', 'O-' -, ('nickname', label)
    or None."""
    import re
    from rdkit import Chem
    text = text.strip().replace('−', '-')
    if text in NICKNAMES:
        return ('nickname', text)
    low = {k.lower(): k for k in NICKNAMES}                 # the fixed nicknames: boc is Boc, not B-O-C
    if text != text.upper() and text.lower() in low:
        return ('nickname', low[text.lower()])
    guess = case_guess(text)                                # all lower case: nh2 is NH2, not nihonium
    if guess and guess != text:
        what = parse_label(guess)
        if what is not None:
            return what
    m = _ELEMENT_LABEL.fullmatch(text)
    if m:
        sym = m.group(1)[0].upper() + m.group(1)[1:].lower()
        ok = sym in _element_symbols()
        if ok:
            signs = m.group(5)
            n = int(m.group(4)) if m.group(4) else len(signs)
            after = label_charge(text)                      # O•+: the sign after the radical dot
            if after is not None and not signs:
                return ('element', sym, after)
            return ('element', sym, 0 if not signs else (n if signs[0] == '+' else -n))
    low = {k.lower(): k for k in NICKNAMES}
    if text.lower() in low:
        return ('nickname', low[text.lower()])
    return None


RADICAL_MARKS = '•·.*'                                     # typed after a label: CH3•, CH2••, O•, N•+ (or . / *)
_ELEMENT_LABEL = re.compile(r'([A-Za-z][a-z]?)(H(\d*))?\s*(\d*)\s*([+\-]*)\s*([•·.*]*)\s*(?:(\d*)\s*([+\-]*))?')


def label_extras(text):
    """(hydrogens typed, or None; unpaired electrons) of a typed element label: 'CH2' -> (2, 0), 'CH3•' -> (3, 1),
    'C••' -> (None, 2), 'O•+' -> (None, 1), 'N' -> (None, 0). Nicknames and formulas: (None, 0)."""
    text = text.strip().replace('−', '-')
    guess = case_guess(text)
    if guess:
        text = guess
    m = _ELEMENT_LABEL.fullmatch(text)
    if not m or (m.group(1)[0].upper() + m.group(1)[1:].lower()) not in _element_symbols():
        return None, 0
    hs = None if m.group(2) is None else int(m.group(3) or 1)
    return hs, min(len(m.group(6)), 2)


def label_charge(text):
    """The charge of a typed element label whose sign follows the radical dot ('O•+' -> 1), else None."""
    m = _ELEMENT_LABEL.fullmatch(text.strip().replace('−', '-'))
    if not m or not m.group(8):
        return None
    n = int(m.group(7)) if m.group(7) else len(m.group(8))
    return n if m.group(8)[0] == '+' else -n


# Two lower-case letters read as one symbol when a label is typed all in lower case (so3k, ch2cl, co2et); any other
# pair is two one-letter elements, so co is C + O, never cobalt
_LOWER_PAIRS = ('cl', 'br', 'na', 'li', 'mg', 'si', 'me', 'et', 'ph', 'ac', 'bn', 'bz', 'bu', 'pr', 'ts', 'tf', 'cy')


def _element_symbols():
    global _SYMBOLS
    if _SYMBOLS is None:
        from rdkit import Chem
        table = Chem.GetPeriodicTable()
        _SYMBOLS = {table.GetElementSymbol(z) for z in range(1, 119)}
    return _SYMBOLS


_SYMBOLS = None


def case_guess(text):
    """An all-lower-case label in capitals ChemDraw-style ('so3k' -> 'SO3K', 'ch2cl' -> 'CH2Cl', 'ome' -> 'OMe'),
    or None when it has capitals already or cannot be read that way."""
    import re
    if not text or text != text.lower() or not re.search('[a-z]', text):
        return None
    element = _element_symbols().__contains__
    out, k = [], 0
    while k < len(text):
        ch = text[k]
        if not ch.isalpha():
            out.append(ch); k += 1
            continue
        pair = text[k:k + 2]
        if len(pair) == 2 and pair.isalpha() and pair in _LOWER_PAIRS:
            out.append(pair.capitalize()); k += 2
        elif element(ch.upper()):
            out.append(ch.upper()); k += 1
        elif len(pair) == 2 and pair.isalpha() and (element(pair.capitalize()) or pair.capitalize() in NICKNAMES):
            out.append(pair.capitalize()); k += 2
        else:
            return None
    return ''.join(out)


@dataclass
class Atom:
    x: float
    y: float
    el: str = 'C'
    charge: int = 0
    label: str = ''                                         # a nickname (OMe, Boc ...) shown instead of the element
    color: str = ''                                         # '' = the drawing's colour
    highlight: str = ''                                     # a coloured disc behind the atom ('' = none)
    charge_pos: list | None = None                          # where its charge sign was dragged: [dx, dy] from the
                                                            # atom in bond lengths; None = placed automatically
    radical: int = 0                                        # unpaired electrons: 1 radical (•), 2 carbene / diradical
    lone_pairs: int = 0                                     # lone pairs drawn as dot pairs: 0 none, -1 all it has,
                                                            # n: that many
    map: int = 0                                            # reaction atom-map number (0: none)
    hcount: int | None = None                               # hydrogens typed with the label (CH2, NH): fixed, whatever
                                                            # the bonds; None = the usual number for the valence


@dataclass
class Bond:
    a: int
    b: int
    order: int = 1                                          # 1, 2, 3
    stereo: str = ''                                        # '' | 'wedge' | 'hash' (the narrow end is atom a)
    side: int | None = None                                 # double bond: None = automatic, 0 = centred, +1 / -1 = side
    style: str = ''                                         # '' | 'bold' | 'dashed' | 'wavy' | 'hashed' (single bonds;
    color: str = ''                                         # 'dashed' on a double bond dashes its second line)
    highlight: str = ''                                     # a coloured band behind the bond ('' = none)


ARROW_KINDS = ('forward', 'equilibrium', 'retro', 'resonance', 'dashed', 'curved', 'curved_half')


@dataclass
class Arrow:
    x1: float                                               # tail
    y1: float
    x2: float                                               # head
    y2: float
    kind: str = 'forward'                                   # one of ARROW_KINDS
    bend: float = 0.0                                       # curved arrows: bulge, as a fraction of the length (+/- side)
    color: str = ''
    tail_on: list | None = None                             # electron-pushing arrows: [kind, id, dx, dy] - the tail
    head_on: list | None = None                             # (head) stays at that atom / bond plus the offset


@dataclass
class Text:
    x: float                                                # left edge, centre or right edge (align); top of the text
    y: float
    text: str
    align: str = 'left'                                     # 'left' | 'center' | 'right'
    size: float = 0.0                                       # pt; 0 = the drawing's font size
    bold: bool = False
    italic: bool = False
    formula: bool = False                                   # numbers after letters drawn as subscripts (H2SO4)
    color: str = ''
    on_arrow: list | None = None                            # reaction conditions: [arrow id, 'above' | 'below']
    live: list = field(default_factory=list)                # a pasted analysis: the property keys it shows, kept up to
                                                            # date with the structure in its group (chem.live_caption)


SHAPE_KINDS = ('rect', 'roundrect', 'oval', 'line', 'brackets', 'parens', 'braces')
BRACKET_KINDS = ('brackets', 'parens', 'braces')             # these can be polymer brackets (params['repeat'])
# ChemDraw's other drawing elements and the orbitals (elements.py draws them): a circle (centre x1,y1, edge x2,y2),
# an arc (chord x1,y1 -> x2,y2, params['angle'] 90 / 120 / 180 / 270), an orbital (node x1,y1, tip x2,y2,
# params['type'] s / sigma / lobe / p / hybrid / dxy / dz2)
ELEMENT_KINDS = ('circle', 'arc', 'orbital')
SHAPE_STYLES = ('', 'dashed', 'shaded', 'filled', 'shadow', 'bold', 'wavy')
# BioDraw-like objects (drawn by render.bio_paint in their box; the "linear" ones along the line x1,y1 -> x2,y2)
BIO_KINDS = ('membrane', 'membrane_arc', 'micelle', 'dna', 'helix', 'receptor', 'channel', 'enzyme', 'gprotein',
             'antibody', 'ribosome', 'trna', 'mrna', 'mito', 'nucleus', 'cell', 'bacterium', 'virus', 'cloud',
             'plasmid')
BIO_LINEAR = ('membrane', 'dna', 'helix', 'mrna')


@dataclass
class Shape:
    x1: float
    y1: float
    x2: float
    y2: float
    kind: str = 'rect'                                      # SHAPE_KINDS, ELEMENT_KINDS or BIO_KINDS (brackets in pairs)
    color: str = ''
    fill: str = ''                                          # bio objects: main fill colour ('' = the kind's default)
    params: dict = field(default_factory=dict)              # bio objects: plasmid {'bp', 'regions', 'markers'} ...
    style: str = ''                                         # one of SHAPE_STYLES ('' plain; lines: dashed / bold /
                                                            # wavy; orbitals: '' / shaded / filled phase)


@dataclass
class TLC:
    """A thin-layer chromatography plate (tlc.py draws it): lanes of spots placed by Rf between the origin line and
    the solvent front."""
    x: float                                                # top-left corner, bond lengths
    y: float
    w: float = 3.0
    h: float = 4.0
    origin: float = 0.12                                    # origin line, fraction of the height from the bottom
    front: float = 0.1                                      # solvent front, fraction of the height from the top
    lanes: list = field(default_factory=list)               # [[{'rf', 'w', 'h', 'style', 'color', 'show_rf'}]]
    show_origin: bool = True
    show_front: bool = True
    show_border: bool = True
    show_ticks: bool = True
    color: str = ''


@dataclass
class Sequence:
    """A peptide / DNA / RNA sequence drawn as a row of monomers (biopolymer.py has the chemistry)."""
    x: float                                                # top-left corner
    y: float
    polymer: str = 'PEPTIDE'                                # 'PEPTIDE' | 'DNA' | 'RNA'
    residues: list = field(default_factory=list)            # one-letter codes
    per_line: int = 10                                      # residues per line (wraps onto the next)
    three_letter: bool = False                              # peptides: Ala Gly ... instead of A G ...
    show_ids: bool = False                                  # residue numbers above the monomers
    cyclic: bool = False                                    # peptides: head-to-tail
    bridges: list = field(default_factory=list)             # peptides: [[i, j]] 0-based (Cys-Cys, Lys-Asp / Glu)
    n_term: str = 'H'                                       # peptides: 'H' | 'Ac'
    c_term: str = 'OH'                                      # peptides: 'OH' | 'NH2'
    p5: bool = False                                        # DNA / RNA: 5'-phosphate
    p3: bool = False                                        # DNA / RNA: 3'-phosphate
    color: str = ''


@dataclass
class Grid:
    """A stoichiometry grid under a reaction: one column per reactant / product. Each column is a dict with the
    compound's 'atoms' (ids) and the values typed into it (stoich.py says which; units g, mL, M, g/mL, mol, %)."""
    x: float                                                # top-left corner
    y: float
    reactants: list = field(default_factory=list)
    products: list = field(default_factory=list)
    hidden: list = field(default_factory=list)              # row keys not shown
    color: str = ''


@dataclass
class Fill:
    """A ring filled with a colour: the ring's atoms in order around it (it follows them when they move)."""
    atoms: list
    color: str = '#FFF59D'


@dataclass
class Group:
    """Objects that select, move and analyse as one (Object > Group): structures by their atoms (their bonds come with
    them) and any arrows, captions, shapes, grids, sequences and plates, as [[kind, id]]."""
    members: list = field(default_factory=list)


KINDS = {'atom': Atom, 'bond': Bond, 'arrow': Arrow, 'text': Text, 'shape': Shape, 'grid': Grid, 'fill': Fill,
         'seq': Sequence, 'tlc': TLC, 'group': Group}             # groups last: they point at the others
PLURAL = {'atom': 'atoms', 'bond': 'bonds', 'arrow': 'arrows', 'text': 'texts', 'shape': 'shapes', 'grid': 'grids',
          'fill': 'fills', 'seq': 'seqs', 'tlc': 'tlcs', 'group': 'groups'}
GROUPABLE = ('atom', 'arrow', 'text', 'shape', 'grid', 'seq', 'tlc')


def angle_of(dx, dy):
    return math.atan2(dy, dx)


def snap_angle(a, step_deg=15):
    s = math.radians(step_deg)
    return round(a / s) * s


class Doc:
    def __init__(self):
        self.atoms = {}
        self.bonds = {}
        self.arrows = {}
        self.texts = {}
        self.shapes = {}
        self.grids = {}
        self.fills = {}
        self.seqs = {}
        self.tlcs = {}
        self.groups = {}
        self.style = {}                                     # overrides of render.DEFAULT_STYLE
        self._next = 1

    # ------------------------------------------------------------------------------------------- items
    def _id(self):
        self._next += 1
        return self._next - 1

    def table(self, kind):
        return {'atom': self.atoms, 'bond': self.bonds, 'arrow': self.arrows, 'text': self.texts,
                'shape': self.shapes, 'grid': self.grids, 'fill': self.fills, 'seq': self.seqs,
                'tlc': self.tlcs, 'group': self.groups}[kind]

    def add_atom(self, x, y, el='C', charge=0):
        i = self._id()
        self.atoms[i] = Atom(x, y, el, charge)
        return i

    def atom_at_or_new(self, x, y, el='C', exclude=()):
        """The existing atom within MERGE of (x, y), else a new one there."""
        hit = self.find_atom(x, y, MERGE, exclude)
        return hit if hit is not None else self.add_atom(x, y, el)

    def add_bond(self, a, b, order=1, stereo=''):
        """Bond a-b (the existing one if the two are already bonded; None if a == b)."""
        if a == b:
            return None
        old = self.bond_between(a, b)
        if old is not None:
            return old
        i = self._id()
        self.bonds[i] = Bond(a, b, order, stereo)
        return i

    def add_arrow(self, x1, y1, x2, y2, kind='forward', bend=0.0):
        i = self._id()
        self.arrows[i] = Arrow(x1, y1, x2, y2, kind, bend)
        return i

    def add_text(self, x, y, text, **props):
        i = self._id()
        self.texts[i] = Text(x, y, text, **props)
        return i

    def add_shape(self, x1, y1, x2, y2, kind='rect'):
        i = self._id()
        self.shapes[i] = Shape(x1, y1, x2, y2, kind)
        return i

    def add_sequence(self, x, y, polymer, residues, **props):
        i = self._id()
        self.seqs[i] = Sequence(x, y, polymer, list(residues), **props)
        return i

    def bond_between(self, a, b):
        for i, bd in self.bonds.items():
            if (bd.a == a and bd.b == b) or (bd.a == b and bd.b == a):
                return i
        return None

    def neighbors(self, a):
        """[(neighbour atom, bond id)]"""
        out = []
        for i, bd in self.bonds.items():
            if bd.a == a:
                out.append((bd.b, i))
            elif bd.b == a:
                out.append((bd.a, i))
        return out

    def remove(self, kind, i):
        if kind == 'atom':
            self.remove_atom(i)
        elif kind == 'bond':
            self.remove_bond(i)
        else:
            self.table(kind).pop(i, None)

    def remove_atom(self, a):
        for _n, b in self.neighbors(a):
            self.remove_bond(b, prune=False)
        self.atoms.pop(a, None)
        for i in [i for i, f in self.fills.items() if a in f.atoms]:   # a ring fill goes with its ring
            del self.fills[i]

    def strip_atom(self, a):
        """Deleting a labelled atom (N, Cl, OMe, a charged carbon ...) turns it back into a plain carbon, keeping its
        bonds, as ChemDraw does. Returns False for an atom that already is a plain carbon (then it goes)."""
        at = self.atoms[a]
        if at.el == 'C' and not at.label and not at.charge:
            return False
        at.el, at.label, at.charge = 'C', '', 0
        return True

    def remove_bond(self, i, prune=True):
        """Remove a bond; a carbon left without bonds goes with it (it would show up as a lone CH4)."""
        bd = self.bonds.pop(i, None)
        if bd and prune:
            for a in (bd.a, bd.b):
                at = self.atoms.get(a)
                if at and at.el == 'C' and not at.charge and not self.neighbors(a):
                    del self.atoms[a]

    def merge_atom(self, src, dst):
        """Atom src is dropped onto dst: its bonds move to dst (a bond to dst itself is dropped) and src goes."""
        for n, b in self.neighbors(src):
            bd = self.bonds[b]
            if n == dst or self.bond_between(n, dst) is not None:
                del self.bonds[b]
                continue
            if bd.a == src:
                bd.a = dst
            else:
                bd.b = dst
        del self.atoms[src]

    def is_empty(self):
        return not (self.atoms or self.arrows or self.texts or self.shapes or self.grids or self.seqs or self.tlcs)

    # ------------------------------------------------------------------------------------------- hit testing
    def find_atom(self, x, y, r, exclude=()):
        best, bd = None, r
        for i, a in self.atoms.items():
            if i in exclude:
                continue
            d = math.hypot(a.x - x, a.y - y)
            if d <= bd:
                best, bd = i, d
        return best

    def find_bond(self, x, y, r):
        best, bd = None, r
        for i, b in self.bonds.items():
            p, q = self.atoms[b.a], self.atoms[b.b]
            d = seg_dist(x, y, p.x, p.y, q.x, q.y)
            if d <= bd:
                best, bd = i, d
        return best

    def find_arrow(self, x, y, r):
        best, bd = None, r
        for i, a in self.arrows.items():
            pts = arrow_points(a)
            d = min(seg_dist(x, y, *pts[k], *pts[k + 1]) for k in range(len(pts) - 1))
            if d <= bd:
                best, bd = i, d
        return best

    def find_shape(self, x, y, r):
        """A shape whose outline (not its inside: it may frame a structure) passes within r of (x, y)."""
        best, bd = None, r
        for i, s in self.shapes.items():
            x1, x2 = sorted((s.x1, s.x2)); y1, y2 = sorted((s.y1, s.y2))
            if s.kind in BIO_KINDS:                         # bio objects: the whole object counts, not its outline
                import biodraw
                if biodraw.contains(s, x, y, r):
                    return i
                continue
            if s.kind in ELEMENT_KINDS or s.kind == 'line':
                import elements
                if elements.contains(s, x, y, r, BL=float(self.style.get('bond_length', 14.4))):
                    return i
                continue
            if (s.style in ('filled', 'shaded', 'shadow') or s.kind == 'image') and \
                    x1 - r <= x <= x2 + r and y1 - r <= y <= y2 + r:       # a picture: anywhere on it
                return i
            if s.kind == 'line':
                d = seg_dist(x, y, s.x1, s.y1, s.x2, s.y2)
            elif s.kind == 'oval':
                rx, ry = max((x2 - x1) / 2, 1e-6), max((y2 - y1) / 2, 1e-6)
                cx, cy = (x1 + x2) / 2, (y1 + y2) / 2
                d = abs(math.hypot((x - cx) / rx, (y - cy) / ry) - 1) * min(rx, ry)
            elif s.kind in ('brackets', 'parens', 'braces'):
                d = min(seg_dist(x, y, x1, y1, x1, y2), seg_dist(x, y, x2, y1, x2, y2))
            else:
                d = min(seg_dist(x, y, x1, y1, x2, y1), seg_dist(x, y, x2, y1, x2, y2), seg_dist(x, y, x2, y2, x1, y2),
                        seg_dist(x, y, x1, y2, x1, y1))
            if d <= bd:
                best, bd = i, d
        return best

    # ------------------------------------------------------------------------------------------- chemistry
    def bond_sum(self, a):
        return sum(self.bonds[b].order for _n, b in self.neighbors(a))

    def _valences(self, a):
        at = self.atoms[a]
        base = VALENCE.get(at.el)
        if base is None:
            return None
        r = at.radical
        if not at.charge:
            return tuple(v - r for v in base if v - r >= 0) if r else base
        if at.el == 'H':                                    # H+ and H- make no bond
            return (0,)
        e = VALENCE_ELECTRONS.get(at.el)
        if e is None:
            return tuple(v - abs(at.charge) for v in base if v - abs(at.charge) >= 0)
        e -= at.charge
        return (max(0, min(e, 8 - e) - r),)

    def implicit_h(self, a):
        """Hydrogens the atom carries without them being drawn (0 for elements without a usual valence and for
        nicknames, whose hydrogens are in the group)."""
        if self.atoms[a].label:
            return 0
        if self.atoms[a].hcount is not None:                # typed with the label: CH2, NH ...
            return self.atoms[a].hcount
        vals = self._valences(a)
        if not vals:
            return 0
        s = self.bond_sum(a)
        for v in vals:
            if v >= s:
                return v - s
        return 0

    def valence_error(self, a):
        return bool(self.valence_problem(a))

    def valence_problem(self, a):
        """'' when the atom's bonds are fine, else what is wrong (for the red box, its tooltip and Check structure):
        more bonds than the element (with its charge) takes, or a nickname bonded more than once."""
        at = self.atoms[a]
        bonds = self.bond_sum(a)
        if at.label:
            if bonds > 1 and at.label in NICKNAMES:
                return f'{at.label} attaches by one bond, but has {bonds}.'
            return ''
        vals = self._valences(a)
        if at.hcount is not None:                           # the user said how many H: only too many is wrong
            if vals and at.el in ('C', 'N', 'O', 'F', 'B') and bonds + at.hcount > 4:
                return f'{at.el} with {at.hcount} H has room for {max(0, 4 - at.hcount)} bonds, but has {bonds}.'
            return ''
        if vals and bonds > max(vals):
            charge = f' {"+" if at.charge > 0 else "−"}{abs(at.charge)}' if at.charge else ''
            return (f'{at.el}{charge} takes at most {max(vals)} bond{"s" if max(vals) != 1 else ""}, but has {bonds}.')
        return ''

    def repeat_units(self, atom_ids=None):
        """Polymer brackets: [(repeat label, atom ids inside the pair)] for each bracket / parenthesis pair given a
        repeat count (Shape.params['repeat']: 'n', '3' ...) with atoms of atom_ids (None: any) inside it."""
        out = []
        for sh in self.shapes.values():
            rep = str(sh.params.get('repeat', '')).strip() if sh.kind in BRACKET_KINDS else ''
            if not rep:
                continue
            x1, x2 = sorted((sh.x1, sh.x2))
            y1, y2 = sorted((sh.y1, sh.y2))
            inside = {a for a, at in self.atoms.items() if (atom_ids is None or a in atom_ids)
                      and x1 < at.x < x2 and y1 < at.y < y2}
            if inside:
                out.append((rep, inside))
        return out

    def fragments(self, atom_ids=None):
        """Connected structures, as sorted lists of atom ids (only atoms in atom_ids, when given)."""
        todo = set(self.atoms if atom_ids is None else (a for a in atom_ids if a in self.atoms))
        out = []
        while todo:
            seed = min(todo)
            comp, stack = {seed}, [seed]
            while stack:
                for n, _b in self.neighbors(stack.pop()):
                    if n in todo and n not in comp:
                        comp.add(n); stack.append(n)
            todo -= comp
            out.append(sorted(comp))
        out.sort(key=lambda c: (self.atoms[c[0]].x, self.atoms[c[0]].y))
        return out

    # ------------------------------------------------------------------------------------------- geometry
    def bond_angles(self, a, skip=None):
        at = self.atoms[a]
        return [angle_of(self.atoms[n].x - at.x, self.atoms[n].y - at.y) for n, b in self.neighbors(a) if b != skip]

    def next_bond_angle(self, a):
        """Direction of a new bond from atom a: up-right for a lone atom, zigzag (away from the crowd) after one bond,
        straight on after a triple bond or two double bonds, else the middle of the widest gap between the bonds."""
        nb = self.neighbors(a)
        at = self.atoms[a]
        if not nb:
            return math.radians(-30)
        angles = self.bond_angles(a)
        if len(nb) == 1:
            base = angles[0]
            if any(self.bonds[b].order == 3 for _n, b in nb):
                return base + math.pi
            cands = [base + math.radians(120), base - math.radians(120)]

            def room(ang):
                x, y = at.x + math.cos(ang), at.y + math.sin(ang)
                others = [math.hypot(o.x - x, o.y - y) for i, o in self.atoms.items() if i != a and i != nb[0][0]]
                return min(others, default=9.0)
            r = [room(c) for c in cands]
            if abs(r[0] - r[1]) > 1e-6:
                return cands[0] if r[0] > r[1] else cands[1]
            return min(cands, key=lambda c: abs(math.sin(c)))   # a free choice: the more horizontal one (a zigzag)
        angles.sort()
        gaps = [((angles[(k + 1) % len(angles)] - angles[k]) % (2 * math.pi) or 2 * math.pi, angles[k])
                for k in range(len(angles))]
        gap, start = max(gaps)
        return start + gap / 2

    def free_direction(self, a):
        """Unit vector pointing away from the bonds of atom a (where its hydrogens and labels can go)."""
        at = self.atoms[a]
        sx = sy = 0.0
        for n, _b in self.neighbors(a):
            o = self.atoms[n]
            d = math.hypot(o.x - at.x, o.y - at.y) or 1.0
            sx += (o.x - at.x) / d; sy += (o.y - at.y) / d
        d = math.hypot(sx, sy)
        if d < 1e-6:
            return (1.0, 0.0) if not self.neighbors(a) else (0.0, -1.0)
        return (-sx / d, -sy / d)

    def ring_of_bond(self, i, limit=8):
        """Atoms of the smallest ring through bond i (in order around the ring), or None."""
        bd = self.bonds[i]
        prev = {bd.a: None}
        frontier = [bd.a]
        for _depth in range(limit - 1):
            nxt = []
            for u in frontier:
                for n, b in self.neighbors(u):
                    if b == i or n in prev:
                        continue
                    prev[n] = u
                    if n == bd.b:
                        path = [n]
                        while prev[path[-1]] is not None:
                            path.append(prev[path[-1]])
                        return path
                    nxt.append(n)
            frontier = nxt
        return None

    def double_side(self, i):
        """Where a double bond's second line goes: +1 / -1 = that side of a->b (left normal = +1), 0 = centred."""
        bd = self.bonds[i]
        if bd.side is not None:
            return bd.side
        p, q = self.atoms[bd.a], self.atoms[bd.b]
        nx, ny = -(q.y - p.y), q.x - p.x                   # normal of a->b

        def side_of(x, y):
            v = (x - p.x) * nx + (y - p.y) * ny
            return 0 if abs(v) < 1e-6 else (1 if v > 0 else -1)
        ring = self.ring_of_bond(i)
        if ring:
            cx = sum(self.atoms[k].x for k in ring) / len(ring)
            cy = sum(self.atoms[k].y for k in ring) / len(ring)
            s = side_of(cx, cy)
            if s:
                return s
        na = [n for n, _b in self.neighbors(bd.a) if n != bd.b]
        nb = [n for n, _b in self.neighbors(bd.b) if n != bd.a]
        if not na and not nb:
            return 0
        if (not na and len(nb) != 1) or (not nb and len(na) != 1):
            return 0                                        # C=O of a carbonyl, =CH2 on a branch point: centred
        total = sum(side_of(self.atoms[n].x, self.atoms[n].y) for n in na + nb)
        if total == 0:
            return 0
        return 1 if total > 0 else -1

    def ring_points(self, n, center, start_angle):
        R = 1 / (2 * math.sin(math.pi / n))
        cx, cy = center
        return [(cx + R * math.cos(start_angle + 2 * math.pi * k / n), cy + R * math.sin(start_angle + 2 * math.pi * k / n))
                for k in range(n)]

    def ring_on_bond(self, i, n):
        """Points of an n-ring fused on bond i, on the side away from the bond's other neighbours."""
        bd = self.bonds[i]
        p, q = self.atoms[bd.a], self.atoms[bd.b]
        mx, my = (p.x + q.x) / 2, (p.y + q.y) / 2
        L = math.hypot(q.x - p.x, q.y - p.y) or 1.0
        nx, ny = -(q.y - p.y) / L, (q.x - p.x) / L
        score = 0.0
        for a in (bd.a, bd.b):
            for o, _b in self.neighbors(a):
                if o in (bd.a, bd.b):
                    continue
                score += (self.atoms[o].x - mx) * nx + (self.atoms[o].y - my) * ny
        s = -1 if score > 0 else 1
        apothem = 1 / (2 * math.tan(math.pi / n))
        cx, cy = mx + s * nx * apothem, my + s * ny * apothem
        return self.ring_points(n, (cx, cy), angle_of(p.x - cx, p.y - cy))

    def ring_on_atom(self, a, n):
        at = self.atoms[a]
        ang = self.next_bond_angle(a) if self.neighbors(a) else math.radians(90)
        R = 1 / (2 * math.sin(math.pi / n))
        cx, cy = at.x + R * math.cos(ang), at.y + R * math.sin(ang)
        return self.ring_points(n, (cx, cy), ang + math.pi)

    def add_ring(self, pts, aromatic=False):
        """Atoms and bonds of a ring through pts (existing atoms are reused); aromatic: alternate double bonds where
        they fit with the double bonds already there. Returns the ring's atom ids."""
        ids = [self.atom_at_or_new(x, y) for x, y in pts]
        n = len(ids)
        edges = [(ids[k], ids[(k + 1) % n]) for k in range(n)]
        bonds = [self.add_bond(u, v) for u, v in edges]
        if aromatic:                                        # an odd ring (cyclopentadiene) gets (n - 1) / 2 doubles
            def doubles(par):
                return list(range(par, n - n % 2, 2))

            def has_double(atom, skip):
                return any(self.bonds[b].order >= 2 for _o, b in self.neighbors(atom) if b not in skip)

            def conflicts(par):
                c = 0
                for k in doubles(par):
                    u, v = edges[k]
                    if self.bonds[bonds[k]].order == 1 and (has_double(u, bonds) or has_double(v, bonds)):
                        c += 1
                for k in set(range(n)) - set(doubles(par)):
                    if self.bonds[bonds[k]].order >= 2:
                        c += 1
                return c
            par = min((0, 1), key=conflicts)
            for k in doubles(par):
                u, v = edges[k]
                b = self.bonds[bonds[k]]
                if b.order == 1 and not has_double(u, bonds) and not has_double(v, bonds):
                    b.order = 2
        return ids

    def put_group(self, a, smiles):
        """Atom a becomes the first atom of the group (a SMILES); the rest is laid out away from a's bonds."""
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        at = self.atoms[a]
        first = mol.GetAtomWithIdx(0)
        at.el, at.charge, at.label = first.GetSymbol(), first.GetFormalCharge(), ''
        if mol.GetNumAtoms() == 1:
            return [a]
        rdDepictor.Compute2DCoords(mol)
        conf = mol.GetConformer()
        pos = [conf.GetAtomPosition(k) for k in range(mol.GetNumAtoms())]
        pts = [(p.x - pos[0].x, -(p.y - pos[0].y)) for p in pos]   # y down, atom 0 at the origin
        L = sorted(math.hypot(pts[b.GetBeginAtomIdx()][0] - pts[b.GetEndAtomIdx()][0],
                              pts[b.GetBeginAtomIdx()][1] - pts[b.GetEndAtomIdx()][1]) for b in mol.GetBonds())
        L = L[len(L) // 2] or 1.0
        cx = sum(x for x, _y in pts[1:]) / (len(pts) - 1); cy = sum(y for _x, y in pts[1:]) / (len(pts) - 1)
        fx, fy = self.free_direction(a) if self.neighbors(a) else (1.0, 0.0)
        rot = angle_of(fx, fy) - angle_of(cx, cy)
        c, s = math.cos(rot), math.sin(rot)
        ids = [a]
        for k in range(1, len(pts)):
            x, y = pts[k][0] / L, pts[k][1] / L
            atom = mol.GetAtomWithIdx(k)
            ids.append(self.add_atom(at.x + x * c - y * s, at.y + x * s + y * c, atom.GetSymbol(), atom.GetFormalCharge()))
        for b in mol.GetBonds():
            order = {Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}.get(b.GetBondType(), 1)
            self.add_bond(ids[b.GetBeginAtomIdx()], ids[b.GetEndAtomIdx()], order)
        return ids

    def set_nickname(self, a, label):
        """Atom a shows label (a key of NICKNAMES); its element becomes the group's attaching atom."""
        from rdkit import Chem
        first = Chem.MolFromSmiles(NICKNAMES[label]).GetAtomWithIdx(0)
        at = self.atoms[a]
        at.label, at.el, at.charge = label, first.GetSymbol(), 0

    def sprout(self, a, smiles='C', order=1, stereo=''):
        """A new bond from atom a (in the next free direction) to the group smiles; returns the new atom."""
        at = self.atoms[a]
        ang = self.next_bond_angle(a)
        b = self.atom_at_or_new(at.x + math.cos(ang), at.y + math.sin(ang), exclude={a})
        bid = self.add_bond(a, b, order, stereo)
        if bid is not None and stereo and self.bonds[bid].a != a:
            bd = self.bonds[bid]; bd.a, bd.b = bd.b, bd.a
        if smiles != 'C' and not self.neighbors(b)[1:]:
            self.put_group(b, smiles)
        return b

    def rings(self, atom_ids=None):
        """The smallest set of smallest rings (RDKit), as atom-id lists in order around each ring; with atom_ids, only
        the rings all of whose atoms are in it."""
        from rdkit import Chem
        block, order = self.to_molblock()
        if not order:
            return []
        mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
        if mol is None:
            return []
        mol.UpdatePropertyCache(strict=False)
        out = []
        for ring in Chem.GetSymmSSSR(mol):                  # the smallest rings (naphthalene: two six-rings)
            ids = [order[k] for k in ring]
            if atom_ids is None or set(ids) <= set(atom_ids):
                out.append(ids)
        return out

    def fill_of_ring(self, ring):
        """The id of the fill on exactly these atoms, or None."""
        key = set(ring)
        for i, f in self.fills.items():
            if set(f.atoms) == key:
                return i
        return None

    def chain_points(self, x0, y0, x1, y1):
        """Zigzag chain from (x0, y0) toward (x1, y1) (direction snapped to 30 degrees): the points after the first."""
        ang = snap_angle(angle_of(x1 - x0, y1 - y0), 30)
        dist = math.hypot(x1 - x0, y1 - y0)
        step = math.sqrt(3) / 2
        n = max(1, round(dist / step))
        ux, uy = math.cos(ang), math.sin(ang)
        vx, vy = -uy, ux
        sgn = -1 if vy > 0 else 1                           # the first zigzag step goes up on the screen
        if abs(uy) > 0.99:                                  # (to the right for a vertical chain)
            sgn = 1 if vx > 0 else -1
        return [(x0 + k * step * ux + (k % 2) * 0.5 * vx * sgn, y0 + k * step * uy + (k % 2) * 0.5 * vy * sgn)
                for k in range(1, n + 1)]

    # ------------------------------------------------------------------------------------------- transforms
    def atoms_of(self, sel):
        """Atom ids moved along with a selection [(kind, id)] (the atoms of selected bonds too)."""
        out = set()
        for kind, i in sel:
            if kind == 'atom':
                out.add(i)
            elif kind == 'bond' and i in self.bonds:
                out.update((self.bonds[i].a, self.bonds[i].b))
        return out

    def points_of(self, sel):
        pts = [(self.atoms[a].x, self.atoms[a].y) for a in self.atoms_of(sel) if a in self.atoms]
        for kind, i in sel:
            if kind == 'arrow' and i in self.arrows:
                ar = self.arrows[i]; pts += [(ar.x1, ar.y1), (ar.x2, ar.y2)]
            elif kind == 'shape' and i in self.shapes:
                sh = self.shapes[i]; pts += [(sh.x1, sh.y1), (sh.x2, sh.y2)]
            elif kind == 'grid' and i in self.grids:
                pts.append((self.grids[i].x, self.grids[i].y))
            elif kind == 'text' and i in self.texts:
                pts.append((self.texts[i].x, self.texts[i].y))
            elif kind == 'seq' and i in self.seqs:
                pts.append((self.seqs[i].x, self.seqs[i].y))
            elif kind == 'tlc' and i in self.tlcs:
                t = self.tlcs[i]; pts += [(t.x, t.y), (t.x + t.w, t.y + t.h)]
        return pts

    def anchor(self, kind, i):
        """Where an arrow attached to atom / bond i sits: the atom, or the bond's middle; None when it is gone."""
        if kind == 'atom' and i in self.atoms:
            return self.atoms[i].x, self.atoms[i].y
        if kind == 'bond' and i in self.bonds:
            a, b = self.atoms[self.bonds[i].a], self.atoms[self.bonds[i].b]
            return (a.x + b.x) / 2, (a.y + b.y) / 2
        return None

    def attach_arrow_end(self, ar, end, item):
        """Ties end ('tail' / 'head') of arrow ar to item (atom / bond), keeping where it is now as the offset."""
        x, y = (ar.x1, ar.y1) if end == 'tail' else (ar.x2, ar.y2)
        p = self.anchor(*item) if item else None
        setattr(ar, end + '_on', [item[0], item[1], round(x - p[0], 4), round(y - p[1], 4)] if p else None)

    def update_attached_arrows(self):
        """Moves the attached ends of arrows with their atoms and bonds; an end whose atom / bond is gone comes loose."""
        for ar in self.arrows.values():
            for end in ('tail', 'head'):
                on = getattr(ar, end + '_on')
                if not on:
                    continue
                p = self.anchor(on[0], on[1])
                if p is None:
                    setattr(ar, end + '_on', None)
                    continue
                x, y = p[0] + on[2], p[1] + on[3]
                if end == 'tail':
                    ar.x1, ar.y1 = x, y
                else:
                    ar.x2, ar.y2 = x, y

    def lone_pair_count(self, a):
        """Lone pairs an atom has: (valence electrons - charge - bonds - H - unpaired) / 2 (O in water: 2)."""
        at = self.atoms[a]
        e = VALENCE_ELECTRONS.get(at.el)
        if e is None or at.label:
            return 0
        free = e - at.charge - self.bond_sum(a) - self.implicit_h(a) - at.radical
        return max(0, free // 2)

    def transform(self, sel, fn):
        """Apply fn(x, y) -> (x, y) to everything in sel."""
        for a in self.atoms_of(sel):
            at = self.atoms[a]; at.x, at.y = fn(at.x, at.y)
        for kind, i in sel:
            if kind == 'arrow' and i in self.arrows:
                ar = self.arrows[i]
                ar.x1, ar.y1 = fn(ar.x1, ar.y1); ar.x2, ar.y2 = fn(ar.x2, ar.y2)
            elif kind == 'shape' and i in self.shapes:
                sh = self.shapes[i]
                sh.x1, sh.y1 = fn(sh.x1, sh.y1); sh.x2, sh.y2 = fn(sh.x2, sh.y2)
            elif kind == 'grid' and i in self.grids:
                g = self.grids[i]; g.x, g.y = fn(g.x, g.y)
            elif kind == 'text' and i in self.texts:
                t = self.texts[i]; t.x, t.y = fn(t.x, t.y)
            elif kind == 'seq' and i in self.seqs:
                q = self.seqs[i]; q.x, q.y = fn(q.x, q.y)
            elif kind == 'tlc' and i in self.tlcs:             # moved / scaled by its corners, kept upright
                t = self.tlcs[i]
                (ax, ay), (bx, by) = fn(t.x, t.y), fn(t.x + t.w, t.y + t.h)
                t.x, t.y, t.w, t.h = min(ax, bx), min(ay, by), abs(bx - ax), abs(by - ay)

    # ------------------------------------------------------------------------------------------- groups
    def group_of(self, item):
        """The id of the group holding item (a bond: through its atoms), or None."""
        kind, i = item
        if kind == 'bond':
            if i not in self.bonds:
                return None
            kind, i = 'atom', self.bonds[i].a
        for g, grp in self.groups.items():
            if [kind, i] in grp.members:
                return g
        return None

    def group_items(self, g):
        """Everything group g holds, as [(kind, id)]: its members and the bonds between its atoms."""
        items = [(k, i) for k, i in self.groups[g].members if i in self.table(k)]
        atoms = {i for k, i in items if k == 'atom'}
        return items + [('bond', b) for b, bd in self.bonds.items() if bd.a in atoms and bd.b in atoms]

    def with_groups(self, items):
        """items, plus everything in the groups they belong to."""
        out = set(items)
        for g in {self.group_of(it) for it in items} - {None}:
            out.update(self.group_items(g))
        return out

    def make_group(self, items):
        """Groups the objects in items (whole structures: every atom of a structure touched), taking in any group
        they already belong to; returns the group's id, or None for fewer than two objects."""
        items = self.with_groups(items)
        touched = self.atoms_of(items)
        frags = [set(f) for f in self.fragments() if touched & set(f)]
        atoms = set().union(*frags)
        others = sorted((k, i) for k, i in items if k in GROUPABLE and k != 'atom' and i in self.table(k))
        if len(frags) + len(others) < 2:
            return None
        for g in {self.group_of(it) for it in items} - {None}:
            del self.groups[g]
        g = self._id()
        self.groups[g] = Group([['atom', a] for a in sorted(atoms)] + [[k, i] for k, i in others])
        return g

    def ungroup(self, items):
        """Dissolves the groups that items belong to; returns how many."""
        found = {self.group_of(it) for it in items} - {None}
        for g in found:
            del self.groups[g]
        return len(found)

    def join_group(self, item, near):
        """Puts item into the group of the structure holding the atoms near (a new group of the two when it has
        none), as a caption pasted under a structure stays with it."""
        near = [a for a in near if a in self.atoms]
        if not near:
            return
        g = self.group_of(('atom', near[0]))
        if g is None:
            atoms = set().union(*[set(f) for f in self.fragments() if set(near) & set(f)])
            g = self._id()
            self.groups[g] = Group([['atom', a] for a in sorted(atoms)])
        if list(item) not in self.groups[g].members:
            self.groups[g].members.append(list(item))

    def prune_groups(self):
        """Forgets deleted members (atoms bonded to a grouped structure join it); drops groups left with fewer than
        two objects."""
        if not self.groups:
            return
        frags = [set(f) for f in self.fragments()]
        for g in list(self.groups):
            grp = self.groups[g]
            members = [[k, i] for k, i in grp.members if i in self.table(k)]
            atoms = {i for k, i in members if k == 'atom'}
            mine = [f for f in frags if f & atoms]
            others = [m for m in members if m[0] != 'atom']
            if len(mine) + len(others) < 2:
                del self.groups[g]
            else:
                grp.members = [['atom', a] for a in sorted(set().union(*mine))] + others

    def all_items(self):
        return ([('atom', i) for i in self.atoms] + [('bond', i) for i in self.bonds] +
                [('arrow', i) for i in self.arrows] + [('text', i) for i in self.texts] +
                [('shape', i) for i in self.shapes] + [('grid', i) for i in self.grids] +
                [('seq', i) for i in self.seqs] + [('tlc', i) for i in self.tlcs])

    # ------------------------------------------------------------------------------------------- files
    def to_dict(self, sel=None):
        """The drawing (or the items in sel, with the bonds between selected atoms) as plain data."""
        if sel is None:
            keep = None
        else:
            keep = set(sel)
            atoms = self.atoms_of(sel)
            keep |= {('atom', a) for a in atoms}
            keep |= {('bond', i) for i, b in self.bonds.items() if b.a in atoms and b.b in atoms}
            keep |= {('fill', i) for i, f in self.fills.items() if set(f.atoms) <= atoms}
            keep |= {('group', i) for i, g in self.groups.items() if all(tuple(m) in keep for m in g.members)}

        def items(kind):
            return [dict(id=i, **asdict(v)) for i, v in self.table(kind).items() if keep is None or (kind, i) in keep]
        return dict(format='structure-editor', version=2, style=dict(self.style),
                    **{PLURAL[k]: items(k) for k in KINDS})

    @classmethod
    def from_dict(cls, d):
        doc = cls()
        doc.insert_dict(d, keep_ids=True)
        doc.style = dict(d.get('style') or {})
        return doc

    def insert_dict(self, d, dx=0.0, dy=0.0, keep_ids=False):
        """Add the items of a to_dict() result (shifted by dx, dy); returns the new items as [(kind, id)]."""
        new, idmap, moved = [], {}, {}
        for kind in KINDS:
            for item in d.get(PLURAL[kind], []):
                item = dict(item)
                old = item.pop('id')
                if kind == 'bond':
                    if item['a'] not in idmap or item['b'] not in idmap:
                        continue
                    item['a'], item['b'] = idmap[item['a']], idmap[item['b']]
                if kind == 'fill':                          # a ring fill follows its atoms, or goes with them
                    if any(a not in idmap for a in item['atoms']):
                        continue
                    item['atoms'] = [idmap[a] for a in item['atoms']]
                if kind == 'group':                         # members follow their new ids; the lost are dropped
                    item['members'] = [[k, moved[(k, j)]] for k, j in item.get('members', []) if (k, j) in moved]
                    if not item['members']:
                        continue
                if kind == 'text' and item.get('on_arrow'):  # conditions follow their arrow's new id
                    on = item['on_arrow']
                    item['on_arrow'] = [moved[('arrow', on[0])], on[1]] if ('arrow', on[0]) in moved else None
                if kind == 'arrow':                         # attached ends follow their atom / bond, or come loose
                    for end in ('tail_on', 'head_on'):
                        on = item.get(end)
                        if on:
                            item[end] = [on[0], moved[(on[0], on[1])], *on[2:]] if (on[0], on[1]) in moved else None
                if kind == 'grid':                          # its columns point at atoms: follow them, drop the lost
                    for side in ('reactants', 'products'):
                        cols = []
                        for col in item.get(side, []):
                            col = dict(col, atoms=[idmap[a] for a in col.get('atoms', []) if a in idmap])
                            if col['atoms']:
                                cols.append(col)
                        item[side] = cols
                for k in ('x', 'x1', 'x2'):
                    if k in item:
                        item[k] += dx
                for k in ('y', 'y1', 'y2'):
                    if k in item:
                        item[k] += dy
                if keep_ids:
                    i = old
                    self._next = max(self._next, i + 1)
                else:
                    i = self._id()
                if kind == 'atom':
                    idmap[old] = i
                moved[(kind, old)] = i
                self.table(kind)[i] = KINDS[kind](**item)
                if kind != 'group':                         # a group is not an object to select
                    new.append((kind, i))
        return new

    def copy(self):
        return copy.deepcopy(self)

    def expanded(self):
        """A copy with every nickname written out as atoms (for the formula, SMILES, MOL): (doc, {atom id: the ids it
        became}). Atom ids of the drawing stay the same in the copy."""
        d = self.copy()
        grown = {}
        for a in [a for a, at in d.atoms.items() if at.label]:
            smiles = NICKNAMES.get(d.atoms[a].label)
            if smiles is None:
                d.atoms[a].label = ''
                continue
            before = set(d.atoms)
            d.put_group(a, smiles)
            grown[a] = set(d.atoms) - before
        return d, grown

    def to_molblock(self, atom_ids=None, title=''):
        """V2000 MOL block of the structures (or of the atoms in atom_ids) with wedges and charges; y up, 1.5 A bonds.
        Returns (text, [atom id in file order])."""
        order = [a for a in self.atoms if atom_ids is None or a in atom_ids]
        index = {a: k + 1 for k, a in enumerate(order)}
        bonds = [b for b in self.bonds.values() if b.a in index and b.b in index]
        lines = [title, '  Kekulex 2D', '',
                 f'{len(order):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000']
        for a in order:
            at = self.atoms[a]
            lines.append(f'{at.x * MOL_BOND:10.4f}{-at.y * MOL_BOND:10.4f}{0:10.4f} {at.el:<3} 0  0  0  0  0  0  0  0  0  0  0  0')
        for b in bonds:
            st = {'wedge': 1, 'hash': 6}.get(b.stereo, 0)
            lines.append(f'{index[b.a]:3d}{index[b.b]:3d}{b.order:3d}{st:3d}')
        charged = [(index[a], self.atoms[a].charge) for a in order if self.atoms[a].charge]
        for k in range(0, len(charged), 8):
            chunk = charged[k:k + 8]
            lines.append(f'M  CHG{len(chunk):3d}' + ''.join(f'{i:4d}{c:4d}' for i, c in chunk))
        radicals = [(index[a], 2 if self.atoms[a].radical == 1 else 3) for a in order if self.atoms[a].radical]
        for k in range(0, len(radicals), 8):                # 2 doublet (•), 3 triplet (carbene, ••)
            chunk = radicals[k:k + 8]
            lines.append(f'M  RAD{len(chunk):3d}' + ''.join(f'{i:4d}{c:4d}' for i, c in chunk))
        lines.append('M  END')
        return '\n'.join(lines) + '\n', order

    def insert_rdkit(self, mol, center=(0.0, 0.0)):
        """Add an RDKit molecule (laid out in 2D when it has no coordinates), centred on center; returns the new items."""
        from rdkit import Chem
        mol = Chem.Mol(mol)
        if mol.GetNumConformers() == 0 or mol.GetConformer().Is3D():
            import layout2d
            mol = layout2d.with_fusion_hydrogens(mol)
            layout2d.compute(mol)
        try:
            Chem.ReapplyMolBlockWedging(mol)
        except (AttributeError, RuntimeError):
            pass
        if not any(b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) for b in mol.GetBonds()):
            Chem.WedgeMolBonds(mol, mol.GetConformer())
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            pass
        conf = mol.GetConformer()
        pos = [conf.GetAtomPosition(k) for k in range(mol.GetNumAtoms())]
        lens = sorted(math.hypot(pos[b.GetBeginAtomIdx()].x - pos[b.GetEndAtomIdx()].x,
                                 pos[b.GetBeginAtomIdx()].y - pos[b.GetEndAtomIdx()].y) for b in mol.GetBonds())
        L = lens[len(lens) // 2] if lens else MOL_BOND
        L = L or MOL_BOND
        cx = sum(p.x for p in pos) / max(len(pos), 1)
        cy = sum(p.y for p in pos) / max(len(pos), 1)
        ids, new = [], []
        for k, atom in enumerate(mol.GetAtoms()):
            i = self.add_atom(center[0] + (pos[k].x - cx) / L, center[1] - (pos[k].y - cy) / L, atom.GetSymbol(),
                              atom.GetFormalCharge())
            if atom.GetNumRadicalElectrons():
                self.atoms[i].radical = min(atom.GetNumRadicalElectrons(), 2)
            ids.append(i); new.append(('atom', i))
        for b in mol.GetBonds():
            bt = b.GetBondType()
            order = 3 if bt == Chem.BondType.TRIPLE else 2 if bt == Chem.BondType.DOUBLE else 1
            stereo = {Chem.BondDir.BEGINWEDGE: 'wedge', Chem.BondDir.BEGINDASH: 'hash'}.get(b.GetBondDir(), '')
            i = self._id()
            self.bonds[i] = Bond(ids[b.GetBeginAtomIdx()], ids[b.GetEndAtomIdx()], order, stereo)
            new.append(('bond', i))
        return new


def arrow_points(ar, n=16):
    """Points along an arrow's shaft (a straight arrow: its two ends; a curved one: n + 1 points on its arc)."""
    if ar.kind not in ('curved', 'curved_half') or not ar.bend:
        return [(ar.x1, ar.y1), (ar.x2, ar.y2)]
    cx, cy = curve_control(ar)
    return [((1 - t) ** 2 * ar.x1 + 2 * (1 - t) * t * cx + t * t * ar.x2,
             (1 - t) ** 2 * ar.y1 + 2 * (1 - t) * t * cy + t * t * ar.y2) for t in (k / n for k in range(n + 1))]


def curve_control(ar):
    """Control point of a curved arrow's quadratic arc: off the middle of tail-head by bend x length (left normal +)."""
    dx, dy = ar.x2 - ar.x1, ar.y2 - ar.y1
    return (ar.x1 + ar.x2) / 2 + dy * ar.bend * 2, (ar.y1 + ar.y2) / 2 - dx * ar.bend * 2


def seg_dist(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / L2))
    return math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))
