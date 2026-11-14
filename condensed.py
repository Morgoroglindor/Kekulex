"""
condensed.py - a typed atom label read as a condensed formula, as ChemDraw reads one: SO3K, CH2CH3, COCH3, N(CH3)2,
CH(CH3)2, PO3H2, OCF3, SO2NH2, CONHMe, CH2OMe, (CH2)3CH3, CH=CH2, CO2-, ... The first atom is the one on the
structure. Known nicknames (Me, Et, Ph, Boc, iPr ...) may appear inside.

    condensed.smiles('SO3K')            # 'S(=O)(=O)[O-].[K+]' (first atom: the S)

How the atoms join, left to right:
  - each atom bonds to the latest atom that still has a free bond (O3, H2 ... all go on the same one while it can
    take them), so CH2CH3 is a chain, OC2H5 an ethoxy and N(CH3)2 a dimethylamino;
  - H fills the latest atoms first (SO3H: the H on the last O);
  - a single O or S on a carbon that still has two free bonds after it is a C=O (COCH3, OCOCH3, CONH2);
  - at the end, bonds between atoms both left with free bonds become double, then triple (CO2H, SO3H, CN, NCO);
  - '=' and '#' (or ≡) set the next bond; a trailing + / - (2+, 3- ...) charges the last atom;
  - Li, Na, K, Rb, Cs on O, S or N become a salt (O- K+).
Returns None for text that is not such a formula, or for a single atom (OH, NH2: element labels).
"""
import re
from functools import lru_cache

from rdkit import Chem

# valence electrons: a neutral atom makes min(e, 8 - e) bonds, P and S may make 5 and 6 (as in model.py)
_ELECTRONS = {'H': 1, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Se': 6, 'Br': 7,
              'I': 7, 'As': 5, 'Te': 6, 'Ge': 4, 'Sn': 4, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1, 'Mg': 2, 'Zn': 2,
              'Al': 3, 'Ca': 2}
_HYPER = {'P': 5, 'S': 6, 'Se': 6, 'As': 5, 'Te': 6, 'I': 7, 'Cl': 7, 'Br': 5}
_SALT_METALS = ('Li', 'Na', 'K', 'Rb', 'Cs')
_ALIASES = {'C6H5': 'Ph', 'C6H11': 'Cy'}                    # ring groups written as formulas


def _valence(sym, charge=0):
    e = _ELECTRONS[sym] - charge
    if charge == 0 and sym in _HYPER:
        return _HYPER[sym]
    return min(e, 8 - e) if e > 0 else 0


def _tokens(text, nicknames):
    """[(kind, value, count)]: ('el', symbol), ('nick', label), ('open', None), ('close', None), ('bond', order),
    ('charge', n); counts from the digits after a symbol or ')'."""
    table = Chem.GetPeriodicTable()
    plain = lambda k: all(t in _ELECTRONS for t in re.findall(r'[A-Z][a-z]?', k)) and re.fullmatch(r'(?:[A-Z][a-z]?)+', k)
    names = sorted((k for k in nicknames if k.isalpha() and len(k) > 1 and not plain(k)), key=len, reverse=True)

    def element(s):
        return s in _ELECTRONS
    out, i = [], 0
    while i < len(text):
        ch = text[i]
        if ch.isspace():
            i += 1
            continue
        if ch in '([':
            out.append(['open', None, 1]); i += 1; continue
        if ch in ')]':
            out.append(['close', None, 1]); i += 1
        elif ch in '=#≡':
            out.append(['bond', 2 if ch == '=' else 3, 1]); i += 1
            continue
        elif ch in '+-':
            m = re.match(r'[+-]+', text[i:])
            signs = m.group(0)
            if out and out[-1][0] == 'num':                 # 2+ written as digits then sign
                n = out.pop()[1]
            else:
                n = len(signs)
            out.append(['charge', n if signs[0] == '+' else -n, 1]); i += len(signs)
            continue
        else:
            nick = next((n for n in names if text.startswith(n, i) and not text[i + len(n):i + len(n) + 1].islower()),
                        None)
            if nick:
                out.append(['nick', nick, 1]); i += len(nick)
            elif ch.isupper() and text[i + 1:i + 2].islower() and element(text[i:i + 2]):
                out.append(['el', text[i:i + 2], 1]); i += 2
            elif ch.isupper() and element(ch):
                out.append(['el', ch, 1]); i += 1
            elif ch.isdigit():
                m = re.match(r'\d+', text[i:])
                out.append(['num', int(m.group(0)), 1]); i += len(m.group(0))
                continue
            else:
                return None
        m = re.match(r'\d+', text[i:])                      # a count after a symbol or ')' (CO2-: two O, then -)
        if m:
            out[-1][2] = int(m.group(0)); i += len(m.group(0))
    if any(t[0] == 'num' for t in out):
        return None
    return out


class _Builder:
    def __init__(self):
        self.mol = Chem.RWMol()
        self.free = []                                      # free bonds left on each atom
        self.hs = []

    def add(self, sym, charge=0):
        atom = Chem.Atom(sym)
        atom.SetFormalCharge(charge)
        a = self.mol.AddAtom(atom)
        self.free.append(_valence(sym, charge) - (1 if a == 0 else 0))   # the first atom keeps a bond for the structure
        self.hs.append(0)
        return a

    def bond(self, a, b, order=1):
        self.mol.AddBond(a, b, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order])
        self.free[a] -= order
        self.free[b] -= order


def _parse(text, nicknames):
    toks = _tokens(text, nicknames)
    if not toks:
        return None
    if toks[-1][0] == 'charge':                             # NMe3+, SO3-: the charge sits on the last atom written
        charge = toks.pop()[1]
        last, depth = None, 0
        for t in reversed(toks):                            # N(CH3)3+: the N, outside the brackets
            depth += {'close': 1, 'open': -1}.get(t[0], 0)
            if depth == 0 and t[0] == 'el' and t[1] != 'H':
                last = t
                break
        if last is None:
            return None
        last.append(charge)
    b = _Builder()
    scopes = [[]]                                           # atoms of each open bracket, outermost first
    group_start = []                                        # where each bracket's tokens begin, to repeat them
    pending_order = 1
    last_heavy = None

    pendant_of = {}                                         # SO2, CO2, PO3: oxygens on their atom

    def latest_free(scope, need=1):
        for a in reversed(scope):
            if b.free[a] >= need:
                host = pendant_of.get(a)
                if host is not None:                        # SO2NH2, SO2Cl: the next goes on the S while it still
                    oxo = sum(b.free[p] for p, h in pendant_of.items() if h == host)   # has room after S=O, S=O
                    if b.free[host] - oxo >= need:
                        return host
                return a
        return None

    entries = [None]                                        # where each open bracket's first atom bonds

    def attach(new, scope, order=1, onto=None):
        target = onto if onto is not None else latest_free(scope, order)
        if target is None and not scope and entries[-1] is not None and b.free[entries[-1]] >= order:
            target = entries[-1]
        if target is None and len(scopes) > 1:              # a bracket's first atom bonds outside it
            for outer in reversed(scopes[:-1]):
                target = latest_free(outer, order)
                if target is not None:
                    break
        if target is None:
            if new != 0:
                raise ValueError('no free bond')
            return
        b.bond(target, new, order)

    def run(tokens):
        nonlocal pending_order, last_heavy
        k = 0
        while k < len(tokens):
            kind, value, count = tokens[k][:3]
            token_charge = tokens[k][3] if len(tokens[k]) > 3 else 0
            if kind == 'open':                              # find its ')' and repeat the inside count times
                depth, j = 1, k + 1
                while j < len(tokens) and depth:
                    depth += {'open': 1, 'close': -1}.get(tokens[j][0], 0)
                    j += 1
                if depth:
                    raise ValueError('unclosed bracket')
                inner, times = tokens[k + 1:j - 1], tokens[j - 1][2]
                # a saturated group branches off the atom before the bracket (N(CH3)2, P(O)(OH)2); one with room
                # left chains on from the copy before ((CH2)3CH3)
                branch = latest_free(scopes[-1])
                if times == 1 and len(inner) == 1 and inner[0][:3] in (['el', 'O', 1], ['el', 'S', 1]) and                         branch is not None and b.free[branch] >= 2:   # C(O)CH3, P(O)(OH)2: (O) is =O
                    a = b.add(inner[0][1])
                    b.bond(branch, a, 2)
                    scopes[-1].append(a)
                    k = j
                    continue
                previous = []
                for _ in range(times):
                    chain = latest_free(previous)
                    entries.append(chain if chain is not None else branch)
                    scopes.append([])
                    run(inner)
                    done = scopes.pop()
                    entries.pop()
                    scopes[-1].extend(done)
                    previous = done
                k = j
                continue
            if kind == 'close':
                raise ValueError('stray bracket')
            if kind == 'bond':
                pending_order = value
            elif kind == 'charge':
                if last_heavy is None:
                    raise ValueError('charge on nothing')
                at = b.mol.GetAtomWithIdx(last_heavy)
                used = _valence(at.GetSymbol()) - b.free[last_heavy]
                at.SetFormalCharge(value)
                b.free[last_heavy] = _valence(at.GetSymbol(), value) - used
            elif kind == 'el' and value == 'H':
                for _ in range(count):                      # the latest atoms with room take the hydrogens
                    a = latest_free(scopes[-1]) if scopes[-1] else None
                    if a is None:
                        a = latest_free([x for s in scopes for x in s])
                    if a is None:
                        raise ValueError('no room for H')
                    b.free[a] -= 1
                    b.hs[a] += 1
            elif kind == 'el':
                scope = scopes[-1]
                onto = latest_free(scope) if count > 1 else None
                for _ in range(count):
                    a = b.add(value, token_charge if _ == count - 1 else 0)
                    here = onto if onto is not None and b.free[onto] >= pending_order else None
                    attach(a, scope, pending_order, here)
                    pending_order = 1
                    if here is not None and value in ('O', 'S'):
                        pendant_of[a] = here
                    host = [n.GetIdx() for n in b.mol.GetAtomWithIdx(a).GetNeighbors()]
                    carbonyl = (count == 1 and value in ('O', 'S') and host and
                                b.mol.GetAtomWithIdx(host[0]).GetSymbol() in ('C', 'P', 'S') and b.free[host[0]] >= 2 and
                                k + 1 < len(tokens) and tokens[k + 1][0] in ('el', 'nick', 'open') and
                                not (tokens[k + 1][0] == 'el' and tokens[k + 1][1] == 'H'))
                    if carbonyl:                            # C=O (P=O, S=O): the chain goes on from that atom
                        bd = b.mol.GetBondBetweenAtoms(host[0], a)
                        bd.SetBondType(Chem.BondType.DOUBLE)
                        b.free[host[0]] -= 1
                        b.free[a] -= 1
                    scope.append(a)
                    last_heavy = a
            elif kind == 'nick':
                sub = Chem.MolFromSmiles(nicknames[value])
                for _ in range(count):
                    base = b.mol.GetNumAtoms()
                    for atom in sub.GetAtoms():
                        new = Chem.Atom(atom.GetSymbol())
                        new.SetFormalCharge(atom.GetFormalCharge())
                        new.SetIsAromatic(atom.GetIsAromatic())
                        b.mol.AddAtom(new)
                        b.free.append(0)
                        b.hs.append(-1)                     # RDKit gives these their hydrogens
                    for bd in sub.GetBonds():
                        b.mol.AddBond(base + bd.GetBeginAtomIdx(), base + bd.GetEndAtomIdx(), bd.GetBondType())
                    b.free[base] = 1
                    attach(base, scopes[-1], pending_order)
                    pending_order = 1
                    b.free[base] = 0
                    scopes[-1].append(base)
            k += 1

    try:
        run(toks)
    except (ValueError, IndexError):
        return None
    n = b.mol.GetNumAtoms()
    if n < 2 or min(b.free) < 0:
        return None
    for order in (2, 3):                                    # leftover free bonds pair up: C=O, S=O, C#N, N=C=O
        changed = True
        while changed:
            changed = False
            for bd in b.mol.GetBonds():
                p, q = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
                sp, sq = bd.GetBeginAtom().GetSymbol(), bd.GetEndAtom().GetSymbol()
                hyper_to_c = (sp in _HYPER and sq == 'C') or (sq in _HYPER and sp == 'C')   # SCN: no S=C
                if bd.GetBondTypeAsDouble() < order and b.free[p] > 0 and b.free[q] > 0 and b.hs[p] >= 0 and                         b.hs[q] >= 0 and not hyper_to_c:
                    bd.SetBondType(Chem.BondType.DOUBLE if bd.GetBondTypeAsDouble() == 1 else Chem.BondType.TRIPLE)
                    b.free[p] -= 1; b.free[q] -= 1
                    changed = True
    for bd in list(b.mol.GetBonds()):                       # SO3K: a potassium sulfonate, O- K+
        p, q = bd.GetBeginAtom(), bd.GetEndAtom()
        for metal, other in ((p, q), (q, p)):
            if metal.GetSymbol() in _SALT_METALS and other.GetSymbol() in ('O', 'S', 'N') and bd.GetBondType() == Chem.BondType.SINGLE:
                b.mol.RemoveBond(metal.GetIdx(), other.GetIdx())
                metal.SetFormalCharge(1); other.SetFormalCharge(other.GetFormalCharge() - 1)
                break
    for a in b.mol.GetAtoms():
        # the hydrogens written; an atom left short gets the rest from RDKit, except P / S / I, which simply take their
        # lower valence (SH, PH2)
        i = a.GetIdx()
        if b.hs[i] >= 0:
            a.SetNumExplicitHs(b.hs[i])
            a.SetNoImplicit(b.free[i] == 0 or a.GetSymbol() in _HYPER)
    mol = b.mol.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol, rootedAtAtom=0, canonical=False)


@lru_cache(maxsize=512)
def _cached(text):
    import model
    return _parse(text, model.FIXED_NICKNAMES)


def smiles(text):
    """The group a condensed label stands for, as a SMILES whose first atom attaches; None when it is not one."""
    text = (text or '').strip().replace('−', '-')
    if not text or len(text) > 40:
        return None
    if text in _ALIASES:
        import model
        return model.FIXED_NICKNAMES[_ALIASES[text]]
    return _cached(text)
