"""
biopolymer.py - peptide, DNA and RNA sequences (like ChemDraw's HELM / biopolymer sequences): the monomers, reading
FASTA, three-letter text (H-Ala-Gly-Pro-NH2) and HELM, writing HELM / FASTA / three-letter text, and the molecule a
sequence stands for (formula, exact mass, m/z of the charge states, expansion into a drawn structure).

    seq = model.Sequence(0, 0, 'PEPTIDE', parse('ACDEF', 'PEPTIDE')[0])
    to_mol(seq)            # RDKit molecule (peptides built residue by residue, so Sec / Pyl, cyclic peptides and
                           # disulfide / lactam bridges work; DNA / RNA through RDKit's MolFromSequence)
    helm(seq), fasta(seq), three_letter(seq)

The monomer colours are this program's own (by residue class), not ChemDraw's palette.
"""
import re

from rdkit import Chem

# one letter: (three letters, name, SMILES N...C(=O)O with N first and the C-terminal OH last, side-chain atom used
# by bridges or None, the side-chain OH that leaves in a lactam or None) - L amino acids
AMINO_ACIDS = {
    'A': ('Ala', 'alanine', 'N[C@@H](C)C(=O)O', None, None),
    'R': ('Arg', 'arginine', 'N[C@@H](CCCNC(=N)N)C(=O)O', None, None),
    'N': ('Asn', 'asparagine', 'N[C@@H](CC(N)=O)C(=O)O', None, None),
    'D': ('Asp', 'aspartic acid', 'N[C@@H](CC(=O)O)C(=O)O', 3, 5),
    'C': ('Cys', 'cysteine', 'N[C@@H](CS)C(=O)O', 3, None),
    'Q': ('Gln', 'glutamine', 'N[C@@H](CCC(N)=O)C(=O)O', None, None),
    'E': ('Glu', 'glutamic acid', 'N[C@@H](CCC(=O)O)C(=O)O', 4, 6),
    'G': ('Gly', 'glycine', 'NCC(=O)O', None, None),
    'H': ('His', 'histidine', 'N[C@@H](Cc1c[nH]cn1)C(=O)O', None, None),
    'I': ('Ile', 'isoleucine', 'N[C@@H]([C@@H](C)CC)C(=O)O', None, None),
    'L': ('Leu', 'leucine', 'N[C@@H](CC(C)C)C(=O)O', None, None),
    'K': ('Lys', 'lysine', 'N[C@@H](CCCCN)C(=O)O', 6, None),
    'M': ('Met', 'methionine', 'N[C@@H](CCSC)C(=O)O', None, None),
    'F': ('Phe', 'phenylalanine', 'N[C@@H](Cc1ccccc1)C(=O)O', None, None),
    'P': ('Pro', 'proline', 'N1[C@@H](CCC1)C(=O)O', None, None),
    'S': ('Ser', 'serine', 'N[C@@H](CO)C(=O)O', None, None),
    'T': ('Thr', 'threonine', 'N[C@@H]([C@@H](C)O)C(=O)O', None, None),
    'W': ('Trp', 'tryptophan', 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O', None, None),
    'Y': ('Tyr', 'tyrosine', 'N[C@@H](Cc1ccc(O)cc1)C(=O)O', None, None),
    'V': ('Val', 'valine', 'N[C@@H](C(C)C)C(=O)O', None, None),
    'U': ('Sec', 'selenocysteine', 'N[C@@H](C[SeH])C(=O)O', 3, None),
    'O': ('Pyl', 'pyrrolysine', 'N[C@@H](CCCCNC(=O)[C@@H]1N=CC[C@H]1C)C(=O)O', None, None),
}
AMBIGUOUS = {'B': 'Asx (Asp or Asn)', 'Z': 'Glx (Glu or Gln)', 'J': 'Xle (Leu or Ile)', 'X': 'Xaa (any amino acid)'}
THREE = {v[0].upper(): k for k, v in AMINO_ACIDS.items()}
THREE.update({'ASX': 'B', 'GLX': 'Z', 'XLE': 'J', 'XAA': 'X', 'UNK': 'X'})

NUCLEOTIDES = {'A': 'adenine', 'C': 'cytosine', 'G': 'guanine', 'T': 'thymine', 'U': 'uracil'}
NUC_AMBIGUOUS = {'N': 'any base', 'R': 'A or G (purine)', 'Y': 'C or T/U (pyrimidine)', 'S': 'G or C', 'W': 'A or T/U',
                 'K': 'G or T/U', 'M': 'A or C', 'B': 'not A', 'D': 'not C', 'H': 'not G', 'V': 'not T/U'}
POLYMERS = ('PEPTIDE', 'DNA', 'RNA')
N_TERMS = {'H': ('H', 'free amine'), 'Ac': ('Ac', 'acetyl')}
C_TERMS = {'OH': ('OH', 'free acid'), 'NH2': ('NH2', 'amide')}

# colours by residue class (this program's own palette)
_CLASS_COLOURS = {'aliphatic': '#E9B872', 'aromatic': '#B49FDC', 'basic': '#7FA7E0', 'acidic': '#E58C8A',
                  'polar': '#8CCB9B', 'cys': '#F2D45C', 'gly': '#C9CDD2', 'pro': '#E8A9C9', 'rare': '#7CC9C4',
                  'unknown': '#D9DDE1'}
_CLASS = {**{k: 'aliphatic' for k in 'AVLIMJ'}, **{k: 'aromatic' for k in 'FWY'}, **{k: 'basic' for k in 'KRH'},
          **{k: 'acidic' for k in 'DEBZ'}, **{k: 'polar' for k in 'STNQ'}, 'C': 'cys', 'G': 'gly', 'P': 'pro',
          'U': 'rare', 'O': 'rare'}
_BASE_COLOURS = {'A': '#8CCB7E', 'C': '#7FA7E0', 'G': '#F0B860', 'T': '#E58C8A', 'U': '#E58C8A'}


def colour(polymer, residue):
    """Fill colour of a residue's shape on the drawing."""
    if polymer == 'PEPTIDE':
        return _CLASS_COLOURS[_CLASS.get(residue, 'unknown')]
    return _BASE_COLOURS.get(residue, _CLASS_COLOURS['unknown'])


def residue_name(polymer, residue):
    if polymer == 'PEPTIDE':
        if residue in AMINO_ACIDS:
            three, name = AMINO_ACIDS[residue][:2]
            return f'{three} ({residue}), {name}'
        return AMBIGUOUS.get(residue, residue)
    sugar = '2′-deoxyribo' if polymer == 'DNA' else 'ribo'
    if residue in NUCLEOTIDES:
        return f'{residue}, {NUCLEOTIDES[residue]} {sugar}nucleotide'
    return f'{residue}, {NUC_AMBIGUOUS.get(residue, "unknown base")}'


def label(seq, k):
    """The text shown on residue k: one letter, or three letters in the three-letter style (peptides)."""
    r = seq.residues[k]
    if seq.polymer == 'PEPTIDE' and seq.three_letter:
        return AMINO_ACIDS[r][0] if r in AMINO_ACIDS else {'B': 'Asx', 'Z': 'Glx', 'J': 'Xle'}.get(r, 'Xaa')
    return r


# ------------------------------------------------------------------------------------------------ reading
def parse(text, polymer=None):
    """(polymer, residues, extras) from FASTA (one letter, '>' header lines and digits / spaces ignored), three-letter
    text with separators ('H-Ala-Gly-Pro-NH2', 'Ala Gly Pro') or HELM ('PEPTIDE1{A.C.D}$$$$', 'RNA1{R(A)P.R(C)}$$$$');
    polymer None guesses (HELM says it; otherwise peptide unless only A C G T / U). extras: n_term, c_term, cyclic,
    bridges (0-based residue pairs) when the text gives them. Raises ValueError with what is wrong."""
    text = text.strip()
    if not text:
        raise ValueError('The sequence is empty.')
    if re.match(r'^(PEPTIDE|RNA|CHEM|BLOB)\d+\{', text, re.I):
        return _parse_helm(text)
    extras = {}
    words = [w for w in re.split(r'[\s\-–·.,;]+', text) if w]
    if len(words) > 1 and words[0] in ('H', 'Ac', 'H2N'):
        extras['n_term'] = 'Ac' if words[0] == 'Ac' else 'H'
        words = words[1:]
    if len(words) > 1 and words[-1] in ('OH', 'NH2', 'COOH'):
        extras['c_term'] = 'NH2' if words[-1] == 'NH2' else 'OH'
        words = words[:-1]
    if (extras or len(words) > 1) and all(len(w) >= 3 for w in words) and any(w.upper() in THREE for w in words):
        residues = []
        for w in words:
            if w.upper() not in THREE:
                raise ValueError(f'"{w}" is not an amino acid (three-letter codes: Ala, Gly, Pro ...).')
            residues.append(THREE[w.upper()])
        return 'PEPTIDE', residues, extras
    lines =[ln for ln in text.splitlines() if not ln.lstrip().startswith(('>', ';'))]
    letters = re.sub(r'[\s\d\-*]', '', ''.join(lines)).upper()
    if not letters:
        raise ValueError('The sequence is empty.')
    if polymer is None:
        polymer = 'PEPTIDE'
        if set(letters) <= set('ACGTN'):
            polymer = 'DNA'
        elif set(letters) <= set('ACGUN'):
            polymer = 'RNA'
    if polymer == 'PEPTIDE':
        bad = sorted(set(letters) - set(AMINO_ACIDS) - set(AMBIGUOUS))
        if bad:
            raise ValueError('Not amino acid letters: ' + ', '.join(bad))
    else:
        swap = {'U': 'T'} if polymer == 'DNA' else {'T': 'U'}  # as ChemDraw: U pasted as DNA becomes T, and back
        letters = ''.join(swap.get(c, c) for c in letters)
        bad = sorted(set(letters) - set(NUCLEOTIDES) - set(NUC_AMBIGUOUS))
        if bad:
            raise ValueError('Not nucleotide letters: ' + ', '.join(bad))
    return polymer, list(letters), extras


def _parse_helm(text):
    """One linear PEPTIDE or RNA polymer of a HELM string, with its self-connections (cyclic, bridges)."""
    parts = text.split('$')
    m = re.match(r'^(PEPTIDE|RNA)(\d+)\{(.*?)\}$', parts[0].strip(), re.I)
    if not m:
        raise ValueError('Only one PEPTIDE or RNA polymer per HELM string is read.')
    kind, name, body = m.group(1).upper(), m.group(1).upper() + m.group(2), m.group(3)
    units = [u for u in body.split('.') if u]
    extras = {}
    if kind == 'PEPTIDE':
        residues = []
        for u in units:
            u = u.strip('[]')
            if u in AMINO_ACIDS or u in AMBIGUOUS:
                residues.append(u)
            elif u.upper() in THREE:
                residues.append(THREE[u.upper()])
            else:
                raise ValueError(f'The monomer "{u}" is not in this program\'s library (natural amino acids only).')
        polymer = 'PEPTIDE'
    else:
        residues, deoxy = [], False
        for u in units:
            bm = re.search(r'\((\w)\)', u)
            if bm:
                residues.append(bm.group(1).upper())
                deoxy = deoxy or 'dR' in u
            elif u.strip('[]') not in ('P', 'R', 'dR'):
                raise ValueError(f'The monomer "{u}" is not in this program\'s library (A, C, G, T, U).')
        polymer = 'DNA' if deoxy or 'T' in residues else 'RNA'
        if units and units[0].strip('[]') == 'P':
            extras['p5'] = True
        if units and units[-1].strip('[]') == 'P':
            extras['p3'] = True
    if len(parts) > 1 and parts[1].strip():                  # connections: NAME,NAME,i:Rx-j:Ry|...
        bridges = []
        for conn in parts[1].split('|'):
            cm = re.match(r'^\s*(\w+),(\w+),(\d+):R(\d)-(\d+):R(\d)\s*$', conn)
            if not cm or cm.group(1).upper() != name or cm.group(2).upper() != name:
                continue
            i, ri, j, rj = int(cm.group(3)) - 1, int(cm.group(4)), int(cm.group(5)) - 1, int(cm.group(6))
            if {ri, rj} == {1, 2} and {i, j} == {0, len(residues) - 1}:
                extras['cyclic'] = True
            elif ri == 3 and rj == 3:
                bridges.append([min(i, j), max(i, j)])
        if bridges:
            extras['bridges'] = bridges
    return polymer, residues, extras


# ------------------------------------------------------------------------------------------------ writing
def fasta(seq, width=60):
    """FASTA text (one-letter codes; the natural analogs), lines of width letters."""
    s = ''.join(seq.residues)
    return '\n'.join(s[k:k + width] for k in range(0, len(s), width))


def three_letter(seq):
    """'H-Ala-Gly-Pro-NH2' (peptides) or "5'-ACGT-3'" (DNA / RNA)."""
    if seq.polymer != 'PEPTIDE':
        return "5′-" + ''.join(seq.residues) + "-3′"
    body = '-'.join(label_three(r) for r in seq.residues)
    text = f'{seq.n_term}-{body}-{seq.c_term}'
    return f'cyclo({body})' if seq.cyclic else text


def label_three(r):
    return AMINO_ACIDS[r][0] if r in AMINO_ACIDS else {'B': 'Asx', 'Z': 'Glx', 'J': 'Xle'}.get(r, 'Xaa')


def helm(seq):
    """HELM string of the sequence (HELM 1 notation, as ChemDraw's Copy As > HELM writes it)."""
    if seq.polymer == 'PEPTIDE':
        name = 'PEPTIDE1'
        body = '.'.join(seq.residues)
        conns = []
        if seq.cyclic and len(seq.residues) > 1:
            conns.append(f'{name},{name},{len(seq.residues)}:R2-1:R1')
        for i, j in seq.bridges:
            conns.append(f'{name},{name},{i + 1}:R3-{j + 1}:R3')
        return f'{name}{{{body}}}$' + '|'.join(conns) + '$$$'
    sugar = '[dR]' if seq.polymer == 'DNA' else 'R'
    units = [f'{sugar}({b})' for b in seq.residues]
    body = 'P.'.join(units)                                 # R(A)P.R(C)P.R(G)
    if seq.p5:
        body = 'P.' + body
    if seq.p3:
        body += 'P'
    return f'RNA1{{{body}}}$$$$'


# ------------------------------------------------------------------------------------------------ the molecule
def to_mol(seq):
    """(RDKit molecule or None, error text) of the sequence with its termini, cyclisation and bridges."""
    if not seq.residues:
        return None, 'The sequence is empty.'
    if seq.polymer == 'PEPTIDE':
        return _peptide(seq)
    unknown = sorted(set(seq.residues) - set(NUCLEOTIDES))
    if unknown:
        return None, 'Ambiguous bases (' + ', '.join(unknown) + ') have no structure.'
    flavor = (6 if seq.polymer == 'DNA' else 2) + (1 if seq.p5 else 0) + (2 if seq.p3 else 0)
    mol = Chem.MolFromSequence(''.join(seq.residues), flavor=flavor)
    if mol is None:
        return None, 'The sequence could not be built.'
    return mol, ''


def _peptide(seq):
    unknown = sorted(set(seq.residues) - set(AMINO_ACIDS))
    if unknown:
        return None, 'Ambiguous residues (' + ', '.join(unknown) + ') have no structure.'
    rw = Chem.RWMol()
    starts, n_atom, c_atom, oh_atom, side, side_oh = [], [], [], [], [], []
    for r in seq.residues:
        _three, _name, smiles, s_atom, s_oh = AMINO_ACIDS[r]
        frag = Chem.MolFromSmiles(smiles)
        base = rw.GetNumAtoms()
        rw.InsertMol(frag)
        n = frag.GetNumAtoms()
        starts.append(base)
        n_atom.append(base)
        c_atom.append(base + n - 3)
        oh_atom.append(base + n - 1)
        side.append(None if s_atom is None else base + s_atom)
        side_oh.append(None if s_oh is None else base + s_oh)
    remove = set()
    for k in range(len(seq.residues) - 1):                  # peptide bonds: C(=O) of k to N of k + 1, its OH leaves
        rw.AddBond(c_atom[k], n_atom[k + 1], Chem.BondType.SINGLE)
        remove.add(oh_atom[k])
    last = len(seq.residues) - 1
    if seq.cyclic and last > 0:
        rw.AddBond(c_atom[last], n_atom[0], Chem.BondType.SINGLE)
        remove.add(oh_atom[last])
    else:
        if seq.c_term == 'NH2':
            rw.GetAtomWithIdx(oh_atom[last]).SetAtomicNum(7)
        if seq.n_term == 'Ac':
            c = rw.AddAtom(Chem.Atom(6)); o = rw.AddAtom(Chem.Atom(8)); me = rw.AddAtom(Chem.Atom(6))
            rw.AddBond(n_atom[0], c, Chem.BondType.SINGLE)
            rw.AddBond(c, o, Chem.BondType.DOUBLE)
            rw.AddBond(c, me, Chem.BondType.SINGLE)
    for i, j in seq.bridges:
        if not (0 <= i < len(seq.residues) and 0 <= j < len(seq.residues)) or i == j:
            return None, f'Bridge {i + 1}-{j + 1} is outside the sequence.'
        a, b = seq.residues[i], seq.residues[j]
        if a in 'CU' and b in 'CU':                         # disulfide (or Se-S / diselenide)
            rw.AddBond(side[i], side[j], Chem.BondType.SINGLE)
        elif {a, b} in ({'K', 'D'}, {'K', 'E'}):            # lactam: Lys NZ to the Asp / Glu side-chain carbonyl
            k_, acid = (i, j) if a == 'K' else (j, i)
            rw.AddBond(side[k_], side[acid], Chem.BondType.SINGLE)
            remove.add(side_oh[acid])
        else:
            return None, (f'A bridge joins Cys-Cys (disulfide) or Lys-Asp / Lys-Glu (lactam), not '
                          f'{label_three(a)}{i + 1}-{label_three(b)}{j + 1}.')
    for idx in sorted(remove, reverse=True):
        rw.RemoveAtom(idx)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception as e:
        return None, str(e).split('\n')[0]
    return mol, ''


def analyse(seq):
    """chem.Info of the sequence, plus the protonated charge states [M+zH]z+ / deprotonated [M-zH]z- that ESI gives
    for peptides and oligonucleotides (z up to the number of basic / acidic sites + 1, at most 10)."""
    import chem
    mol, err = to_mol(seq)
    if mol is None:
        info = chem.Info()
        info.errors.append(err)
        return info
    info = chem.analyse_mol(mol)
    if seq.polymer == 'PEPTIDE':
        zmax = min(10, 1 + sum(r in 'KRH' for r in seq.residues) + (0 if seq.cyclic or seq.n_term == 'Ac' else 1))
        sign = 1
    else:
        zmax = min(10, len(seq.residues))
        sign = -1
    have = {ion.label for ion in info.ions}
    for z in range(2, zmax + 1):
        q = sign * z
        lab = chem.ion_label(f'M{"+" if q > 0 else "-"}{z}H', q)
        if lab in have:
            continue
        mz = (info.exact + q * chem.H - q * chem.ELECTRON) / z
        counts = dict(info.counts)
        counts['H'] = counts.get('H', 0) + q
        info.ions.append(chem.Ion(lab, q, mz, counts))
    return info


def summary(seq):
    """'Peptide, 12 residues' / 'DNA, 20 nt' for lists and tooltips."""
    n = len(seq.residues)
    if seq.polymer == 'PEPTIDE':
        return f'Peptide, {n} residue{"s" if n != 1 else ""}' + (', cyclic' if seq.cyclic else '')
    return f'{seq.polymer}, {n} nt'
