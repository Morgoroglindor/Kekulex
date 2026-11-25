"""
fragments.py - the mass-spectrum fragments a structure can logically give, from the textbook fragmentation rules:

    pre, frags = fragments.predict(doc, atom_ids, mode='ESI+', steps=2)

mode: 'EI' (the radical cation M•+), 'ESI+' ([M+H]+) or 'ESI-' ([M−H]−). A structure drawn as an ion is used as it is.

EI (odd-electron M•+): α-cleavage next to N, O, S (iminium / oxonium ions) and next to C=O (acylium ions), benzylic and
allylic cleavage (C7H7+ tropylium), inductive cleavage at C–halogen and C–O, loss of HX, the McLafferty rearrangement
(the γ-H moves to the C=O oxygen), retro-Diels–Alder of a six-membered ring with a double bond, and the usual neutral
losses (H2O, CH3•, ketene, NO2•, CO, HCN ...).
ESI (even-electron [M+H]+ / [M−H]−): heterolytic cleavage of C–O, C–N and C–S bonds, either the carbon side as a cation
(acylium / b-type, benzylic, sugar oxonium) or the heteroatom side with the charge and an H moved to it (aglycone Y ions,
y-type amines, [acid+H]+, carboxylates, quinate ...), retro-Diels–Alder (flavonoid 1,3A ions), and neutral losses
(H2O, NH3, CH3OH, HCOOH, CO2, CO, C4H8, ketene, SO3, H3PO4, CH3• from aryl methoxy groups, C-glycoside cross-ring 90 / 120).
With steps > 1 each fragment loses the neutral molecules its own atoms allow (and an acylium ion CO).

Electron bookkeeping: a fragment is a set of atoms, the bonds to the rest broken and capped with H, then marks on some
atoms (charge, unpaired electron, H gained or lost); its formula and m/z follow, and the lost neutral is the precursor
minus the fragment. The same marks draw it (draw_doc), so the drawn fragment has the listed formula.
Likelihood is a rough plausibility from the rule (and the stability of the ion it gives), not a predicted intensity.
"""
from dataclasses import dataclass, field

from rdkit import Chem

import chem

MODES = (('ESI+', 'ESI+  [M+H]+'), ('ESI-', 'ESI−  [M−H]−'), ('EI', 'EI  M•+'))

REFS = {
    'ei': 'McLafferty FW, Tureček F. Interpretation of Mass Spectra, 4th ed. University Science Books, 1993.',
    'mclafferty': 'McLafferty FW. Anal. Chem. 1959, 31, 82-87; McLafferty & Tureček, Interpretation of Mass '
                  'Spectra, 1993.',
    'esi': 'Demarque DP et al. Fragmentation reactions using electrospray ionization mass spectrometry. Nat. Prod. '
           'Rep. 2016, 33, 432-455 (doi:10.1039/C5NP00073D).',
    'even': 'Karni M, Mandelbaum A. The even-electron rule. Org. Mass Spectrom. 1980, 15, 53-64.',
    'rda': 'Fabre N et al. J. Am. Soc. Mass Spectrom. 2001, 12, 707-715 (flavonoid 1,3A / 1,3B ions); '
           'Ma YL et al. Rapid Commun. Mass Spectrom. 1997, 11, 1357-1364.',
    'glyco': 'Domon B, Costello CE. Glycoconjugate J. 1988, 5, 397-409 (B / Y ions); Cuyckens F, Claeys M. '
             'J. Mass Spectrom. 2004, 39, 1-15.',
    'radical': 'Hvattum E, Ekeberg D. J. Mass Spectrom. 2003, 38, 43-49 (radical aglycone [Y0-H]•−).',
    'peptide': 'Roepstorff P, Fohlman J. Biomed. Mass Spectrom. 1984, 11, 601 (b / y ions); Paizs B, Suhai S. '
               'Mass Spectrom. Rev. 2005, 24, 508-548 (mobile proton, immonium ions).',
    'lipid': 'Hsu FF, Turk J. J. Am. Soc. Mass Spectrom. 2003, 14, 352-363 (phosphocholine m/z 184; PE loss of 141).',
    'tms': 'Harvey DJ, Vouros P. Mass Spectrom. Rev. 2020, 39, 105-211 (TMS derivatives: m/z 73, 147, M−15).',
    'chlorogenic': 'Clifford MN et al. J. Agric. Food Chem. 2003, 51, 2900-2911 (quinate m/z 191, caffeate 179).',
}

LOSS_TEXT = {
    'alcohol': ('Protonation of the OH makes water a good leaving group: the C–O bond pair goes to O and H2O '
                'leaves, giving the carbocation (or, with a neighbouring H, a double bond).',
                'The radical cation loses water: a hydrogen (often from C-4, through a six-membered ring) moves to '
                'the OH as the C–O bond breaks (1,2- / 1,4-elimination).'),
    'carboxylic acid': ('The protonated acid loses water from its OH (giving the acylium ion), or CO2 / HCOOH; '
                        'the carboxylate anion loses CO2, its C–C bond pair staying on the carbon.', None),
    'primary amine': ('The protonated NH3+ group leaves as ammonia: the C–N bond pair goes to N, leaving a '
                      'carbocation (or an alkene after a 1,2-H shift).', None),
    'primary amide': ('The protonated amide loses ammonia, giving the acylium ion.', None),
    'methyl ester': ('Protonation at the ester O-methyl oxygen; the acyl–O bond pair goes to O and methanol '
                     'leaves, giving the acylium ion.', None),
    'methyl ether': ('The protonated methoxy group leaves as methanol, giving the carbocation.', None),
    'aryl methoxy group': ('Homolysis of the O–CH3 bond: a methyl radical leaves and the stable phenoxyl-type '
                           'radical ion remains (an exception to the even-electron rule, typical of methoxylated '
                           'flavonoids and phenolic acids).', None),
    'tert-butyl group': ('Protonation at the heteroatom; the O/N–C(CH3)3 bond pair goes to the heteroatom and the '
                         'tert-butyl cation leaves as isobutene after losing a proton (E1).', None),
    'acetyl group, as ketene': ('Four-centre elimination: an H of the acetyl CH3 moves to the O / N as the acyl–'
                                'O / N bond breaks, and ketene (CH2=C=O) leaves.', None),
    'alkyl chloride': ('A hydrogen on the next carbon leaves with the halogen (1,2-elimination of HX).', None),
    'alkyl bromide': ('A hydrogen on the next carbon leaves with the halogen (1,2-elimination of HX).', None),
    'alkyl fluoride': ('A hydrogen on the next carbon leaves with the halogen (1,2-elimination of HX).', None),
    'nitro group': ('The C–NO2 bond breaks homolytically: NO2• leaves and the aryl cation (or anion) remains.', None),
    'thiol': ('The protonated SH leaves as H2S, as water leaves an alcohol.', None),
    'sulfate / sulfonic acid': ('The S–O bond pair stays on the oxygen and SO3 leaves (sulfates: [M−H−80]−).', None),
    'phosphate ester': ('The protonated phosphate leaves as H3PO4 (with a neighbouring H, an alkene forms), or HPO3 '
                        'is expelled and the alcohol remains (phosphopeptides: −98 / −80).', None),
    'trimethylsilyl group': ('A methyl radical leaves from silicon (α-cleavage at Si): the [M−CH3]+ ion of TMS '
                             'derivatives, often stronger than M•+.', None),
    'tert-butyldimethylsilyl group': ('The tert-butyl radical leaves from silicon: [M−57]+, usually the base peak '
                                      'of TBDMS derivatives.', None),
}

LIKELY, POSSIBLE = 0.65, 0.40


@dataclass
class Fragment:
    mz: float
    z: int
    odd: bool                                               # odd-electron (radical) ion
    counts: dict
    ion: str                                                # 'C7H7+', 'C3H6O•+'
    lost: str                                               # 'C3H7•', 'H2O + CO' (from the precursor)
    how: str
    score: float
    atoms: frozenset = frozenset()                          # drawing atom ids in the ion
    lost_atoms: frozenset = frozenset()                     # drawing atom ids that leave
    step: int = 1
    also: list = field(default_factory=list)                # other rules giving the same ion
    S: frozenset = None                                     # expanded-drawing atom ids (None: formula only)
    marks: dict = None                                      # {expanded atom id: (charge, radicals, H)}
    mech: str = ''                                          # the mechanism of this step, in words
    arrows: list = field(default_factory=list)              # its curved arrows: ('full' | 'half', from, to), each end
                                                            # ('atom', id) | ('bond', (id, id)) | ('H', id) | ('bondH', id)
    pre: list = field(default_factory=list)                 # where the charge sits before it: [(atom, charge, •, H)]
    ref: str = ''                                           # literature for the rule
    parent: object = None                                   # the fragment it comes from (None: the precursor ion)

    @property
    def likelihood(self):
        return 'likely' if self.score >= LIKELY else 'possible' if self.score >= POSSIBLE else 'minor'

    @property
    def formula(self):
        return chem.formula_text(self.counts)


def ion_text(counts, z, odd):
    return chem.formula_text(counts) + ('•' if odd else '') + chem.charge_text(z)


def _sub(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) - v
    return out


def _clean(c):
    return {k: v for k, v in c.items() if v}


# (name, SMARTS, atoms that leave, attach atom, H that leave with them, radical, {mode: likelihood})
# formula-only losses (no atoms): (name, test, formula, radical, {mode: likelihood})
GROUP_LOSSES = (
    ('H2O', '[CX4;!$(C=O)][OX2H1]', (1,), 0, 2, False, {'EI': 0.45, 'ESI+': 0.75, 'ESI-': 0.35}, 'alcohol'),
    ('H2O', '[CX3](=O)[OX2H1]', (2,), 0, 2, False, {'ESI+': 0.55, 'ESI-': 0.25}, 'carboxylic acid'),
    ('CO2', '[#6][CX3](=O)[OX2H1,OX1-]', (1, 2, 3), 0, 0, False, {'EI': 0.15, 'ESI+': 0.2, 'ESI-': 0.9},
     'carboxylic acid'),
    ('HCOOH', '[#6][CX3](=O)[OX2H1]', (1, 2, 3), 0, 2, False, {'ESI+': 0.45, 'EI': 0.15}, 'carboxylic acid'),
    ('NH3', '[CX4][NX3H2;!$(N[#6]=[O,S,N])]', (1,), 0, 3, False, {'EI': 0.2, 'ESI+': 0.75}, 'primary amine'),
    ('NH3', '[CX3](=O)[NX3H2]', (2,), 0, 3, False, {'ESI+': 0.5}, 'primary amide'),
    ('CH3OH', '[CX3](=O)[OX2][CH3]', (2, 3), 0, 4, False, {'ESI+': 0.65, 'ESI-': 0.3}, 'methyl ester'),
    ('CH3OH', '[CX4][OX2][CH3]', (1, 2), 0, 4, False, {'ESI+': 0.4}, 'methyl ether'),
    ('CH3•', 'c[OX2][CH3]', (2,), 1, 3, True, {'EI': 0.55, 'ESI+': 0.55, 'ESI-': 0.6}, 'aryl methoxy group'),
    ('C4H8', '[O,N,S]C([CH3])([CH3])[CH3]', (1, 2, 3, 4), 0, 8, False, {'EI': 0.5, 'ESI+': 0.8, 'ESI-': 0.3},
     'tert-butyl group'),
    ('C2H2O', '[O,N;!$([O,N]C(=O)C(=O))]C(=O)[CH3]', (1, 2, 3), 0, 2, False, {'EI': 0.75, 'ESI+': 0.7, 'ESI-': 0.4},
     'acetyl group, as ketene'),
    ('HCl', '[CX4;!H0,$(C[C;!H0])][Cl]', (1,), 0, 1, False, {'EI': 0.45, 'ESI+': 0.4, 'ESI-': 0.2}, 'alkyl chloride'),
    ('HBr', '[CX4;!H0,$(C[C;!H0])][Br]', (1,), 0, 1, False, {'EI': 0.3, 'ESI+': 0.3}, 'alkyl bromide'),
    ('HF', '[CX4;!H0,$(C[C;!H0])][F]', (1,), 0, 1, False, {'EI': 0.35, 'ESI+': 0.4}, 'alkyl fluoride'),
    ('NO2•', '[#6][$([NX3](=O)=O),$([NX3+](=O)[O-])](~[OX1])~[OX1]', (1, 2, 3), 0, 0, True,
     {'EI': 0.6, 'ESI-': 0.5}, 'nitro group'),
    ('H2S', '[CX4][SX2H1]', (1,), 0, 2, False, {'EI': 0.35, 'ESI+': 0.6}, 'thiol'),
    ('SO3', '[O,N,#6]S(=O)(=O)[OX2H1,OX1-]', (1, 2, 3, 4), 0, 0, False, {'ESI+': 0.6, 'ESI-': 0.85},
     'sulfate / sulfonic acid'),
    ('H3PO4', '[#6]O[PX4](=O)([OX2H1,OX1-])[OX2H1,OX1-]', (1, 2, 3, 4, 5), 0, 3, False, {'ESI+': 0.8, 'ESI-': 0.4},
     'phosphate ester'),
    ('HPO3', '[#6]O[PX4](=O)([OX2H1,OX1-])[OX2H1,OX1-]', (2, 3, 4, 5), 1, 1, False, {'ESI-': 0.6, 'ESI+': 0.3},
     'phosphate ester'),
    ('CH3•', '[Si]([CH3])([CH3])[CH3]', (1,), 0, 3, True, {'EI': 0.75}, 'trimethylsilyl group'),
    ('C4H9•', '[Si]([CH3])([CH3])C([CH3])([CH3])[CH3]', (3, 4, 5, 6), 0, 9, True, {'EI': 0.9},
     'tert-butyldimethylsilyl group'),
)

# diagnostic ions: (ion formula, odd electron, SMARTS that must match (count), {mode: likelihood}, what, text, ref)
DIAGNOSTIC = (
    ({'C': 3, 'H': 9, 'Si': 1}, False, ('[O,N,S,#6][Si]([CH3])([CH3])[CH3]', 1), {'EI': 0.85},
     'trimethylsilyl ion Si(CH3)3+ (m/z 73)', 'Cleavage of the X–Si bond with the charge on silicon: the '
     'trimethylsilyl cation, present in the spectrum of every TMS derivative.', 'tms'),
    ({'C': 5, 'H': 15, 'O': 1, 'Si': 2}, False, ('O[Si]([CH3])([CH3])[CH3]', 2), {'EI': 0.6},
     'pentamethyldisiloxanyl ion (m/z 147)', 'Two TMS groups close in space: a silyl group migrates to the other '
     'TMS oxygen and (CH3)2Si=O+–Si(CH3)3 forms; it shows two or more TMS–O groups.', 'tms'),
    ({'C': 5, 'H': 15, 'N': 1, 'O': 4, 'P': 1}, False, ('C[N+](C)(C)CCOP(=O)([O-,OX2H1])O', 1), {'ESI+': 0.95},
     'phosphocholine ion (m/z 184)', 'The protonated phosphate of a phosphatidylcholine / sphingomyelin breaks at '
     'the glycerol C–O bond; the charge stays on the phosphocholine head group. Precursor scans of m/z 184 find '
     'choline lipids.', 'lipid'),
    ({'C': 8, 'H': 14, 'N': 1, 'O': 5}, False, ('[CX4;R1]1([OX2])[CX4;R1](NC(=O)[CH3])[CX4;R1][CX4;R1][CX4;R1]O1', 1),
     {'ESI+': 0.8}, 'HexNAc oxonium ion (m/z 204)', 'Glycosidic cleavage with the charge on the N-acetylhexosamine '
     '(B1 oxonium ion); it then loses water (186, 168) and C2H4O2 (138, 126) - the marker of glycopeptides.', 'glyco'),
)


FORMULA_TEXT = {
    'ring C=O (from the ring)': 'The ring opens next to the C=O and carbon monoxide is expelled (ring contraction): '
                                'flavones, chromones, xanthones, quinones.',
    'phenol (from the ring)': 'The phenol ring keto-tautomerises, then expels CO (the cyclopentadiene ion '
                              'remains).',
    'nitro group (as the nitrite)': 'The nitro group rearranges to the nitrite Ar–O–N=O, which loses NO• (then '
                                    'often CO).',
    'aromatic N': 'The ring opens and HCN is expelled (pyridines, anilines, indoles).',
    'C-glycoside (cross-ring, 0,3X)': 'C-glycosides cannot lose the sugar whole; the sugar ring breaks across two '
                                      'bonds instead (0,3X: −90; 0,2X: −120).',
    'C-glycoside (cross-ring, 0,2X)': 'C-glycosides cannot lose the sugar whole; the sugar ring breaks across two '
                                      'bonds instead (0,2X: −120; 0,3X: −90).',
    'lactone / coumarin (from the ring)': 'The lactone ring opens and CO2 leaves (coumarins also lose CO first).',
    'chromone C ring (from the ring)': 'The C ring of the flavone / chromone opens and loses CO2.',
    'chromone C ring, as ketene (from the ring)': 'The C ring loses ketene (C2H2O, −42), as flavones do in CID.',
    'ortho effect': 'Ortho effect: the neighbouring OH / NH2 gives its hydrogen to the ester / acid oxygen through a '
                    'six-membered ring, and the alcohol (or water) leaves; meta and para isomers do not do this.',
    'aryl methyl ether': 'The methyl H moves to the ring and formaldehyde (CH2O) leaves (anisoles: M−30).',
    'phenol, as formyl': 'After ring opening the phenol loses the formyl radical (M−29).',
    'phosphocholine': 'The phosphocholine head group leaves as a neutral (−183), the diacylglycerol ion remains.',
    'phosphoethanolamine': 'Protonated phosphatidylethanolamines lose the phosphoethanolamine head group (−141).',
}

FORMULA_LOSSES = (
    ('CO', '[#6;R][#6X3;R](=O)[#6;R]', {'C': 1, 'O': 1}, False, {'EI': 0.45, 'ESI+': 0.55, 'ESI-': 0.4},
     'ring C=O (from the ring)'),
    ('CO2', '[#6;R](=O)[O,o;R]', {'C': 1, 'O': 2}, False, {'EI': 0.45, 'ESI+': 0.35, 'ESI-': 0.5},
     'lactone / coumarin (from the ring)'),
    ('CO2', 'O=[#6;R]1[#6;R]~[#6;R][O,o;R][#6;R]~[#6;R]1', {'C': 1, 'O': 2}, False, {'ESI+': 0.2, 'ESI-': 0.4},
     'chromone C ring (from the ring)'),
    ('C2H2O', 'O=[#6;R]1[#6;R]~[#6;R][O,o;R][#6;R]~[#6;R]1', {'C': 2, 'H': 2, 'O': 1}, False, {'ESI+': 0.35},
     'chromone C ring, as ketene (from the ring)'),
    ('CO', 'c[OX2H1]', {'C': 1, 'O': 1}, False, {'EI': 0.4, 'ESI+': 0.2, 'ESI-': 0.3}, 'phenol (from the ring)'),
    ('NO•', '[#6][$([NX3](=O)=O),$([NX3+](=O)[O-])]', {'N': 1, 'O': 1}, True, {'EI': 0.4, 'ESI-': 0.35},
     'nitro group (as the nitrite)'),
    ('HCN', '[n,$([NX3H2]c)]', {'C': 1, 'N': 1, 'H': 1}, False, {'EI': 0.35, 'ESI+': 0.25}, 'aromatic N'),
    ('C3H6O3', '[c][CX4;R1]1[OX2;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1', {'C': 3, 'H': 6, 'O': 3}, False,
     {'ESI+': 0.55, 'ESI-': 0.7}, 'C-glycoside (cross-ring, 0,3X)'),
    ('C4H8O4', '[c][CX4;R1]1[OX2;R1][CX4;R1][CX4;R1][CX4;R1][CX4;R1]1', {'C': 4, 'H': 8, 'O': 4},
     False, {'ESI+': 0.55, 'ESI-': 0.75}, 'C-glycoside (cross-ring, 0,2X)'),
    ('CH3OH', 'c([OX2H1,NX3H2,NX3H1])c[CX3](=O)[OX2][CH3]', {'C': 1, 'H': 4, 'O': 1}, False, {'EI': 0.8},
     'ortho effect'),
    ('C2H5OH', 'c([OX2H1,NX3H2,NX3H1])c[CX3](=O)[OX2][CH2][CH3]', {'C': 2, 'H': 6, 'O': 1}, False, {'EI': 0.8},
     'ortho effect'),
    ('H2O', 'c([OX2H1,NX3H2,NX3H1,CH3])c[CX3](=O)[OX2H1]', {'H': 2, 'O': 1}, False, {'EI': 0.75}, 'ortho effect'),
    ('CH2O', 'c[OX2][CH3]', {'C': 1, 'H': 2, 'O': 1}, False, {'EI': 0.4}, 'aryl methyl ether'),
    ('CHO•', 'c[OX2H1]', {'C': 1, 'H': 1, 'O': 1}, True, {'EI': 0.3}, 'phenol, as formyl'),
    ('C5H14NO4P', 'C[N+](C)(C)CCOP(=O)([O-,OX2H1])OC', {'C': 5, 'H': 14, 'N': 1, 'O': 4, 'P': 1}, False,
     {'ESI+': 0.5}, 'phosphocholine'),
    ('C2H8NO4P', '[NX3H2]CCOP(=O)([OX2H1,O-])OC', {'C': 2, 'H': 8, 'N': 1, 'O': 4, 'P': 1}, False, {'ESI+': 0.85},
     'phosphoethanolamine'),
)


class _Mol:
    """The drawing's structure as an RDKit molecule, its atoms mapped back to the drawing."""

    def __init__(self, doc, atom_ids=None):
        full, grown = doc.expanded()
        self.full = full
        self.grown = grown
        self.owner = {a: a for a in doc.atoms}
        for a, more in grown.items():
            for g in more:
                self.owner[g] = a
        ids = [a for a in doc.atoms if atom_ids is None or a in atom_ids]
        ids = {a for a in chem.expand_ids(grown, ids) if a in full.atoms}
        self.mol, self.order, self.error = chem.to_mol(full, ids)
        if self.mol is None:
            return
        m = self.mol
        self.n = m.GetNumAtoms()
        self.el = [a.GetSymbol() for a in m.GetAtoms()]
        self.key = []
        for a in m.GetAtoms():
            iso = a.GetIsotope()
            self.key.append(chem.iso_key(a.GetSymbol(), iso) if iso else a.GetSymbol())
        self.mass = [chem.key_mass(k) or chem.key_mass(e) for k, e in zip(self.key, self.el)]
        self.hs = [a.GetTotalNumHs() for a in m.GetAtoms()]
        self.q = [a.GetFormalCharge() for a in m.GetAtoms()]
        self.rad = [a.GetNumRadicalElectrons() for a in m.GetAtoms()]
        self.nb = [[n.GetIdx() for n in a.GetNeighbors()] for a in m.GetAtoms()]
        self.kek = Chem.Mol(m)                              # bond orders as single / double (pyranone rings ...)
        try:
            Chem.Kekulize(self.kek, clearAromaticFlags=False)
        except Exception:
            self.kek = m
        # aromatic as a chemist sees it: rings of aromatic atoms without an exocyclic C=O (benzene, pyridine, furan;
        # not the pyranone ring of a flavone or coumarin, whose bonds break like ordinary ones)
        self.arobond = set()
        aro = set()
        for ring in m.GetRingInfo().AtomRings():
            if not all(m.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                continue
            if any(self.el[j] in ('O', 'S', 'N') and j not in ring and
                   self.kek.GetBondBetweenAtoms(i, j).GetBondType() == Chem.BondType.DOUBLE
                   for i in ring for j in self.nb[i]):
                continue
            aro |= set(ring)
            for k in range(len(ring)):
                self.arobond.add(frozenset((ring[k], ring[(k + 1) % len(ring)])))
        self.aro = [i in aro for i in range(self.n)]
        self._sides = {}

    # ------------------------------------------------------------------------------------------- bonds and pieces
    def bond(self, i, j):
        return self.kek.GetBondBetweenAtoms(i, j)

    def aromatic_bond(self, i, j):
        return frozenset((i, j)) in self.arobond

    def single(self, i, j, ring=False):
        b = self.bond(i, j)
        return (b is not None and b.GetBondType() == Chem.BondType.SINGLE and (ring or not b.IsInRing())
                and not self.aromatic_bond(i, j))

    def double(self, i, j):
        b = self.bond(i, j)
        return b is not None and b.GetBondType() == Chem.BondType.DOUBLE and not self.aromatic_bond(i, j)

    def bond_order(self, i, j):
        if self.aromatic_bond(i, j):
            return None
        b = self.bond(i, j)
        t = b.GetBondType()
        return {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}.get(t)

    def h(self, i):
        return self.hs[i] + sum(1 for j in self.nb[i] if self.el[j] == 'H')

    def heavy(self, i):
        return [j for j in self.nb[i] if self.el[j] != 'H']

    def side(self, i, j, cut=()):
        """The atoms reached from i without crossing bond i-j (and the bonds in cut); None when j is reached."""
        key = (i, j, tuple(sorted(cut)))
        if key in self._sides:
            return self._sides[key]
        blocked = {frozenset((i, j))} | {frozenset(c) for c in cut}
        seen, todo = {i}, [i]
        while todo:
            x = todo.pop()
            for y in self.nb[x]:
                if y not in seen and frozenset((x, y)) not in blocked:
                    seen.add(y); todo.append(y)
        out = None if j in seen else frozenset(seen)
        self._sides[key] = out
        return out

    def cut_orders(self, i, S):
        out = 0
        for j in self.nb[i]:
            if j not in S:
                o = self.bond_order(i, j)
                if o is None:
                    return None                             # an aromatic bond: no capped piece
                out += o
        return out

    def hcap(self, S):
        total = 0
        for i in S:
            c = self.cut_orders(i, S)
            if c is None:
                return None
            total += self.hs[i] + c
        return total

    def counts(self, S, h):
        out = {}
        for i in S:
            out[self.key[i]] = out.get(self.key[i], 0) + 1
        if h:
            out['H'] = out.get('H', 0) + h
        return out

    # ------------------------------------------------------------------------------------------- atom types
    def carbonyl_o(self, c):
        """The O (or S) double-bonded to carbon c, else None."""
        if self.el[c] != 'C' or self.aro[c]:
            return None
        return next((j for j in self.nb[c] if self.el[j] in ('O', 'S') and self.double(c, j)), None)

    def saturated(self, c):
        return self.el[c] == 'C' and not self.aro[c] and all(self.single(c, j, ring=True) for j in self.nb[c])

    def vinyl(self, c):
        return any(self.el[j] == 'C' and self.double(c, j) and not self.aro[j] for j in self.nb[c])

    def anomeric(self, c, x):
        """c is a sugar's anomeric carbon (in a ring, bonded to the ring O) and x an atom outside that ring."""
        if self.el[c] != 'C' or self.aro[c]:
            return False
        ring = self.mol.GetRingInfo()
        for o in self.nb[c]:
            if o != x and self.el[o] == 'O' and self.bond(c, o).IsInRing():
                for r in ring.AtomRings():
                    if c in r and o in r and x not in r and len(r) in (5, 6):
                        return True
        return False

    def cation(self, c, S):
        """0..1: how well carbon c carries a positive charge in piece S."""
        if self.el[c] != 'C':
            return 0.3
        if self.aro[c]:
            return 0.15
        if self.carbonyl_o(c) is not None:
            return 0.85
        if not self.saturated(c):
            return 0.1                                      # vinyl cation
        near = [j for j in self.nb[c] if j in S and self.el[j] != 'H']
        if any(self.el[j] in ('N', 'O', 'S') and self.single(c, j, ring=True) and self.carbonyl_o(j) is None
               for j in near):
            return 0.7                                      # iminium / oxocarbenium
        if any(self.aro[j] for j in near):
            return 0.75                                     # benzylic
        if any(self.vinyl(j) for j in near) or self.vinyl(c):
            return 0.55                                     # allylic
        return {0: 0.05, 1: 0.15, 2: 0.35}.get(len(near), 0.6)


class _Engine:
    def __init__(self, m, mode):
        self.m = m
        self.mode = mode
        allh = sum(m.hs) + sum(1 for e in m.el if e == 'H')
        self.M = m.counts(range(m.n), sum(m.hs))
        drawn_z = sum(m.q)
        drawn_rad = sum(m.rad)
        self.charged = bool(drawn_z)
        if self.charged:                                    # an ion drawn as such: used as it is
            self.z, self.odd, self.pre = drawn_z, bool(drawn_rad % 2), dict(self.M)
            self.label = 'M' + ('•' if self.odd else '') + chem.charge_text(drawn_z)
            self.mode = 'ESI+' if drawn_z > 0 else 'ESI-'
            if self.odd:
                self.mode = 'EI'
        elif mode == 'EI':
            self.z, self.odd, self.pre, self.label = 1, not drawn_rad % 2, dict(self.M), '[M]•+'
        elif mode == 'ESI+':
            self.z, self.odd, self.label = 1, bool(drawn_rad % 2), '[M+H]+'
            self.pre = dict(self.M, H=self.M.get('H', 0) + 1)
        else:
            self.z, self.odd, self.label = -1, bool(drawn_rad % 2), '[M−H]−'
            self.pre = dict(self.M, H=self.M.get('H', 0) - 1)
        self.allh = allh
        self.out = []

    # ------------------------------------------------------------------------------------------- making fragments
    def mz(self, counts, z):
        mass = sum(chem.key_mass(k) * n for k, n in counts.items())
        return (mass - z * chem.ELECTRON) / abs(z)

    def place(self, S, base, site, dq, rad, dh):
        """Marks: base ones plus dq / rad / dh at site; H that the site cannot give comes from its neighbours."""
        m = self.m
        marks = {k: list(v) for k, v in (base or {}).items() if k in S}
        marks.setdefault(site, [0, 0, 0])
        marks[site][0] += dq

        def room(i):
            c = m.cut_orders(i, S)
            return None if c is None else m.hs[i] + c + marks.get(i, [0, 0, 0])[2]
        if rad and dh < 0:                                  # distonic: the unpaired electron on a neighbour, −H
            j = next((j for j in m.nb[site] if j in S and j not in marks and m.el[j] == 'C' and room(j)), None)
            if j is not None:
                marks[j] = [0, rad, -1]
                rad, dh = 0, dh + 1
        marks[site][1] += rad
        r = room(site)
        if r is None:
            return None
        if dh >= 0 or r + dh >= 0:
            marks[site][2] += dh
        else:
            marks[site][2] -= r
            need = -(r + dh)
            seen, todo = {site}, [site]
            while need and todo:
                x = todo.pop(0)
                for y in m.nb[x]:
                    if y in S and y not in seen:
                        seen.add(y); todo.append(y)
                        ry = room(y)
                        if ry and not m.aro[y]:             # never an H off a benzene ring (no arynes)
                            take = min(ry, need, 1)
                            marks.setdefault(y, [0, 0, 0])[2] -= take
                            need -= take
                            if not need:
                                break
            if need:
                return None
        for k, v in marks.items():                          # a ring carbon losing an H with nothing to show for it
            if m.aro[k] and v[2] < 0 and not v[0] and not v[1]:     # would make an aryne
                return None
        return {k: tuple(v) for k, v in marks.items() if any(v)}

    def acid_site(self, S, avoid=()):
        """The O–H in S most likely to lose its proton in negative mode: a carboxylic acid, then a phenol, then an
        alcohol (None when there is none)."""
        m = self.m
        best = None
        for o in S:
            if o in avoid or m.el[o] != 'O' or not m.hs[o] or len(m.heavy(o)) != 1:
                continue
            c = m.heavy(o)[0]
            rank = 0 if m.carbonyl_o(c) is not None else 1 if m.aro[c] else 2
            if best is None or rank < best[0]:
                best = (rank, o)
        return best[1] if best else None

    def base_site(self, S, avoid=()):
        """The atom in S most likely to hold the proton in positive mode: an amine N, a pyridine-type N, an aniline
        N, a C=O oxygen, another O (None when there is none)."""
        m = self.m
        best = None
        for i in S:
            if i in avoid or m.q[i] or m.el[i] not in ('N', 'O'):
                continue
            hv = m.heavy(i)
            if m.el[i] == 'N':
                if m.aro[i]:
                    rank = 1 if len(hv) == 2 and not m.hs[i] else 9
                elif any(m.carbonyl_o(j) is not None for j in hv) or not all(m.single(i, j, True) for j in m.nb[i]):
                    rank = 9
                else:
                    rank = 2 if any(m.aro[j] for j in hv) else 0
            elif any(m.double(i, j) for j in m.nb[i]):
                rank = 3
            else:
                rank = 4
            if rank < 9 and (best is None or rank < best[0]):
                best = (rank, i)
        return best[1] if best else None

    def make(self, S, marks, odd, how, score, step=1, lost_name=None, parent=None, mech=None):
        """A fragment from atom set S (molecule indices) and marks, or None when the bookkeeping fails."""
        m = self.m
        S = frozenset(S)
        if not S or len(S) == m.n and not marks:
            return None
        hc = m.hcap(S)
        if hc is None or marks is None:
            return None
        h = hc + sum(v[2] for v in marks.values())
        for i, v in marks.items():
            if i not in S or m.hs[i] + m.cut_orders(i, S) + v[2] < 0:
                return None
        z = sum(m.q[i] for i in S) + sum(v[0] for v in marks.values())
        if z != self.z or h < 0:
            return None
        counts = _clean(m.counts(S, h))
        if parent is not None and (not counts.get('C') or sum(v for k, v in counts.items() if k != 'H') < 2):
            return None                                     # nothing sensible is left
        start = parent.counts if parent else self.pre
        lost = _sub(start, counts)
        if any(v < 0 for v in lost.values()) or not any(lost.values()):
            return None
        lost = _clean(lost)
        was_odd = parent.odd if parent else self.odd
        name = lost_name or (chem.formula_text(lost) + ('•' if was_odd != odd else ''))
        total = name if parent is None else f'{parent.lost} + {name}'
        return self._frag(counts, odd, total, how, score, step, S, marks, mech, parent)

    def formula_only(self, lost, radical, how, score, highlight, parent=None, step=1, mech=None):
        start = parent.counts if parent else self.pre
        counts = _sub(start, lost)
        if any(v < 0 for v in counts.values()):
            return None
        counts = _clean(counts)
        if not counts.get('C') and not counts.get('N'):
            return None
        odd = (parent.odd if parent else self.odd) != radical
        name = chem.formula_text(lost) + ('•' if radical else '')
        total = name if parent is None else f'{parent.lost} + {name}'
        f = self._frag(counts, odd, total, how, score, step, None, None, mech, parent)
        f.atoms = frozenset(highlight)
        return f

    def fixed_ion(self, counts, odd, how, score, highlight, mech=None):
        """A diagnostic ion known from the literature (Si(CH3)3+, phosphocholine ...): its formula only."""
        lost = _sub(self.pre, counts)
        if any(v < 0 for v in lost.values()) or not any(lost.values()):
            return None
        name = chem.formula_text(_clean(lost)) + ('•' if odd != self.odd else '')
        f = self._frag(dict(counts), odd, name, how, score, 1, None, None, mech, None)
        f.atoms = frozenset(highlight)
        return f

    def _frag(self, counts, odd, lost, how, score, step, S, marks, mech=None, parent=None):
        m = self.m
        f = Fragment(self.mz(counts, self.z), self.z, odd, counts, ion_text(counts, self.z, odd), lost, how,
                     max(0.0, min(1.0, score)), step=step)
        f.parent = parent
        if mech:
            text, arrows, pre, ref = (tuple(mech) + ((), (), ''))[:4]
            f.mech, f.ref = text, REFS.get(ref, ref)

            def end(e):
                kind, v = e
                if kind == 'bond':
                    return ('bond', (m.order[v[0]], m.order[v[1]]))
                return (kind, m.order[v])
            f.arrows = [(k, end(a), end(b)) for k, a, b in arrows]
            f.pre = [(m.order[a], dq, dr, dh) for a, dq, dr, dh in pre]
        if S is not None:
            f.S = frozenset(m.order[i] for i in S)
            f.marks = {m.order[i]: v for i, v in marks.items()}
            f._idx = S
            f._imarks = marks
            f.atoms = frozenset(m.owner[m.order[i]] for i in S)
        else:
            f._idx = None
        return f

    def add(self, f):
        if f is not None:
            self.out.append(f)
        return f

    def proton(self, x):
        """Where the charge sits before an ESI step at atom x: the proton on x ([M+H]+), or the deprotonated acidic
        OH ([M−H]−), as [(atom, charge, •, H)]."""
        if self.z > 0:
            return [(x, 1, 0, 1)] if x is not None else []
        acid = self.acid_site(frozenset(range(self.m.n)), ())
        return [(acid, -1, 0, -1)] if acid is not None else []

    def beta_h(self, c, S, avoid=()):
        """A neighbour of c in S carrying an H (for eliminations), or None."""
        m = self.m
        return next((j for j in m.heavy(c) if j in S and j not in avoid and m.h(j)), None)

    # ------------------------------------------------------------------------------------------- the rules
    def run(self, steps=2):
        m = self.m
        if self.charged:
            self.charge_remote()
        elif self.mode == 'EI':
            self.ei_cleavages()
            self.mclafferty()
        else:
            self.esi_cleavages()
        if not self.charged:
            self.rda()
        self.group_losses()
        frontier = [f for f in self.out if f._idx is not None]
        for step in range(2, steps + 1):
            frontier = sorted(frontier, key=lambda f: -f.score)[:60]
            new = []
            for f in frontier:
                new += self.next_losses(f, step)
            self.out += new
            frontier = [f for f in new if f._idx is not None]
        return self.merged()

    def ei_cleavages(self):
        m = self.m
        alpha = {'N': 0.85, 'S': 0.7, 'O': 0.65}
        for y in range(m.n):                                # α-cleavage next to a heteroatom
            if m.el[y] not in alpha or m.aro[y] or m.q[y] > 0:
                continue
            acyl = any(m.carbonyl_o(j) is not None for j in m.nb[y])
            kind = {'N': 'iminium', 'O': 'oxonium', 'S': 'sulfonium'}[m.el[y]]
            Y = m.el[y]
            for ca in m.nb[y]:
                if not m.single(ca, y, ring=True) or not m.saturated(ca):
                    continue
                if m.h(ca):                                 # [M−H]+: the α-hydrogen leaves as H•
                    self.add(self.make(range(m.n), {ca: (1, 0, -1)}, False, f'loss of H• α to {Y} ({kind} ion)',
                                       {'N': 0.45}.get(Y, 0.3) - (0.2 if acyl else 0), mech=(
                        f'Ionisation at {Y} ({Y}•+). One electron of an α C–H bond pairs with the unpaired electron '
                        f'on {Y} (C={Y}+ forms), the other leaves with the hydrogen atom: [M−H]+.',
                        [('half', ('bondH', ca), ('bond', (y, ca))), ('half', ('atom', y), ('bond', (y, ca))),
                         ('half', ('bondH', ca), ('H', ca))], [(y, 1, 1, 0)], 'ei')))
                for r in m.heavy(ca):
                    if r == y or not m.single(ca, r):
                        continue
                    S = m.side(ca, r)
                    if S is None:
                        continue
                    score = alpha[m.el[y]] - (0.3 if acyl else 0) + (0.05 if len(m.side(r, ca)) >= 2 else 0)
                    self.add(self.make(S, {ca: (1, 0, -1)}, False, f'α-cleavage at {Y} ({kind} ion)', score, mech=(
                        f'Ionisation removes an electron from a lone pair of {Y} ({Y}•+). One electron of the '
                        f'Cα–{m.el[r]} bond pairs with the unpaired electron on {Y}, making the C={Y}+ double bond of '
                        f'the {kind} ion; the other electron leaves with that side as a radical. The largest radical '
                        "is lost most easily (Stevenson's rule).",
                        [('half', ('bond', (ca, r)), ('bond', (y, ca))), ('half', ('atom', y), ('bond', (y, ca))),
                         ('half', ('bond', (ca, r)), ('atom', r))], [(y, 1, 1, 0)], 'ei')))
        for c in range(m.n):                                # α-cleavage at C=O: acylium ions
            o = m.carbonyl_o(c)
            if o is None:
                continue
            if m.hs[c]:                                     # an aldehyde loses its H•
                self.add(self.make(range(m.n), {c: (1, 0, -1)}, False, 'loss of H• from the aldehyde (acylium ion)',
                                   0.55, mech=(
                    'Ionisation at the C=O oxygen. The aldehyde C–H bond gives one electron to the oxygen\'s unpaired '
                    'electron (C≡O+ forms) and the hydrogen atom leaves: [M−H]+, strong for aldehydes.',
                    [('half', ('bondH', c), ('bond', (c, o))), ('half', ('atom', o), ('bond', (c, o))),
                     ('half', ('bondH', c), ('H', c))], [(o, 1, 1, 0)], 'ei')))
            for r in m.heavy(c):
                if r == o or not m.single(c, r):
                    continue
                S = m.side(c, r)
                if S is None:
                    continue
                if m.el[r] == 'C':
                    score = 0.85
                elif m.el[r] == 'O':
                    score = 0.45 if m.h(r) else 0.75
                elif m.el[r] == 'N':
                    score = 0.55
                else:
                    score = 0.4
                self.add(self.make(S, {c: (1, 0, -1)}, False, 'α-cleavage at C=O (acylium ion)', score, mech=(
                    f'Ionisation removes a lone-pair electron of the C=O oxygen (O•+). One electron of the '
                    f'C(=O)–{m.el[r]} bond pairs with it, making the C≡O+ triple bond of the acylium ion; the other '
                    'electron leaves with the group as a radical.',
                    [('half', ('bond', (c, r)), ('bond', (c, o))), ('half', ('atom', o), ('bond', (c, o))),
                     ('half', ('bond', (c, r)), ('atom', r))], [(o, 1, 1, 0)], 'ei')))
                if m.el[r] == 'C':
                    S2 = m.side(r, c)
                    self.add(self.make(S2, {r: (1, 0, -1)}, False, 'α-cleavage at C=O (charge on the other side)',
                                       0.6 * m.cation(r, S2), mech=(
                        'Inductive cleavage: the charged C=O oxygen pulls the C–C bond pair toward the carbonyl '
                        'carbon; the acyl radical leaves and the alkyl / aryl cation remains.',
                        [('full', ('bond', (c, r)), ('atom', c))], [(o, 1, 1, 0)], 'ei')))
        for b in range(m.n):                                # benzylic / allylic cleavage
            if not m.saturated(b):
                continue
            aro = [j for j in m.nb[b] if m.aro[j]]
            allyl = [j for j in m.nb[b] if m.vinyl(j)]
            if not aro and not allyl:
                continue
            for r in m.heavy(b):
                if r in aro or r in allyl or not m.single(b, r):
                    continue
                S = m.side(b, r)
                if S is None:
                    continue
                pi = (aro or allyl)[0]
                f = self.make(S, {b: (1, 0, -1)}, False, 'benzylic cleavage' if aro else 'allylic cleavage',
                              0.9 if aro else 0.6, mech=(
                    ('Ionisation of the aromatic π system (ring•+). One electron of the benzylic C–R bond pairs with '
                     'the ring\'s unpaired electron, the other leaves with R•. The benzyl cation formed rearranges '
                     'to the seven-membered aromatic tropylium ion.') if aro else
                    ('Ionisation of the C=C π bond. One electron of the allylic C–R bond pairs with the π radical, the '
                     'other leaves with R•: the resonance-stabilised allyl cation remains.'),
                    [('half', ('bond', (b, r)), ('bond', (b, pi))), ('half', ('bond', (b, r)), ('atom', r))],
                    [(pi, 1, 1, 0)], 'ei'))
                if f is not None and f.counts == {'C': 7, 'H': 7}:
                    f.how = 'benzylic cleavage (tropylium ion C7H7+)'
                self.add(f)
        base = {'I': 0.75, 'Br': 0.7, 'Cl': 0.45, 'F': 0.15}
        for x in range(m.n):                                # inductive cleavage at C–X
            for c in m.nb[x]:
                if m.el[c] != 'C' or not m.single(c, x):
                    continue
                S = m.side(c, x)
                if S is None:
                    continue
                X = m.el[x]
                mech = (f'Ionisation at {X} ({X}•+). The C–{X} bond pair moves to {X} (inductive cleavage, i): '
                        f'{X} leaves as a radical and the carbon keeps the positive charge.',
                        [('full', ('bond', (c, x)), ('atom', x))], [(x, 1, 1, 0)], 'ei')
                if X in base:
                    s = base[X] * (0.5 + 0.5 * m.cation(c, S)) * (0.6 if m.aro[c] else 1)
                    self.add(self.make(S, {c: (1, 0, -1)}, False, f'inductive cleavage at C–{X}', s, mech=mech))
                elif X in ('O', 'S') and not m.h(x) and m.saturated(c):
                    s = 0.3 * (0.5 + m.cation(c, S))
                    self.add(self.make(S, {c: (1, 0, -1)}, False, f'inductive cleavage at C–{X}', s, mech=mech))
        for a in range(m.n):                                # σ-cleavage of C–C bonds (branch points)
            if not m.saturated(a):
                continue
            for b in m.nb[a]:
                if b < a or not m.saturated(b) or not m.single(a, b):
                    continue
                for c, other in ((a, b), (b, a)):
                    S = m.side(c, other)
                    if S is None or len(S) < 2:
                        continue
                    s = m.cation(c, S)
                    if s >= 0.35:
                        self.add(self.make(S, {c: (1, 0, -1)}, False, 'σ-cleavage of C–C (alkyl cation)',
                                           0.1 + 0.4 * s, mech=(
                            'σ-Ionisation of the C–C bond (no lone pair or π bond nearby). The one-electron bond '
                            'breaks; the more substituted carbon keeps the charge (the most stable carbocation) and '
                            'the other side leaves as a radical - branch points break most.',
                            [('half', ('bond', (c, other)), ('atom', other))], [(c, 1, 1, 0)], 'ei')))
        for x in range(m.n):                                # alkene loss with H transfer from ethers / esters
            if m.el[x] not in ('O', 'N') or m.aro[x] or m.q[x]:
                continue
            for c in m.heavy(x):
                a = next((j for j in m.heavy(x) if j != c and (m.aro[j] or m.carbonyl_o(j) is not None)), None)
                if a is None or not m.saturated(c) or not m.single(c, x):
                    continue
                Sx, Sc = m.side(x, c), m.side(c, x)
                if Sx is None:
                    continue
                bh = self.beta_h(c, Sc, (x,))
                if bh is None:
                    continue
                ester = m.carbonyl_o(a) is not None
                arrows = [('half', ('bondH', bh), ('atom', x)), ('half', ('bondH', bh), ('bond', (bh, c))),
                          ('half', ('bond', (c, x)), ('bond', (bh, c)))]
                o_ion = m.carbonyl_o(a) if ester else x
                self.add(self.make(Sx, {x: (1, 1, 0)}, True, 'alkene loss, H to ' + ('the ester O' if ester else
                                                                                     m.el[x] + ' (four-centre)'),
                                   0.4 if ester else 0.6, mech=(
                    'A β-hydrogen of the alkyl group moves to the heteroatom as the C–' + m.el[x] + ' bond breaks '
                    '(four-membered transition state); the alkene leaves and the ' +
                    ('acid' if ester else 'phenol / amine') + ' radical cation remains (phenetole → phenol•+, m/z 94).',
                    arrows, [(o_ion, 1, 1, 0)], 'ei')))
                if ester and len([j for j in Sc if m.el[j] == 'C']) >= 2:
                    self.add(self.make(Sx, {x: (1, 0, 1)}, False, 'double hydrogen rearrangement (McLafferty+1)',
                                       0.5, mech=(
                        'Two hydrogens of the alkoxy group move to the ester oxygens (two consecutive H transfers) '
                        'and an allyl / alkenyl radical leaves: the protonated acid [RC(OH)2]+ (ethyl acetate → '
                        'm/z 61).', arrows, [(m.carbonyl_o(a), 1, 1, 0)], 'mclafferty')))

    def mclafferty(self):
        m = self.m
        done = set()
        for c1 in range(m.n):
            ys = [y for y in m.nb[c1] if m.double(c1, y) and m.el[y] in ('O', 'N', 'S', 'C') and not m.aro[y]]
            ortho = [y for y in m.nb[c1] if m.aro[y] and m.el[y] == 'C'] if m.aro[c1] else []
            if m.el[c1] != 'C' or (not ys and not ortho):
                continue
            for xa in m.heavy(c1):
                if xa in ys or m.aro[xa] or m.el[xa] not in ('C', 'O', 'N') or not m.single(c1, xa, ring=True):
                    continue
                if m.el[xa] == 'C' and not m.saturated(xa):
                    continue
                for cb in m.heavy(xa):
                    if cb == c1 or not m.saturated(cb) or not m.single(xa, cb):
                        continue
                    S = m.side(xa, cb)
                    if S is None:
                        continue
                    gs = [g for g in m.nb[cb] if g not in S and m.el[g] == 'C' and m.h(g) and not m.aro[g]]
                    if not gs:
                        continue
                    for y in ys:
                        if (y, xa, cb) in done:
                            continue
                        done.add((y, xa, cb))
                        s = {'O': 0.85, 'S': 0.6, 'N': 0.6, 'C': 0.35}[m.el[y]] - (0.15 if m.el[xa] == 'O' else 0)
                        g = gs[0]
                        arrows = [('half', ('atom', y), ('H', g)), ('half', ('bondH', g), ('bond', (cb, g))),
                                  ('half', ('bond', (xa, cb)), ('bond', (cb, g))),
                                  ('half', ('bond', (xa, cb)), ('bond', (c1, xa)))]
                        Y = m.el[y]
                        self.add(self.make(S, {y: (1, 1, 0)}, True, 'McLafferty rearrangement', s, mech=(
                            f'Ionisation at {Y} ({Y}•+). Through a six-membered ring the unpaired electron on {Y} '
                            'takes a γ-hydrogen; the Cα–Cβ bond then breaks: one of its electrons pairs with the '
                            'one left on Cγ (the neutral alkene Cβ=Cγ leaves), the other makes the C=C of the enol '
                            'radical cation.', arrows, [(y, 1, 1, 0)], 'mclafferty')))
                        S2 = m.side(cb, xa)
                        self.add(self.make(S2, {cb: (1, 0, -1), g: (0, 1, -1)}, True,
                                           'McLafferty rearrangement (charge on the alkene)', 0.3, mech=(
                            'The same γ-hydrogen transfer and Cα–Cβ cleavage, but the charge stays on the alkene '
                            'side (when its ionisation energy is the lower).', arrows, [(y, 1, 1, 0)],
                            'mclafferty')))
                    if ortho and not ys and m.el[xa] == 'C':
                        g = gs[0]
                        self.add(self.make(S, {ortho[0]: (1, 1, 0)}, True, 'McLafferty-type H shift to the ring',
                                           0.55, mech=(
                            'The ring radical cation takes a γ-hydrogen at its ortho carbon (six-membered '
                            'transition state) and the benzylic Cα–Cβ bond breaks: an alkene leaves and the '
                            'methylene-cyclohexadiene radical cation remains (propylbenzene → m/z 92).',
                            [('half', ('atom', ortho[0]), ('H', g)), ('half', ('bondH', g), ('bond', (cb, g))),
                             ('half', ('bond', (xa, cb)), ('bond', (cb, g)))], [(ortho[0], 1, 1, 0)],
                            'mclafferty')))

    def rda(self):
        """Retro-Diels–Alder: a six-membered ring with a double (or aromatic fusion) bond a=b splits at the two
        bonds across from it into the diene and the dienophile."""
        m = self.m
        seen = set()
        for ring in m.mol.GetRingInfo().AtomRings():
            if len(ring) != 6:
                continue
            for k in range(6):
                p, a, b, q, u, w = (ring[(k + d) % 6] for d in (-1, 0, 1, 2, 3, 4))
                if not (m.double(a, b) or m.aromatic_bond(a, b)):
                    continue
                if not (m.single(p, a, ring=True) and m.single(b, q, ring=True) and m.single(q, u, ring=True)
                        and m.single(p, w, ring=True)):
                    continue
                if m.aromatic_bond(u, w):
                    continue
                key = frozenset((frozenset((q, u)), frozenset((p, w))))
                if key in seen:
                    continue
                seen.add(key)
                D = m.side(a, -1, cut=((q, u), (p, w)))
                if D is None or u in D or w in D:
                    continue
                E = frozenset(range(m.n)) - D
                arrows = [('full', ('bond', (a, b)), ('bond', (a, p))), ('full', ('bond', (p, w)), ('bond', (w, u))),
                          ('full', ('bond', (q, u)), ('bond', (b, q)))]
                text = ('The six-membered ring opens at the two bonds across from its double (or fused aromatic) bond '
                        '- a retro-Diels–Alder reaction: the three electron pairs of the ring shift together, giving a '
                        'diene and an alkene (the dienophile). ')
                if self.mode == 'EI':
                    site = p if m.el[p] == 'C' else q
                    other = q if site == p else p
                    mech = (text + 'In the radical cation the charge usually stays on the diene.', arrows,
                            [(a, 1, 1, 0)], 'ei')
                    self.add(self.make(D, {site: (1, 0, -1), other: (0, 1, -1)}, True, 'retro-Diels–Alder (diene)',
                                       0.7, mech=mech))
                    self.add(self.make(E, {w: (1, 0, -1), u: (0, 1, -1)}, True, 'retro-Diels–Alder (dienophile)',
                                       0.3, mech=mech))
                else:
                    sign = 1 if self.z > 0 else -1
                    for piece, ends, s, part in ((D, (p, q), 0.65, '1,3A-type diene'),
                                                 (E, (w, u), 0.35, 'dienophile')):
                        if sign > 0:            # the proton on the piece: one end keeps its cap H, the other + and −H
                            c = sorted(ends, key=lambda i: (m.el[i] != 'C', m.carbonyl_o(i) is None))[0]
                            marks = self.place(piece, {}, c, 1, 0, -1)
                        else:                   # an acidic OH gives the proton; the two ends close the diene
                            acid = self.acid_site(piece, ends)
                            if acid is not None:
                                marks = self.place(piece, {}, acid, -1, 0, -1)
                                for e in ends:
                                    marks = self.place(piece, marks, e, 0, 0, -1) if marks else None
                            else:
                                marks = None            # nothing to deprotonate: no [M−H]− ion of this piece
                        if sign > 0:                    # the proton on the most basic atom of the molecule
                            pre = self.proton(self.base_site(frozenset(range(m.n))))
                        else:
                            pre = self.proton(None)
                        self.add(self.make(piece, marks, False, f'retro-Diels–Alder ({part})', s, mech=(
                            text + ('The protonated (or deprotonated) molecule splits the same way; for flavonoids '
                                    'the A-ring piece is the 1,3A ion and the B-ring piece the 1,3B ion.'),
                            arrows, pre, 'rda')))

    def esi_cleavages(self):
        """Heterolysis of C–O, C–N, C–S bonds outside rings, as in CID of [M+H]+ / [M−H]−."""
        m = self.m
        pos = self.z > 0
        for c in range(m.n):
            if m.el[c] != 'C':
                continue
            for x in m.nb[c]:
                if m.el[x] not in ('N', 'O', 'S') or not m.single(c, x):
                    continue
                Sc, Sx = m.side(c, x), m.side(x, c)
                if Sc is None:
                    continue
                acyl = m.carbonyl_o(c) is not None
                glyco = m.anomeric(c, x) and m.el[x] in ('O', 'N', 'S')
                pendant = len(m.heavy(x)) == 1                  # OH, NH2, SH: their losses are group losses
                x_acyl = any(m.carbonyl_o(j) is not None for j in m.heavy(x) if j != c)
                x_aryl = any(m.aro[j] for j in m.heavy(x) if j != c)
                beta_h = any(m.h(j) for j in m.heavy(c) if j in Sc and j != x) and not m.aro[c]
                if pos and not pendant:                         # the carbon side as the cation, XH leaves
                    if glyco:
                        s, how = 0.55, 'glycosidic cleavage (sugar oxonium ion, B)'
                    elif acyl:
                        s, how = (0.7, 'amide cleavage (acylium, b-type)') if m.el[x] == 'N' else \
                            (0.8, 'ester cleavage (acylium ion)')
                    else:
                        s = 0.9 * m.cation(c, Sc)
                        how = f'C–{m.el[x]} cleavage (carbocation)'
                    X = m.el[x]
                    if s >= 0.15:
                        self.add(self.make(Sc, {c: (1, 0, -1)}, False, how, s, mech=(
                            f'Protonation of {X} turns it into a good leaving group. The C–{X} bond pair goes to '
                            f'{X}, which leaves as a neutral molecule, and the ' +
                            ('acylium (b) ion' if acyl else 'oxocarbenium (B) ion' if glyco else 'carbocation') +
                            ' remains (charge-migration fragmentation).',
                            [('full', ('bond', (c, x)), ('atom', x))], [(x, 1, 0, 1)],
                            'peptide' if acyl and X == 'N' else 'glyco' if glyco else 'esi')))
                if not pos and glyco and (x_aryl or any(m.vinyl(j) for j in m.heavy(x) if j != c)):
                    acid = self.acid_site(Sx, (x,))                 # flavonol 3-O-glycosides: radical aglycone
                    if acid is not None:
                        marks = {acid: (-1, 0, -1), x: (0, 1, -1)}
                        flavonol = any(m.vinyl(j) for j in m.heavy(x) if j != c)
                        self.add(self.make(Sx, marks, True, 'homolytic glycosidic cleavage (radical aglycone '
                                           '[Y0−H]•−)', 0.75 if flavonol else 0.4, mech=(
                            'The glycosidic O–C bond breaks homolytically: the aglycone keeps the oxygen as a '
                            'radical next to the deprotonated phenol ([Y0−H]•−) and the sugar leaves as a radical. '
                            'Strong for flavonol 3-O-glycosides, where the radical is conjugated with the 2,3-double '
                            'bond and C-4 carbonyl.', [('half', ('bond', (c, x)), ('atom', x)),
                                                        ('half', ('bond', (c, x)), ('atom', c))],
                            self.proton(None), 'radical')))
                if not pos and acyl and m.el[x] == 'O' and not pendant:   # the acyl side as the anion (ketene)
                    acid = self.acid_site(Sc)
                    alpha = self.beta_h(c, Sc, (x,))
                    if acid is not None and alpha is not None:
                        marks = self.place(Sc, {acid: (-1, 0, -1)}, c, 0, 0, -1)
                        marks = self.place(Sc, marks, alpha, 0, 0, -1) if marks else None
                        self.add(self.make(Sc, marks, False, 'ester cleavage (acyl side as the anion, ketene)', 0.45,
                                           mech=(
                            'The deprotonated phenol / acid on the acyl side stays charged while the ester O–C(=O) '
                            'bond breaks with an α-hydrogen moving to the alkoxy oxygen: the alcohol leaves as a '
                            'neutral molecule and the acyl part remains as a ketene anion (chlorogenic acids: '
                            'caffeoyl m/z 161).',
                            [('full', ('bondH', alpha), ('bond', (alpha, c))), ('full', ('bond', (c, x)), ('atom', x))],
                            self.proton(None), 'chlorogenic')))
                if not beta_h or pendant and not glyco:
                    continue
                bh = self.beta_h(c, Sc, (x,))
                elim = [('full', ('bondH', bh), ('bond', (bh, c))), ('full', ('bond', (c, x)), ('atom', x))]
                if pos:                                         # the heteroatom side takes the proton
                    if glyco:
                        s, how = 0.9, 'glycosidic cleavage (aglycone or sugar loss, Y ion)'
                    elif acyl and m.el[x] == 'N':
                        s, how = 0.65, 'amide cleavage, charge on the amine (y-type)'
                    elif acyl:
                        s, how = 0.35, 'ester cleavage, charge on the alcohol (ketene lost)'
                    elif x_acyl and m.el[x] == 'O':
                        s, how = 0.5, 'ester cleavage: [acid+H]+, alkene lost'
                    elif x_acyl:
                        s, how = 0.35, 'C–N cleavage: [amide+H]+, alkene lost'
                    elif m.el[x] == 'N':
                        s, how = 0.45, 'C–N cleavage: amine, alkene lost'
                    else:
                        s, how = 0.3, f'C–{m.el[x]} cleavage, alkene lost'
                    self.add(self.make(Sx, {x: (1, 0, 1)}, False, how, s, mech=(
                        f'Proton on {m.el[x]}. As the C–{m.el[x]} bond breaks, a hydrogen next to the carbon '
                        f'goes with the electrons into a C=C (or C=C=O) bond and is handed over to {m.el[x]} '
                        '(1,2-elimination, charge retention on the heteroatom side): the neutral ' +
                        ('anhydro sugar leaves and the aglycone (Y) ion remains.' if glyco else
                         'ketene leaves.' if acyl else 'alkene leaves.'),
                        elim, [(x, 1, 0, 1)], 'glyco' if glyco else 'peptide' if acyl and m.el[x] == 'N' else 'esi')))
                else:                                           # the heteroatom side as the anion
                    if glyco:
                        s, how = 0.85, 'glycosidic cleavage (aglycone or sugar loss, Y⁻ ion)'
                    elif x_acyl and m.el[x] == 'O':
                        s, how = 0.7, 'ester cleavage (carboxylate ion)'
                    elif acyl and m.el[x] == 'O':
                        s, how = 0.6, 'ester cleavage (alkoxide ion, ketene lost)'
                    elif x_aryl:
                        s, how = 0.5, 'C–O cleavage (phenoxide ion)'
                    elif acyl:
                        s, how = 0.3, 'amide cleavage (amide anion)'
                    else:
                        s, how = 0.2, f'C–{m.el[x]} cleavage (anion)'
                    self.add(self.make(Sx, {x: (-1, 0, -1)}, False, how, s, mech=(
                        f'The negative charge (the deprotonated acidic OH) takes a hydrogen next to the carbon as the '
                        f'C–{m.el[x]} bond pair goes to {m.el[x]}: {m.el[x]} leaves with the charge and the carbon side '
                        'becomes a neutral ' + ('anhydro sugar.' if glyco else 'ketene.' if acyl else 'alkene.'),
                        elim, self.proton(None), 'chlorogenic' if x_acyl or acyl else 'glyco' if glyco else 'esi')))

    def charge_remote(self):
        """A structure drawn as an ion: each single bond outside a ring splits, the charged side keeping an H from
        the side that leaves as a neutral molecule."""
        m = self.m
        for i in range(m.n):
            for j in m.heavy(i):
                if j < i or not m.single(i, j):
                    continue
                for a, b in ((i, j), (j, i)):
                    S = m.side(a, b)
                    if S is None or not any(m.q[k] for k in S):
                        continue
                    if not any(m.h(k) for k in m.heavy(b) if k not in S) and not m.h(b):
                        continue
                    self.add(self.make(S, {}, self.odd, 'cleavage, neutral molecule lost (charge-remote)', 0.45,
                                       mech=('The charge stays where it is drawn: a hydrogen moves to the charged side '
                                             'as the bond breaks and the other side leaves as a neutral molecule '
                                             '(charge-remote fragmentation).', [], [], 'esi')))

    def _matches(self):
        if not hasattr(self, '_gm'):
            self._gm = []
            for row in GROUP_LOSSES:
                patt = Chem.MolFromSmarts(row[1])
                self._gm.append((row, self.m.mol.GetSubstructMatches(patt) if patt is not None else ()))
        return self._gm

    def group_losses(self):
        m = self.m
        everything = frozenset(range(m.n))
        for row, matches in self._matches():
            name, _s, rem, att, _lh, radical, scores, what = row
            s = scores.get(self.mode)
            if not s:
                continue
            for match in matches:
                gone = {match[k] for k in rem}
                gone |= {j for g in gone for j in m.nb[g] if m.el[j] == 'H'}
                S = everything - gone
                attach = match[att]
                if attach not in S or m.side(attach, -1, cut=[(attach, g) for g in gone]) != S:
                    continue
                f = self._loss(S, None, attach, radical, name, _lh, s, f'loss of {name} ({what})', None, 1,
                               mech=self.loss_mech(row, match, attach, None))
                self.add(f)
        for name, smarts, lost, radical, scores, what in FORMULA_LOSSES:
            s = scores.get(self.mode)
            patt = Chem.MolFromSmarts(smarts)
            if not s or patt is None:
                continue
            hits = m.mol.GetSubstructMatches(patt)
            if hits:
                hl = {m.owner[m.order[i]] for i in range(m.n)}
                ref = 'lipid' if 'phospho' in what else 'rda' if 'chromone' in what else 'glyco' if 'glycoside' in \
                    what else 'ei'
                self.add(self.formula_only(lost, radical, f'loss of {name} ({what})', s, hl,
                                           mech=(FORMULA_TEXT.get(what, ''), [], [], ref)))
        for counts, odd, (smarts, most), scores, what, text, ref in DIAGNOSTIC:
            s = scores.get(self.mode)
            patt = Chem.MolFromSmarts(smarts)
            if not s or patt is None or len(m.mol.GetSubstructMatches(patt)) < most:
                continue
            self.add(self.fixed_ion(counts, odd, what, s, {m.owner[m.order[i]] for i in range(m.n)},
                                    mech=(text, [], [], ref)))
        if self.z > 0 and not self.odd:
            self.immonium()

    def immonium(self):
        """Peptides / amino acids: the immonium ion H2N+=CH–R of each residue."""
        m = self.m
        aromatic = 0
        for ca in range(m.n):
            if not m.saturated(ca):
                continue
            n = next((j for j in m.heavy(ca) if m.el[j] == 'N' and not m.aro[j] and m.single(ca, j, True)), None)
            co = next((j for j in m.heavy(ca) if m.carbonyl_o(j) is not None and m.single(ca, j)), None)
            if n is None or co is None:
                continue
            prev = [j for j in m.heavy(n) if j != ca and m.carbonyl_o(j) is not None and m.single(n, j)]
            S = m.side(ca, co, cut=[(n, j) for j in prev])
            if S is None or any(j in S for j in prev):
                continue
            aromatic = any(m.aro[j] for j in S)
            self.add(self.make(S, {ca: (1, 0, -1)}, False, 'immonium ion of the residue',
                               0.7 if aromatic else 0.45, mech=(
                'An internal fragment (b → a ion of one residue): the amide bond before the residue and its Cα–C(=O) '
                'bond break, CO leaves, and the iminium ion H2N+=CH–R remains. It tells which amino acid is present '
                '(Phe 120.08, Tyr 136.08, Trp 159.09, His 110.07, Leu / Ile 86.10, Pro 70.07).',
                [('full', ('bond', (ca, co)), ('atom', co))], [(n, 1, 0, 1)], 'peptide')))

    def loss_mech(self, row, match, attach, parent):
        """(text, arrows, pre, ref) of losing a group: heterolysis toward the leaving group, or toward the ion for
        CO2, SO3, ketene, isobutene; homolysis for a radical."""
        m = self.m
        name, _s, rem, _att, _lh, radical, _sc, what = row
        g0 = next((match[k] for k in rem if attach in m.nb[match[k]]), match[rem[0]])
        texts = LOSS_TEXT.get(what, ('', None))
        text = texts[1] if self.mode == 'EI' and texts[1] else texts[0]
        if radical:
            arrows = [('half', ('bond', (attach, g0)), ('atom', g0)), ('half', ('bond', (attach, g0)), ('atom', attach))]
        elif name in ('CO2', 'SO3', 'HPO3', 'C2H2O', 'C4H8'):
            arrows = [('full', ('bond', (attach, g0)), ('atom', attach))]
        else:
            arrows = [('full', ('bond', (attach, g0)), ('atom', g0))]
        if parent is not None:
            pre = []
        elif self.mode == 'EI':
            pre = [(g0 if m.el[g0] in ('O', 'N', 'S', 'Cl', 'Br', 'I', 'F') else attach, 1, 1, 0)]
        elif self.z > 0:
            site = g0 if m.el[g0] in ('O', 'N', 'S') and name not in ('C2H2O', 'C4H8', 'CO2') else attach \
                if m.el[attach] in ('O', 'N') else self.base_site(frozenset(range(m.n)))
            pre = self.proton(site)
        else:
            pre = self.proton(None)
        ref = 'tms' if 'silyl' in what else 'lipid' if 'phosph' in what else 'esi' if self.mode != 'EI' else 'ei'
        return (text, arrows, pre, ref)

    def _loss(self, S, base, attach, radical, name, lh, s, how, parent, step, mech=None):
        """S after a group left from attach: the charge / radical and the H balance placed at attach."""
        m = self.m
        start = parent if parent is not None else None
        odd = (start.odd if start else self.odd) != radical
        lost = {}
        if start is None:
            gone = frozenset(range(m.n)) - S
        else:
            gone = start._idx - S
        for g in gone:
            if m.el[g] != 'H':
                lost[m.key[g]] = lost.get(m.key[g], 0) + 1
        lost['H'] = lh
        target = _sub(start.counts if start else self.pre, lost)
        hc = m.hcap(S)
        if hc is None or any(v < 0 for v in target.values()):
            return None
        have = hc + sum(v[2] for k, v in (base or {}).items() if k in S)
        dh = target.get('H', 0) - have
        if start is None:
            dq = self.z - sum(m.q[i] for i in S)
            rad = 1 if odd else 0
            if dq < 0 and (lh or radical):                  # [M−H]−: an acidic OH elsewhere carries the charge
                acid = self.acid_site(S, (attach,))
                if acid is not None:
                    base = {acid: (-1, 0, -1)}
                    dq = 0
                    dh += 1
            elif dq > 0 and self.mode != 'EI' and (radical or not lh):  # [M+H]+: the proton stays on a basic atom
                site = self.base_site(S, (attach,))
                if site is None:
                    return None
                base = {site: (1, 0, 1)}
                dq = 0
                dh -= 1
        else:                                               # the charge stays where it was
            if self.z > 0 and lh and not radical:           # ... unless the proton leaves with the molecule
                proton = next((k for k, v in base.items() if k in S and v[0] > 0 and v[2] > 0 and m.el[k] in 'NOS'),
                              None)
                if proton is not None:
                    v = base[proton]
                    base = dict(base)
                    base[proton] = (v[0] - 1, v[1], v[2] - 1)
                    dh += 1
            dq = self.z - sum(m.q[i] for i in S) - sum(v[0] for k, v in base.items() if k in S)
            rad = (1 if odd else 0) - sum(v[1] for k, v in base.items() if k in S)
            if rad < 0:                                     # a radical left: the unpaired electron goes with it
                base = {k: list(v) for k, v in base.items()}
                for k, v in base.items():
                    if v[1] > 0 and rad < 0:
                        v[1] -= 1; rad += 1
                base = {k: tuple(v) for k, v in base.items()}
        marks = self.place(S, base, attach, dq, rad, dh)
        if marks is None:
            return None
        return self.make(S, marks, odd, how, s, step, lost_name=name + ('' if not radical or name.endswith('•')
                                                                            else '•'), parent=parent, mech=mech)

    def next_losses(self, f, step):
        m = self.m
        out = []
        S0 = f._idx
        for row, matches in self._matches():
            name, _s, rem, att, _lh, radical, scores, what = row
            s = scores.get(self.mode)
            if not s or (radical and not f.odd):          # even-electron ions lose molecules, not radicals
                continue
            for match in matches:
                gone = {match[k] for k in rem}
                gone |= {j for g in gone for j in m.nb[g] if m.el[j] == 'H'}
                attach = match[att]
                if not gone <= S0 or attach not in S0 or any(g in f._imarks for g in gone):
                    continue
                if any(f._imarks.get(j, (0,))[0] for g in gone | {attach} for j in m.nb[g] + [g]):
                    continue                                # not the group that carries (or sits by) the charge
                S = S0 - gone
                if m.side(attach, -1, cut=[(attach, g) for g in gone] +
                          [(x, y) for x in S0 for y in m.nb[x] if y not in S0]) != S:
                    continue
                g = self._loss(S, f._imarks, attach, radical, name, _lh, f.score * s * 0.85,
                               f'{f.how}, then −{name}', f, step, mech=self.loss_mech(row, match, attach, f))
                if g is not None:
                    out.append(g)
        for c, v in f._imarks.items():                      # an acylium ion loses CO
            o = m.carbonyl_o(c)
            if v[0] <= 0 or o is None or o not in S0:
                continue
            rest = [j for j in m.heavy(c) if j in S0 and j != o]
            if len(rest) != 1 or m.el[rest[0]] != 'C':     # R–C≡O+ → R+ needs a carbon to take the charge
                continue
            S = S0 - {c, o}
            base = {k: w for k, w in f._imarks.items() if k != c}
            hc = m.hcap(S)
            if hc is None:
                continue
            target = _sub(f.counts, {'C': 1, 'O': 1})
            dh = target.get('H', 0) - hc - sum(w[2] for k, w in base.items() if k in S)
            marks = self.place(S, base, rest[0], v[0], 0, dh)
            g = self.make(S, marks, f.odd, f'{f.how}, then −CO', f.score * 0.75, step, lost_name='CO', parent=f,
                          mech=('The acylium ion R–C≡O+ expels carbon monoxide: the R–C bond pair goes with the CO '
                                'carbon and the carbocation R+ remains (benzoyl 105 → phenyl 77; peptide b → a ions).',
                                [('full', ('bond', (rest[0], c)), ('atom', c))], [], 'ei'))
            if g is not None:
                out.append(g)
        for x, v in f._imarks.items():                      # a carboxylate ion loses CO2
            c = next((j for j in m.heavy(x) if j in S0), None) if v[0] < 0 and m.el[x] == 'O' else None
            o = m.carbonyl_o(c) if c is not None else None
            if o is None or o not in S0:
                continue
            rest = [j for j in m.heavy(c) if j in S0 and j not in (o, x)]
            if len(rest) != 1:
                continue
            S = S0 - {c, o, x}
            base = {k: w for k, w in f._imarks.items() if k != x}
            hc = m.hcap(S)
            if hc is None:
                continue
            target = _sub(f.counts, {'C': 1, 'O': 2})
            dh = target.get('H', 0) - hc - sum(w[2] for k, w in base.items() if k in S)
            marks = self.place(S, base, rest[0], v[0], 0, dh)
            g = self.make(S, marks, f.odd, f'{f.how}, then −CO2', f.score * 0.8, step, lost_name='CO2', parent=f,
                          mech=('The carboxylate COO− expels CO2: the C–C bond pair stays on the carbon, which keeps '
                                'the negative charge (caffeate 179 → 135).',
                                [('full', ('bond', (rest[0], c)), ('atom', rest[0]))], [], 'chlorogenic'))
            if g is not None:
                out.append(g)
        if f.z > 0 and not f.odd:                           # an onium ion loses an alkene, its H moving to N / O
            for y in S0:
                if m.el[y] not in ('N', 'O') or m.aro[y]:
                    continue
                for c in m.heavy(y):
                    if c not in S0 or not m.saturated(c) or not m.single(c, y):
                        continue
                    R = m.side(c, y)
                    if R is None or not R <= S0 or any(k in f._imarks for k in R):
                        continue
                    if not any(m.h(j) for j in m.heavy(c) if j in R):
                        continue
                    S = S0 - R
                    hc = m.hcap(S)
                    if hc is None:
                        continue
                    kept = _clean(m.counts(S, hc + sum(w[2] for k, w in f._imarks.items() if k in S)))
                    lost = _clean(_sub(f.counts, kept))
                    bh = self.beta_h(c, R)
                    g = self.make(S, {k: w for k, w in f._imarks.items() if k in S}, f.odd,
                                  f'{f.how}, then −{chem.formula_text(lost)} (alkene; H to {m.el[y]})',
                                  f.score * 0.6, step, parent=f, mech=(
                        f'The onium ion loses an alkene: a hydrogen of the alkyl group moves to {m.el[y]} '
                        f'(four-membered transition state) as the {m.el[y]}–C bond breaks (triethylamine 86 → 58).',
                        [('full', ('bondH', bh), ('bond', (bh, c))), ('full', ('bond', (c, y)), ('atom', y))]
                        if bh is not None else [], [], 'ei'))
                    if g is not None:
                        out.append(g)
        if f.counts in ({'C': 7, 'H': 7}, {'C': 6, 'H': 5}) and f.z > 0 and not f.odd:
            g = self.formula_only({'C': 2, 'H': 2}, False, f'{f.how}, then −C2H2', f.score * 0.5, f.atoms, f, step,
                                  mech=('The tropylium (or phenyl) ion breaks its ring and loses acetylene: 91 → 65, '
                                        '77 → 51.', [], [], 'ei'))
            if g is not None:
                out.append(g)
        return out

    def merged(self):
        best = {}
        for f in self.out:
            key = (f.formula, f.z, f.odd)
            old = best.get(key)
            if old is None:
                best[key] = f
                continue
            if f.score > old.score:
                f.also = [old.how] + [h for h in old.also if h != f.how]
                best[key] = f
            elif f.how != old.how and f.how not in old.also:
                old.also.append(f.how)
        out = list(best.values())
        for f in out:                                       # agreeing rules make an ion more likely
            f.score = min(1.0, f.score + 0.05 * min(len(f.also), 3))
        out.sort(key=lambda f: -f.mz)
        return out[:500]


def predict(doc, atom_ids=None, mode='ESI+', steps=2):
    """(precursor Fragment, [Fragment] highest m/z first) of the structure(s) with these atoms; (None, []) when
    the structure cannot be read. A structure drawn with a charge is used as it is (mode is ignored)."""
    m = _Mol(doc, atom_ids)
    if m.mol is None or not m.n:
        return None, []
    e = _Engine(m, mode)
    if e.pre.get('H', 0) < 0:
        return None, []
    pre = Fragment(e.mz(e.pre, e.z), e.z, e.odd, e.pre, ion_text(e.pre, e.z, e.odd), '', f'precursor ion {e.label}',
                   1.0, atoms=frozenset(m.owner[m.order[i]] for i in range(m.n)), step=0)
    pre.label = e.label
    pre.drawn_ion = e.charged
    frags = e.run(steps)
    everything = pre.atoms
    for f in frags:
        f.lost_atoms = frozenset(everything - f.atoms) if f.S is not None else frozenset()
    return pre, frags


# ------------------------------------------------------------------------------------------- drawing a fragment
def draw_doc(doc, frag):
    """(Doc holding the fragment ion, '') with the drawing's atom ids and places, its charge, unpaired electron and
    H as worked out; (None, why) when it cannot be drawn."""
    from model import Doc
    if frag.S is None:
        return None, 'This fragment is worked out from the formula only, so it has no structure to draw.'
    full, grown = doc.expanded()
    owners = {}
    for a, more in grown.items():
        for g in more:
            owners[g] = a
    piece = set()
    for fid in frag.S:
        o = owners.get(fid, fid)
        if o not in doc.atoms:
            return None, 'The fragment could not be matched to the drawing.'
        piece.add(o)
    for o in piece:
        group = {o} | grown.get(o, set())
        if not group <= frag.S:
            return None, f'This fragment splits the {doc.atoms[o].label} label: expand it first, then draw it.'
    for fid in frag.marks:
        o = owners.get(fid, fid)
        if fid != o or doc.atoms[o].label:
            return None, f'The charge of this fragment sits inside the {doc.atoms[o].label} label: expand it first.'
    d = Doc()
    d.style = dict(doc.style)
    d.insert_dict(doc.to_dict([('atom', a) for a in piece]), keep_ids=True)
    target = {}
    for a in d.atoms:
        if d.atoms[a].label:
            continue
        cut = sum(doc.bonds[b].order for _n, b in doc.neighbors(a) if _n not in piece)
        target[a] = doc.implicit_h(a) + cut + frag.marks.get(a, (0, 0, 0))[2]
    for a, (dq, rad, _dh) in frag.marks.items():
        d.atoms[a].charge += dq
        d.atoms[a].radical = min(2, d.atoms[a].radical + rad)
    for a in target:
        d.atoms[a].hcount = None

    def natural(a):
        return d.implicit_h(a)
    deficit = {a: natural(a) - t for a, t in target.items() if t < natural(a)}
    for a in sorted(deficit, key=lambda k: -deficit[k]):        # too few H for the bonds: make a double bond
        guard = 0
        while deficit.get(a, 0) > 0 and guard < 6:
            guard += 1
            nbrs = [(n, b) for n, b in d.neighbors(a) if not d.atoms[n].label and d.bonds[b].order < 3]
            pair = next(((n, b) for n, b in nbrs if deficit.get(n, 0) > 0), None)
            if pair is None and deficit[a] >= 2:
                spare = next(((n, b) for n, b in nbrs if target.get(n, 0) >= 1 and n not in deficit), None)
                if spare is not None:                       # one of a's missing H moves to the neighbour
                    target[spare[0]] -= 1
                    target[a] += 1
                    deficit[spare[0]] = 1
                    deficit[a] -= 1
                    pair = spare
            if pair is None:
                break
            n, b = pair
            d.bonds[b].order += 1
            deficit[a] -= 1
            deficit[n] -= 1
    for a in [a for a in deficit if deficit[a] > 0]:      # the two ends of a new diene: a–b=c–e becomes a=b–c=e
        if deficit.get(a, 0) <= 0:
            continue
        for b, ab in d.neighbors(a):
            if deficit.get(b, 0) > 0 or d.bonds[ab].order >= 3:
                continue
            for c, bc in d.neighbors(b):
                if c == a:
                    continue
                if d.bonds[bc].order != 2:
                    _flip_ring(d, bc)
                if d.bonds[bc].order != 2:
                    continue
                e = next(((e, ce) for e, ce in d.neighbors(c) if e not in (a, b) and deficit.get(e, 0) > 0
                          and d.bonds[ce].order < 3), None)
                if e is None:
                    continue
                d.bonds[bc].order = 1
                d.bonds[ab].order += 1
                d.bonds[e[1]].order += 1
                deficit[a] -= 1
                deficit[e[0]] -= 1
                break
            if deficit.get(a, 0) <= 0:
                break
    for a, t in target.items():
        if natural(a) != t:
            d.atoms[a].hcount = max(0, t)
    if not lewis(d, frag.z, frag.odd):
        return None, ('Drawing this fragment needs bonds moved round (a quinoid or rearranged structure); it is '
                      'listed, but not drawn automatically.')
    return d, ''


# ------------------------------------------------------------------------------------------- Lewis structures
def _h(d, a):
    return d.implicit_h(a)


def _set_h(d, a, h):
    """Atom a keeps h hydrogens, written out only when its bonds and charge would give another number."""
    at = d.atoms[a]
    at.hcount = None
    if d.implicit_h(a) != h:
        at.hcount = max(0, h)


def _valid(d):
    """Every atom within its valence, and RDKit reads the structure."""
    import chem
    if any(d.valence_problem(a) for a in d.atoms):
        return False
    mol, _order, _err = chem.to_mol(d)
    return mol is not None


def _lone_pair(d, x):
    at = d.atoms[x]
    if at.label or at.charge or at.radical or at.el not in ('N', 'O', 'S'):
        return False
    return d.bond_sum(x) + d.implicit_h(x) <= {'N': 3, 'O': 2, 'S': 2}[at.el]


def lewis(d, z, odd):
    """Redraws a fragment ion as chemists write it, keeping each atom's H: a structure that is not valid is solved
    again by RDKit (bond orders and charges, the unpaired electron put back afterwards); then a carbocation next to
    N, O or S becomes the onium C=X+ (an acylium C≡O+), a protonated OH of an acid or ester moves to its C=O
    oxygen, a carbanion gives its charge to an acidic OH (carboxylate, phenolate), and a radical next to a charge
    becomes a C=C with •+ (•−) written beside it. False when no valid structure is found."""
    if not d.atoms:
        return False
    if not _valid(d) and not _solve(d, z, odd):
        return False
    for _round in range(len(d.atoms)):
        if not (_onium(d) or _proton_to_carbonyl(d) or _carbanion_to_acid(d) or _radical_pairs(d)):
            break
    return _valid(d)


def _solve(d, z, odd):
    """Bond orders and formal charges from the atoms, their bonds and H (RDKit DetermineBondOrders); an odd-electron
    ion gets a temporary H on its radical atom first."""
    from rdkit import Chem
    try:
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        return False
    ids = [a for a in d.atoms if not d.atoms[a].label]
    if len(ids) != len(d.atoms):
        return False
    hs = {a: _h(d, a) for a in ids}
    rad = None
    if odd:
        rad = next((a for a in ids if d.atoms[a].radical), None) or next(
            (a for a in ids if d.atoms[a].el == 'C' and d.atoms[a].charge), None) or next(
            (a for a in ids if d.atoms[a].el == 'C'), None)
        if rad is None:
            return False
    rw = Chem.RWMol()
    index = {}
    for a in ids:
        atom = Chem.Atom(d.atoms[a].el)
        atom.SetNoImplicit(True)
        index[a] = rw.AddAtom(atom)
    for a in ids:
        for _k in range(hs[a] + (1 if a == rad else 0)):
            h = rw.AddAtom(Chem.Atom(1))
            rw.AddBond(index[a], h, Chem.BondType.SINGLE)
    for i, b in d.bonds.items():
        rw.AddBond(index[b.a], index[b.b], Chem.BondType.SINGLE)
    try:
        rdDetermineBonds.DetermineBondOrders(rw, charge=z, allowChargedFragments=True, embedChiral=False)
    except Exception:
        return False
    order = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
    for i, b in d.bonds.items():
        bd = rw.GetBondBetweenAtoms(index[b.a], index[b.b])
        if bd.GetBondType() == Chem.BondType.AROMATIC:
            return False
        b.order = order.get(bd.GetBondType(), 1)
    for a in ids:
        at = d.atoms[a]
        at.charge = rw.GetAtomWithIdx(index[a]).GetFormalCharge()
        at.radical = 1 if a == rad else 0
    for a in ids:
        _set_h(d, a, hs[a])
    return _valid(d)


def _onium(d):
    """C+ next to a lone pair: C=N+, C=O+, C=S+; an acylium C≡O+."""
    for a, at in d.atoms.items():
        if at.label or at.el != 'C' or at.charge != 1 or at.radical:
            continue
        o = next(((n, b) for n, b in d.neighbors(a) if d.atoms[n].el == 'O' and d.bonds[b].order == 2
                  and not d.atoms[n].charge), None)
        if o is None:
            ranked = sorted(((n, b) for n, b in d.neighbors(a) if d.bonds[b].order == 1 and _lone_pair(d, n)),
                            key=lambda nb: 'NSO'.index(d.atoms[nb[0]].el))
            o = ranked[0] if ranked else None
        if o is None:
            continue
        n, b = o
        ha, hn = _h(d, a), _h(d, n)
        d.bonds[b].order += 1
        at.charge, d.atoms[n].charge = 0, 1
        _set_h(d, a, ha); _set_h(d, n, hn)
        if _valid(d):
            return True
        d.bonds[b].order -= 1                               # not possible here: as it was
        at.charge, d.atoms[n].charge = 1, 0
        _set_h(d, a, ha); _set_h(d, n, hn)
    return False


def _proton_to_carbonyl(d):
    """R–C(=O)–OH2+ is written R–C(OH)=OH+ (the proton on the carbonyl oxygen)."""
    for x, at in d.atoms.items():
        if at.label or at.el != 'O' or at.charge != 1 or _h(d, x) < 1:
            continue
        for c, b in d.neighbors(x):
            if d.bonds[b].order != 1 or d.atoms[c].el != 'C':
                continue
            o = next((n for n, bb in d.neighbors(c) if n != x and d.atoms[n].el == 'O' and d.bonds[bb].order == 2
                      and not d.atoms[n].charge), None)
            if o is None:
                continue
            hx, ho = _h(d, x), _h(d, o)
            at.charge, d.atoms[o].charge = 0, 1
            _set_h(d, x, hx - 1); _set_h(d, o, ho + 1)
            if _valid(d):
                return True
            at.charge, d.atoms[o].charge = 1, 0
            _set_h(d, x, hx); _set_h(d, o, ho)
    return False


def _carbanion_to_acid(d):
    """A carbanion takes the proton of an acidic OH in the same ion: the charge sits on the carboxylate / phenolate."""
    for a, at in d.atoms.items():
        if at.label or at.el != 'C' or at.charge != -1 or at.radical:
            continue
        best = None
        for o, oa in d.atoms.items():
            if oa.label or oa.el != 'O' or oa.charge or _h(d, o) < 1 or len(d.neighbors(o)) != 1:
                continue
            c = d.neighbors(o)[0][0]
            carboxyl = any(d.atoms[n].el == 'O' and d.bonds[bb].order == 2 for n, bb in d.neighbors(c))
            aromatic = any(d.bonds[bb].order == 2 for _n, bb in d.neighbors(c))
            rank = 0 if carboxyl else 1 if aromatic else None
            if rank is not None and (best is None or rank < best[0]):
                best = (rank, o)
        if best is None:
            continue
        o = best[1]
        ha, ho = _h(d, a), _h(d, o)
        at.charge, d.atoms[o].charge = 0, -1
        _set_h(d, a, ha + 1); _set_h(d, o, ho - 1)
        if _valid(d):
            return True
        at.charge, d.atoms[o].charge = -1, 0
        _set_h(d, a, ha); _set_h(d, o, ho)
    return False


def _radical_pairs(d):
    """An unpaired electron next to a charged carbon: one more bond between them, •+ (or •−) written beside it."""
    for a, at in d.atoms.items():
        if at.label or not at.radical:
            continue
        for n, b in d.neighbors(a):
            na = d.atoms[n]
            if na.label or not na.charge or na.radical or d.bonds[b].order >= 3 or na.el not in ('C', 'N', 'O'):
                continue
            ha, hn, q = _h(d, a), _h(d, n), na.charge
            d.bonds[b].order += 1
            at.radical, na.charge = 0, 0
            _set_h(d, a, ha); _set_h(d, n, hn)
            if _valid(d):
                mx, my = (at.x + na.x) / 2, (at.y + na.y) / 2
                d.add_text(mx + 0.3, my - 0.9, '•' + ('+' if q > 0 else '−'), size=0)
                return True
            d.bonds[b].order -= 1
            at.radical, na.charge = 1, q
            _set_h(d, a, ha); _set_h(d, n, hn)
        for b1, bb1 in d.neighbors(a):                      # 1,4: a–b=c–n becomes a=b–c=n (a diene radical ion)
            if d.bonds[bb1].order != 1:
                continue
            for c1, bb2 in d.neighbors(b1):
                if c1 == a or d.bonds[bb2].order != 2:
                    continue
                for n, bb3 in d.neighbors(c1):
                    na = d.atoms[n]
                    if n in (a, b1) or d.bonds[bb3].order != 1 or na.label or not na.charge or na.radical:
                        continue
                    hs = {k: _h(d, k) for k in (a, b1, c1, n)}
                    q = na.charge
                    d.bonds[bb1].order, d.bonds[bb2].order, d.bonds[bb3].order = 2, 1, 2
                    at.radical, na.charge = 0, 0
                    for k, h in hs.items():
                        _set_h(d, k, h)
                    if _valid(d):
                        xs = [d.atoms[k].x for k in (a, b1, c1, n)]
                        ys = [d.atoms[k].y for k in (a, b1, c1, n)]
                        d.add_text(max(xs) + 0.3, min(ys) - 0.9, '•' + ('+' if q > 0 else '−'), size=0)
                        return True
                    d.bonds[bb1].order, d.bonds[bb2].order, d.bonds[bb3].order = 1, 2, 1
                    at.radical, na.charge = 1, q
                    for k, h in hs.items():
                        _set_h(d, k, h)
    return False


def _flip_ring(d, bond):
    """Swaps single and double bonds round a six-membered ring of alternating bonds holding bond (the other Kekulé
    structure of a benzene ring), so that bond becomes double."""
    a0, b0 = d.bonds[bond].a, d.bonds[bond].b
    for path in _rings6(d, a0, b0):
        ids = [d.bond_between(path[k], path[(k + 1) % 6]) for k in range(6)]
        orders = [d.bonds[i].order for i in ids]
        if sorted(orders) == [1, 1, 1, 2, 2, 2] and all(orders[k] != orders[(k + 1) % 6] for k in range(6)):
            for i in ids:
                d.bonds[i].order = 3 - d.bonds[i].order
            return


def _rings6(d, a, b):
    """Six-membered rings through bond a–b, as atom lists starting a, b."""
    out = []

    def walk(path):
        if len(path) == 6:
            if d.bond_between(path[-1], a) is not None:
                out.append(list(path))
            return
        for n, _bd in d.neighbors(path[-1]):
            if n not in path:
                walk(path + [n])
    walk([a, b])
    return out


# ------------------------------------------------------------------------------------------- drawing a mechanism
def chain_of(frag):
    """The fragments from the precursor's first product to frag, in order."""
    out = []
    while frag is not None:
        out.insert(0, frag)
        frag = frag.parent
    return out


def step_loss(f):
    """What one step loses ('CH3•', 'H2O'), from its fragment's running list."""
    return f.lost.split(' + ')[-1] if f.lost else ''


def mechanism(doc, frag, atoms):
    """The fragmentation drawn step by step: ([(Doc, loss after it)], why the chain stops early or ''). The first
    Doc is the precursor (the structure with atom ids atoms, its charge where the first step starts), the others the
    ions along the way; each holds, as curved arrows tied to its atoms and bonds (fish-hooks for one electron), the
    step that follows it - with explicit H atoms where a hydrogen moves."""
    from model import Doc
    chain = chain_of(frag)
    first = Doc()
    first.style = dict(doc.style)
    first.insert_dict(doc.to_dict([('atom', a) for a in atoms if a in doc.atoms]), keep_ids=True)
    _mark(first, chain[0].pre)
    _arrows(first, chain[0].arrows)
    states = [(first, step_loss(chain[0]))]
    for k, f in enumerate(chain):
        piece, why = draw_doc(doc, f)
        if piece is None:
            return states, why
        if k + 1 < len(chain):
            _mark(piece, chain[k + 1].pre)
            _arrows(piece, chain[k + 1].arrows)
        states.append((piece, step_loss(chain[k + 1]) if k + 1 < len(chain) else ''))
    return states, ''


def _mark(d, pre):
    """The charge (and unpaired electron, and the extra or missing H) where a step starts, on its atom - or, when
    the atom cannot carry it as drawn (an aromatic carbon with •+), written beside it."""
    for a, dq, dr, dh in pre or []:
        at = d.atoms.get(a)
        if at is None or at.label:
            continue
        keep = (at.charge, at.radical, at.hcount)
        h = d.implicit_h(a) + dh
        at.charge += dq
        at.radical = min(2, at.radical + dr)
        at.hcount = None
        if d.implicit_h(a) != h:
            at.hcount = max(0, h)
        if not _valid(d):
            at.charge, at.radical, at.hcount = keep
            sign = ('•' if dr else '') + ('+' if dq > 0 else '−' if dq < 0 else '')
            if sign:
                d.add_text(at.x + 0.35, at.y - 0.95, sign, size=0)


def _h_atom(d, a, made):
    """An explicit H on atom a (made once), pointing away from its neighbours."""
    import math
    if a in made:
        return made[a]
    at = d.atoms[a]
    vx = vy = 0.0
    for n, _b in d.neighbors(a):
        vx += d.atoms[n].x - at.x
        vy += d.atoms[n].y - at.y
    L = math.hypot(vx, vy)
    ux, uy = (-vx / L, -vy / L) if L > 1e-6 else (0.0, -1.0)
    if L <= 1e-6 and d.neighbors(a):
        n = d.atoms[d.neighbors(a)[0][0]]
        ux, uy = -(n.y - at.y), n.x - at.x
    h = d.add_atom(at.x + 0.8 * ux, at.y + 0.8 * uy, 'H')
    if at.hcount:
        at.hcount -= 1
    d.add_bond(a, h)
    made[a] = h
    return h


def _arrows(d, arrows):
    """Curved arrows (full or fish-hook) tied to the atoms and bonds of d; ends that are not in d are skipped."""
    import model
    made = {}

    def item(e):
        kind, v = e
        if kind == 'atom':
            return ('atom', v) if v in d.atoms and not d.atoms[v].label else None
        if kind == 'bond':
            b = d.bond_between(*v) if all(x in d.atoms for x in v) else None
            return ('bond', b) if b is not None else None
        if v not in d.atoms or d.atoms[v].label:
            return None
        h = _h_atom(d, v, made)
        return ('atom', h) if kind == 'H' else ('bond', d.bond_between(v, h))
    ends = [(k, item(a), item(b)) for k, a, b in arrows]
    if not d.atoms:
        return
    cx = sum(at.x for at in d.atoms.values()) / len(d.atoms)
    cy = sum(at.y for at in d.atoms.values()) / len(d.atoms)
    for kind, a, b in ends:
        if a is None or b is None or a == b:
            continue
        (x1, y1), (x2, y2) = d.anchor(*a), d.anchor(*b)
        if a[0] == 'atom':                                  # a lone pair: start just off the atom, away from the head
            dx, dy = x1 - x2, y1 - y2
            L = max((dx * dx + dy * dy) ** 0.5, 1e-6)
            x1, y1 = x1 + 0.22 * dx / L, y1 + 0.22 * dy / L
        if b[0] == 'atom' and d.atoms[b[1]].el != 'H':      # to an atom: stop just short of its label
            dx, dy = x1 - x2, y1 - y2
            L = max((dx * dx + dy * dy) ** 0.5, 1e-6)
            x2, y2 = x2 + 0.2 * dx / L, y2 + 0.2 * dy / L
        if ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5 < 0.15:
            continue
        span = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5
        ar = model.Arrow(x1, y1, x2, y2, 'curved' if kind == 'full' else 'curved_half', 0.3 if span < 1.2 else 0.16)
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        px, py = model.curve_control(ar)
        if (px - cx) ** 2 + (py - cy) ** 2 < (mx - cx) ** 2 + (my - cy) ** 2:
            ar.bend = -ar.bend                              # bulge away from the middle of the structure
        i = d._id()
        d.arrows[i] = ar
        d.attach_arrow_end(ar, 'tail', a)
        d.attach_arrow_end(ar, 'head', b)
