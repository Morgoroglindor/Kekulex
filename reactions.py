"""Chemistry on the drawing beyond one structure: reaction atom maps and balance, tautomers and resonance forms,
substructure search and stereo checks (RDKit)."""
from rdkit import Chem
from rdkit.Chem import rdFMCS

import chem
import model


# ------------------------------------------------------------------------------------------- molecules with their ids
def mol_of(doc, atom_ids, expand=True):
    """(RDKit molecule, [drawing atom id for each molecule atom]) of these atoms. With expand, nicknames are written
    out and their atoms stand for the nickname's own atom id."""
    if expand and any(doc.atoms[a].label for a in atom_ids if a in doc.atoms):
        full, grown = doc.expanded()
        back = {}
        for a in atom_ids:
            back[a] = a
            for g in grown.get(a, ()):
                back[g] = a
        ids = [a for a in full.atoms if a in back]
        mol, order, _err = chem.to_mol(full, set(ids))
        return mol, [back.get(a, a) for a in order]
    mol, order, _err = chem.to_mol(doc, set(atom_ids))
    return mol, order


def _with_coords(doc, mol, order, src=None):
    """mol with a 2D conformer from the drawing's positions (y up, 1.5 A bonds) so it keeps its look when inserted."""
    src = src or doc
    conf = Chem.Conformer(mol.GetNumAtoms())
    for k, a in enumerate(order):
        at = src.atoms[a]
        conf.SetAtomPosition(k, (at.x * model.MOL_BOND, -at.y * model.MOL_BOND, 0.0))
    conf.Set3D(False)
    mol = Chem.Mol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


# ------------------------------------------------------------------------------------------- tautomers / resonance
def _laid_out(doc, atom_ids):
    """(molecule, order, the doc whose positions it has) - nicknames written out, positions as drawn."""
    if any(doc.atoms[a].label for a in atom_ids):
        full, grown = doc.expanded()
        ids = set(chem.expand_ids(grown, atom_ids))
        mol, order, _err = chem.to_mol(full, ids)
        return mol, order, full
    mol, order, _err = chem.to_mol(doc, set(atom_ids))
    return mol, order, doc


def tautomers(doc, atom_ids, most=8):
    """Other tautomers of the structure (RDKit's enumerator), as molecules laid out like the drawing; the drawn form
    is left out."""
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol, order, src = _laid_out(doc, atom_ids)
    if mol is None:
        return []
    en = rdMolStandardize.TautomerEnumerator()
    en.SetMaxTautomers(64)
    drawn = Chem.MolToSmiles(mol)
    out, seen = [], {drawn}
    for t in en.Enumerate(mol):
        smi = Chem.MolToSmiles(t)
        if smi in seen:
            continue
        seen.add(smi)
        out.append(_with_coords(doc, t, order, src))
        if len(out) >= most:
            break
    return out


def resonance_forms(doc, atom_ids, most=8):
    """Other resonance structures (RDKit's ResonanceMolSupplier: charges and lone pairs moved along conjugation),
    laid out like the drawing; the drawn form is left out."""
    mol, order, src = _laid_out(doc, atom_ids)
    if mol is None:
        return []
    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception:
        pass
    flags = (Chem.ALLOW_CHARGE_SEPARATION | Chem.UNCONSTRAINED_CATIONS | Chem.UNCONSTRAINED_ANIONS |
             Chem.KEKULE_ALL)
    sup = Chem.ResonanceMolSupplier(mol, flags)
    seen = {_resonance_key(mol)}
    hs = [a.GetTotalNumHs() for a in mol.GetAtoms()]
    charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
    out, kekule = [], []
    for k in range(min(len(sup), 96)):
        r = sup[k]
        if r is None or not _octets_ok(r):
            continue
        try:
            r = Chem.Mol(r)
            Chem.SanitizeMol(r, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception:
            continue
        if any(r.GetAtomWithIdx(i).GetTotalNumHs() != hs[i] for i in range(len(hs))):
            continue                                        # resonance moves electrons, never hydrogens
        new = [(a.GetAtomicNum(), a.GetFormalCharge()) for a in r.GetAtoms()
               if a.GetFormalCharge() != charges[a.GetIdx()]]
        if ((6, 1) in new and (6, -1) in new) or new.count((6, 1)) > 1 or new.count((6, -1)) > 1:
            continue                                        # C+ and C- made on carbon alone (benzene): not drawn
        key = _resonance_key(r)
        if key in seen:                                     # the same charges and the same non-ring bonds: at most
            if _bond_key(r) != _bond_key(mol) and not kekule:  # the other Kekulé form of a ring (benzene)
                kekule.append(r)
            continue
        seen.add(key)
        out.append(_with_coords(doc, r, order, src))
        if len(out) >= most:
            break
    if not out and kekule:
        out.append(_with_coords(doc, kekule[0], order, src))
    return out


def _octets_ok(mol):
    """No nitrogen, oxygen, sulfur or halogen short of its octet (an O+ with one bond and no H is not a form a
    chemist draws); carbocations are fine."""
    for a in mol.GetAtoms():
        if a.GetAtomicNum() in (1, 6) or a.GetAtomicNum() > 18 and a.GetAtomicNum() not in (35, 53):
            continue
        ve = Chem.GetPeriodicTable().GetNOuterElecs(a.GetAtomicNum())
        shared = sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) + a.GetTotalNumHs()
        if ve - a.GetFormalCharge() + shared - a.GetNumRadicalElectrons() < 8 and a.GetAtomicNum() not in (5, 13):
            return False
    return True


def _bond_key(mol):
    m = Chem.Mol(mol)
    try:
        Chem.Kekulize(m, clearAromaticFlags=True)
    except Exception:
        pass
    return tuple(sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                         b.GetBondTypeAsDouble()) for b in m.GetBonds()))


def _resonance_key(mol):
    """Charges atom by atom and the orders of the bonds outside aromatic rings: Kekulé forms of one ring count as one."""
    m = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    ring = {b.GetIdx() for b in m.GetBonds() if b.GetIsAromatic()}
    m2 = Chem.Mol(m)
    try:
        Chem.Kekulize(m2, clearAromaticFlags=True)
    except Exception:
        pass
    return (tuple(a.GetFormalCharge() for a in m2.GetAtoms()),
            tuple(sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                          b.GetBondTypeAsDouble()) for b in m2.GetBonds() if b.GetIdx() not in ring)))


# ------------------------------------------------------------------------------------------- reactions
def sides(doc, arrow=None):
    """([reactant atom sets], [reagent atom sets above / below the arrow], [product atom sets]) or None."""
    arrows = [i for i, a in doc.arrows.items()
              if not a.kind.startswith('curved') and abs(a.x2 - a.x1) > abs(a.y2 - a.y1)]   # across, not down
    if arrow is None:
        if not arrows:
            return None
        arrow = arrows[0]
    ar = doc.arrows[arrow]
    lo, hi = sorted((ar.x1, ar.x2))
    ay = (ar.y1 + ar.y2) / 2
    forward = ar.x2 >= ar.x1
    # this step only: its row (structures reaching within 4 bond lengths of the arrow's line), between the arrows
    # before and after it on that row
    row = [o for i, o in doc.arrows.items() if i != arrow and i in arrows and abs((o.y1 + o.y2) / 2 - ay) < 2.5]
    before = max((max(o.x1, o.x2) for o in row if max(o.x1, o.x2) <= lo), default=-1e9)
    after = min((min(o.x1, o.x2) for o in row if min(o.x1, o.x2) >= hi), default=1e9)
    left, over, prods = [], [], []
    for f in doc.fragments():
        ys = [doc.atoms[a].y for a in f]
        if max(ys) < ay - 4 or min(ys) > ay + 4:
            continue
        cx = sum(doc.atoms[a].x for a in f) / len(f)
        if lo <= cx <= hi:
            over.append(set(f))
        elif (cx < lo) == forward:
            if (cx > before) if forward else (cx < after):
                left.append(set(f))
        elif (cx < after) if forward else (cx > before):
            prods.append(set(f))
    if not prods and not left:
        return None
    return left, over, prods


def map_atoms(doc, arrow=None, timeout=2):
    """Numbers matching atoms on both sides of the reaction (Atom.map), from the largest common substructure of each
    product with each reactant (elements must match, ring atoms only to ring atoms). Returns how many were mapped."""
    found = sides(doc, arrow)
    if not found:
        return 0
    left, over, prods = found
    for a in doc.atoms.values():
        a.map = 0
    reactants = []
    for f in left + over:
        mol, order = mol_of(doc, f, expand=False)
        if mol is not None:
            reactants.append((mol, order, set()))
    nxt = 1
    for f in sorted(prods, key=len, reverse=True):
        pmol, porder = mol_of(doc, f, expand=False)
        if pmol is None:
            continue
        taken = set()
        while len(taken) < pmol.GetNumAtoms():
            psub, pkeep = _submol(pmol, taken)              # what is still unmapped in the product
            best = None
            for rmol, rorder, used in reactants:
                if len(used) >= rmol.GetNumAtoms():
                    continue
                rsub, rkeep = _submol(rmol, used)
                res = rdFMCS.FindMCS([rsub, psub], atomCompare=rdFMCS.AtomCompare.CompareElements,
                                     bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=True,
                                     timeout=timeout)
                if res.numAtoms < 2 or (best is not None and res.numAtoms <= best[0]):
                    continue
                q = Chem.MolFromSmarts(res.smartsString)
                if q is None:
                    continue
                rms = rsub.GetSubstructMatches(q, uniquify=False, maxMatches=24)
                pms = psub.GetSubstructMatches(q, uniquify=False, maxMatches=24)
                pair = None
                for rm in rms:                              # the pairing whose atoms keep their number of
                    for pm in pms:                          # neighbours best (CH3 to CH3, CH2-O to CH2-O)
                        cost = sum(abs(rmol.GetAtomWithIdx(rkeep[i]).GetDegree() -
                                       pmol.GetAtomWithIdx(pkeep[j]).GetDegree()) for i, j in zip(rm, pm))
                        if pair is None or cost < pair[0]:
                            pair = (cost, rm, pm)
                if pair:
                    _c, rm, pm = pair
                    best = (res.numAtoms, [rkeep[k] for k in rm], [pkeep[k] for k in pm], rorder, used)
            if best is None:
                break
            _n, rmatch, pmatch, rorder, used = best
            for ri, pi in zip(rmatch, pmatch):
                doc.atoms[rorder[ri]].map = nxt
                doc.atoms[porder[pi]].map = nxt
                nxt += 1
            used |= set(rmatch)
            taken |= set(pmatch)
    return nxt - 1


def _submol(mol, leave_out):
    """(mol without the atoms in leave_out, [original index of each atom kept])."""
    rw = Chem.RWMol(mol)
    for k in sorted(leave_out, reverse=True):
        rw.RemoveAtom(k)
    sub = rw.GetMol()
    sub.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(sub)
    keep = [k for k in range(mol.GetNumAtoms()) if k not in leave_out]
    return sub, keep


SMALL = {'H2O': {'H': 2, 'O': 1}, 'H2': {'H': 2}, 'HCl': {'H': 1, 'Cl': 1}, 'HBr': {'H': 1, 'Br': 1},
         'HI': {'H': 1, 'I': 1}, 'HF': {'H': 1, 'F': 1}, 'CO2': {'C': 1, 'O': 2}, 'CO': {'C': 1, 'O': 1},
         'NH3': {'N': 1, 'H': 3}, 'N2': {'N': 2}, 'O2': {'O': 2}, 'CH3OH': {'C': 1, 'H': 4, 'O': 1},
         'C2H5OH': {'C': 2, 'H': 6, 'O': 1}, 'CH3COOH': {'C': 2, 'H': 4, 'O': 2}, 'SO2': {'S': 1, 'O': 2},
         'NaCl': {'Na': 1, 'Cl': 1}, 'NaBr': {'Na': 1, 'Br': 1}, 'KBr': {'K': 1, 'Br': 1}, 'HOAc': {}}


def balance(doc, arrow=None):
    """(text, {element: products - reactants}) for the reaction: the reactants left of the arrow against the
    products (reagents over / under the arrow are left out). Names a by-product that would balance it."""
    found = sides(doc, arrow)
    if not found:
        return 'Draw a reaction first: structures left and right of a reaction arrow.', {}
    left, over, prods = found

    def total(frags):
        out = {}
        for f in frags:
            for e, n in chem.analyse(doc, f).counts.items():
                out[e] = out.get(e, 0) + n
        return out
    r, p = total(left), total(prods)
    diff = {e: p.get(e, 0) - r.get(e, 0) for e in set(r) | set(p)}
    diff = {e: n for e, n in diff.items() if n}
    if not diff:
        return 'Balanced: the products hold the same atoms as the reactants.', {}
    lost = {e: -n for e, n in diff.items() if n < 0}
    gained = {e: n for e, n in diff.items() if n > 0}
    parts = []
    if lost:
        parts.append('the products have ' + chem.formula_text(lost) + ' fewer')
    if gained:
        parts.append('the products have ' + chem.formula_text(gained) + ' more')
    text = 'Not balanced: ' + ' and '.join(parts) + '.'
    for name, f in SMALL.items():
        if not f:
            continue
        for k in (1, 2, 3):
            if not gained and lost == {e: n * k for e, n in f.items()}:
                text += f' A by-product of {k if k > 1 else ""}{name} would balance it.'
                break
            if not lost and gained == {e: n * k for e, n in f.items()}:
                text += f' {k if k > 1 else ""}{name} more on the left would balance it.'
                break
    if over:
        text += ' (Reagents over / under the arrow are not counted.)'
    return text, diff


# ------------------------------------------------------------------------------------------- find
def find(doc, text):
    """Items matching text: captions holding it, atoms with that label or element, else a SMILES / SMARTS
    substructure in the structures (with the bonds between the atoms found). [(kind, id)]."""
    t = text.strip()
    if not t:
        return []
    low = t.lower()
    hits = [('text', i) for i, x in doc.texts.items() if low in x.text.lower()]
    for a, at in doc.atoms.items():
        if (at.label and at.label.lower() == low) or (not at.label and at.el.lower() == low and len(t) <= 2):
            hits.append(('atom', a))
    if hits:
        return hits
    q = Chem.MolFromSmarts(t)
    if q is None or q.GetNumAtoms() == 0:
        m = Chem.MolFromSmiles(t)
        q = m if m is not None else None
    if q is None:
        return []
    found = set()
    for frag in doc.fragments():
        mol, order = mol_of(doc, frag, expand=True)
        if mol is None:
            continue
        for match in mol.GetSubstructMatches(q, maxMatches=200):
            found |= {order[k] for k in match}
    out = [('atom', a) for a in found if a in doc.atoms]
    out += [('bond', b) for b, bd in doc.bonds.items() if bd.a in found and bd.b in found]
    return out


# ------------------------------------------------------------------------------------------- stereo checks
def stereo_problems(doc, atom_ids):
    """Readable stereo problems: stereocentres with no wedge, double bonds that could be E or Z drawn straight
    through (both ends with two different groups) are fine; wedges on atoms that are not stereocentres, and more
    than one wedge from an atom pointing the same way are flagged."""
    out = []
    mol, order = mol_of(doc, atom_ids, expand=False)
    if mol is None:
        return out
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    undefined = [k for k, code in centres if code == '?']
    names = {}
    for k in range(mol.GetNumAtoms()):
        at = doc.atoms[order[k]]
        names[k] = (at.label or at.el) + str(k + 1)
    if undefined:
        out.append(f'{len(undefined)} stereocentre{"s" if len(undefined) > 1 else ""} with no wedge or hash: '
                   + ', '.join(names[k] for k in undefined[:8]) + ('…' if len(undefined) > 8 else '') +
                   ' (draw one, or leave it if the configuration is unknown or racemic).')
    chiral = {k for k, _c in centres}
    idx = {a: k for k, a in enumerate(order)}
    wedges = {}
    for b, bd in doc.bonds.items():
        if bd.stereo and bd.a in idx and bd.b in idx:
            wedges.setdefault(bd.a, []).append(bd.stereo)
            k = idx[bd.a]
            if k not in chiral:
                atom = mol.GetAtomWithIdx(k)
                if not (atom.GetDegree() >= 3 and any(bb.GetBondType() == Chem.BondType.DOUBLE for bb in atom.GetBonds())):
                    out.append(f'A {bd.stereo} bond starts at {names[k]}, which is not a stereocentre: it means '
                               'nothing there (start it at the stereocentre).')
    for a, kinds in wedges.items():
        if a in idx and (kinds.count('wedge') > 2 or kinds.count('hash') > 2):
            out.append(f'{names[idx[a]]} has {len(kinds)} wedged / hashed bonds the same way: its configuration is '
                       'ambiguous.')
    return out
