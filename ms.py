"""Mass spectrometry helpers: breaking a bond into fragment ions (as in EI / CID fragmentation schemes) and matching a
measured m/z to the structures drawn or to elemental compositions.

Fragments are worked out by electron bookkeeping on the two pieces left when a bond breaks:
  - a radical cation (an atom carrying both • and +, as in M•+): the bond's electrons go one to each side and pair
    with the radical site, so one piece keeps the charge (an even-electron cation) and the other leaves as a radical -
    "A+ + B•" or "A• + B+" (alpha cleavage of acetone M•+: CH3C≡O+ m/z 43 + CH3•);
  - an even-electron cation (+ without •, as in [M+H]+): inductive cleavage - the charged site leaves neutral and the
    other piece takes the + ("A+ + B");
  - a neutral molecule: homolysis "A• + B•" or heterolysis "A+ + B−" / "A− + B+".
A bond in a ring leaves one piece, so it is not offered.
"""
import chem
from model import Doc


def sides(doc, b):
    """The atoms on each side of bond b ([set, set]) or None when it is in a ring."""
    bd = doc.bonds[b]
    seen, todo = {bd.a}, [bd.a]
    while todo:
        x = todo.pop()
        for n, bb in doc.neighbors(x):
            if bb != b and n not in seen:
                seen.add(n); todo.append(n)
    if bd.b in seen:
        return None
    other, todo = {bd.b}, [bd.b]
    while todo:
        x = todo.pop()
        for n, bb in doc.neighbors(x):
            if bb != b and n not in other:
                other.add(n); todo.append(n)
    return [seen, other]


def options(doc, b):
    """[(spec, description)] for breaking bond b: spec is (mark on bd.a's side, mark on bd.b's side, cleared atoms),
    a mark being '+', '-', '•' or ''. Empty when the bond is in a ring."""
    parts = sides(doc, b)
    if parts is None:
        return []
    atoms = parts[0] | parts[1]
    radical_sites = [a for a in atoms if doc.atoms[a].radical and doc.atoms[a].charge]
    charge = sum(doc.atoms[a].charge for a in atoms)
    unpaired = sum(doc.atoms[a].radical for a in atoms)
    out = []
    if unpaired % 2 and charge:                             # an odd-electron ion such as M•+
        sign = '+' if charge > 0 else '-'
        site = radical_sites[:1]
        out.append(((sign, '•', site), f'charge stays left: A{sign} + B•'))
        out.append((('•', sign, site), f'charge stays right: A• + B{sign}'))
        return out
    if charge:                                              # even-electron ion, e.g. [M+H]+: inductive cleavage
        sign = '+' if charge > 0 else '-'
        for k, (mine, theirs) in enumerate(((parts[0], parts[1]), (parts[1], parts[0]))):
            sites = [a for a in theirs if doc.atoms[a].charge and (doc.atoms[a].charge > 0) == (charge > 0)]
            if sites:
                spec = (sign, '', sites[:1]) if k == 0 else ('', sign, sites[:1])
                out.append((spec, 'inductive: ' + ('A' + sign + ' + B (neutral)' if k == 0 else 'A (neutral) + B' + sign)))
        out.append((('•', '•', []), 'homolysis: A• + B•'))
        return out
    out.append((('•', '•', []), 'homolysis: A• + B•'))
    out.append((('+', '-', []), 'heterolysis: A+ + B−'))
    out.append((('-', '+', []), 'heterolysis: A− + B+'))
    return out


def pieces(doc, b, spec):
    """The two fragments of breaking bond b as new Docs (atom ids kept): [(Doc, cut atom)] for bd.a's and bd.b's side."""
    parts = sides(doc, b)
    bd = doc.bonds[b]
    mark_a, mark_b, cleared = spec
    out = []
    for side, cut, mark in ((parts[0], bd.a, mark_a), (parts[1], bd.b, mark_b)):
        d = Doc()
        d.style = dict(doc.style)
        d.insert_dict(doc.to_dict([('atom', a) for a in side]), keep_ids=True)
        for a in cleared:                                   # the old ion site gives up its charge and radical
            if a in d.atoms:
                d.atoms[a].charge = 0
                d.atoms[a].radical = 0
        at = d.atoms[cut]
        if mark == '•':
            at.radical = min(at.radical + 1, 2)
        elif mark == '+':
            at.charge += 1
        elif mark == '-':
            at.charge -= 1
        out.append((d, cut))
    return out


def describe(doc, atom_ids=None, sep='  '):
    """'C2H3O+  m/z 43.0178' for an ion, 'CH3•  15.0235 (neutral)' for a neutral piece (sep between the two)."""
    info = chem.analyse(doc, atom_ids)
    if not info.atoms:
        return ''
    dots = '•' * (info.radicals % 3)
    if info.charge:
        z = abs(info.charge)
        mz = (info.exact - info.charge * chem.ELECTRON) / z
        return f'{info.formula}{dots}{chem.charge_text(info.charge)}{sep}m/z {mz:.4f}'
    return f'{info.formula}{dots}{sep}{info.exact:.4f} (neutral)'


# ------------------------------------------------------------------------------------------- matching a measured m/z
def ppm(measured, calc):
    return (measured - calc) / calc * 1e6


def drawn_matches(doc, units, mz, tol_ppm):
    """[(where, ion label, formula, calculated m/z, ppm, atom ids)] of the drawn structures' ions (adducts included,
    charged structures as themselves) within tol_ppm of mz. units: [(name, atom ids)]."""
    out = []
    for name, ids in units:
        info = chem.analyse(doc, ids)
        if not info.atoms:
            continue
        for ion in info.ions:
            d = ppm(mz, ion.mz)
            if abs(d) <= tol_ppm:
                out.append((name, ion.label, chem.formula_text(ion.counts), ion.mz, d, ids))
    out.sort(key=lambda r: abs(r[4]))
    return out


DEFAULT_LIMITS = {'C': 80, 'H': 160, 'N': 10, 'O': 25, 'S': 3}
MAX_LIMITS = {'C': 120, 'H': 240, 'N': 20, 'O': 40, 'S': 6, 'P': 4, 'F': 20, 'Cl': 6, 'Br': 4, 'I': 3, 'Si': 6, 'B': 4,
              'Na': 1, 'K': 1}


def limits_from_text(text):
    """'C H N O S Cl' or 'C0-40 N0-4 O ...' -> {element: most}; unknown symbols are skipped."""
    import re
    out = {}
    for tok in text.replace(',', ' ').split():
        m = re.fullmatch(r'([A-Z][a-z]?)(?:(\d+)\s*-\s*)?(\d+)?', tok)
        if not m or m.group(1) not in MAX_LIMITS:
            continue
        out[m.group(1)] = int(m.group(3)) if m.group(3) else DEFAULT_LIMITS.get(m.group(1), MAX_LIMITS[m.group(1)])
    if 'H' not in out and out:
        out['H'] = MAX_LIMITS['H']
    return out


def compositions(mz, z=1, tol_ppm=5.0, limits=None, max_results=50):
    """Elemental compositions of an ion of charge z whose m/z is within tol_ppm of mz:
    [(formula, calculated m/z, ppm, RDBE, 'even' / 'odd' electron)], best first. RDBE (rings plus double bonds) of
    the ion must be at least -0.5; H follows from the other elements."""
    limits = dict(limits or DEFAULT_LIMITS)
    mass = {e: chem.PT.GetMostCommonIsotopeMass(e) for e in limits}
    target = mz * abs(z) + z * chem.ELECTRON                 # the neutral atoms' mass
    tol = target * tol_ppm * 1e-6
    heavy = [e for e in limits if e != 'H']
    heavy.sort(key=lambda e: -mass[e])
    hmax = limits.get('H', 0)
    out = []

    def valence(e):
        return {'C': 4, 'Si': 4, 'N': 3, 'P': 3, 'B': 3, 'O': 2, 'S': 2}.get(e, 1)

    def rec(k, m, counts):
        if m > target + tol:
            return
        if k == len(heavy):
            rest = target - m
            h = round(rest / mass['H']) if 'H' in mass else 0
            if h < 0 or h > hmax:
                return
            total = m + h * mass.get('H', 0)
            if abs(total - target) > tol:
                return
            c = dict(counts)
            if h:
                c['H'] = h
            if not c:
                return
            # RDBE = 1 + sum(n (v - 2)) / 2 over all atoms (H and halogens: v = 1)
            rdbe = 1 + sum(n * (valence(e) - 2) for e, n in c.items()) / 2
            if rdbe < -0.5 or (rdbe > 0.5 * sum(c.values()) + 1):
                return
            nc = c.get('C', 0)
            if nc and not (0.2 <= h / nc <= 3.1 and c.get('N', 0) / nc <= 1.3 and c.get('O', 0) / nc <= 1.2
                           and c.get('S', 0) / nc <= 0.8 and c.get('P', 0) / nc <= 0.3):
                return                                      # outside the usual ratios of real compounds
            if c.get('Na') and c.get('K'):
                return
            electrons = 'odd' if abs(rdbe - round(rdbe)) < 1e-9 else 'even'
            calc = (total - z * chem.ELECTRON) / abs(z)
            out.append((chem.formula_text(c), calc, ppm(mz, calc), rdbe, electrons))
            return
        e = heavy[k]
        for n in range(0, limits[e] + 1):
            mm = m + n * mass[e]
            if mm > target + tol:
                break
            if n:
                counts[e] = n
            else:
                counts.pop(e, None)
            rec(k + 1, mm, counts)
        counts.pop(e, None)
    rec(0, 0.0, {})
    out.sort(key=lambda r: (abs(r[2]), -r[3]))
    return out[:max_results]


def rdbe_of(counts):
    v = {'C': 4, 'Si': 4, 'N': 3, 'P': 3, 'B': 3, 'O': 2, 'S': 2}
    return 1 + sum(n * (v.get(e, 1) - 2) for e, n in counts.items()) / 2



# ------------------------------------------------------------------------------------------- annotating a spectrum
def candidate_ions(doc, units, positive=True, cleave=True, rules=True):
    """[(where, ion label, formula, calculated m/z, atom ids)] the drawn structures could give: their own ions and
    adducts (a drawn fragment ion as itself); with rules, the fragments of fragments.predict (ESI+ and EI, or ESI−,
    three steps, all likelihoods: where is '<name>, fragment: <rule>'); and, with cleave, each piece of breaking one
    bond outside a ring, as the piece capped with H (PH) protonated [PH+H]+, as a radical cation [PH]•+ and as the
    even-electron cation [PH−H]+ (P+); in negative mode [PH−H]−."""
    out, seen = [], set()

    def add(where, label, counts, mz, ids):
        key = (label, round(mz, 4))
        if key not in seen:
            seen.add(key)
            out.append((where, label, chem.formula_text(counts), mz, ids))
    for name, ids in units:
        info = chem.analyse(doc, ids)
        if not info.atoms:
            continue
        for ion in info.ions:
            if (ion.z > 0) == positive:
                add(name, ion.label, ion.counts, ion.mz, ids)
        if positive and not info.charge and not info.radicals:  # EI: the molecular ion M•+
            add(name, '[M]•+', info.counts, info.exact - chem.ELECTRON, ids)
        if rules and info.fragments == 1 and not info.errors:
            import fragments
            for mode in (('ESI+', 'ESI+Na', 'ESI+NH4', 'EI') if positive else ('ESI-',)):
                _pre, frs = fragments.predict(doc, set(ids), mode, 3)
                if _pre is None or (_pre.z > 0) != positive:
                    continue
                for f in frs:
                    add(f'{name}, fragment: {f.how}', f.ion, f.counts, f.mz, set(f.atoms))
        if not cleave or info.charge:
            continue
        for b, bd in doc.bonds.items():
            if bd.a not in ids or doc.atoms[bd.a].label or doc.atoms[bd.b].label:
                continue
            parts = sides(doc, b)
            if parts is None:
                continue
            for d, cut in pieces(doc, b, ('•', '•', [])):
                p = chem.analyse(d)
                if not p.atoms or not p.counts.get('C'):
                    continue
                ph = dict(p.counts)
                ph['H'] = ph.get('H', 0) + 1
                m_ph = p.exact + chem.H
                piece = set(d.atoms)
                f = chem.formula_text(ph)
                if positive:
                    plus = dict(ph, H=ph['H'] + 1)
                    add(f'{name}, piece', f'[{f}+H]+', plus, m_ph + chem.H - chem.ELECTRON, piece)
                    add(f'{name}, piece', f'[{f}]•+', ph, m_ph - chem.ELECTRON, piece)
                    if ph['H'] >= 1:
                        minus = dict(ph, H=ph['H'] - 1)
                        add(f'{name}, piece', chem.formula_text(minus) + '+', minus, m_ph - chem.H - chem.ELECTRON,
                            piece)
                elif ph['H'] >= 1:
                    minus = dict(ph, H=ph['H'] - 1)
                    add(f'{name}, piece', f'[{f}−H]−', minus, m_ph - chem.H + chem.ELECTRON, piece)
    return out


def annotate(peaks, candidates, tol, unit='ppm'):
    """For each peak (m/z, intensity): the best candidate within tol (ppm or Da) or None:
    [(m/z, intensity, (where, label, formula, m/z calc., ids) or None, error)]. At the same error a drawn
    structure's own ion wins over a rule-based fragment, and that over a plain cleavage piece."""
    def rank(c):
        return 2 if ', piece' in c[0] else 1 if ', fragment' in c[0] else 0
    out = []
    for mz, y in peaks:
        best, err = None, None
        for c in candidates:
            e = ppm(mz, c[3]) if unit == 'ppm' else mz - c[3]
            if abs(e) <= tol:
                score = (round(abs(e), 6), rank(c))
                if best is None or score < (round(abs(err), 6), rank(best)):
                    best, err = c, e
        out.append((mz, y, best, err))
    return out


def all_matches(peaks, candidates, tol, unit='ppm'):
    """{m/z: [(candidate, error)]} - every candidate within tol of each peak, closest first."""
    out = {}
    for mz, _y in peaks:
        hits = []
        for c in candidates:
            e = ppm(mz, c[3]) if unit == 'ppm' else mz - c[3]
            if abs(e) <= tol:
                hits.append((c, e))
        hits.sort(key=lambda h: abs(h[1]))
        out[mz] = hits
    return out



# ------------------------------------------------------------------------------------------- isotopes
def parse_formula(text):
    """{element: n} of a plain formula such as C7H11O6 (isotope labels and charges ignored)."""
    import re
    out = {}
    for el, n in re.findall(r'([A-Z][a-z]?)(\d*)', text or ''):
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def ion_charge(label):
    """The charge of an ion label: '[M+2H]2+' → 2, 'C7H11O6−' → -1, '[M]•+' → 1 (0 when it shows none)."""
    import re
    text = (label or '').replace('•', '')
    m = re.search(r'\](\d*)([+−-])$', text)                 # [M+2H]2+: the count after the bracket
    if m:
        n = int(m.group(1)) if m.group(1) else 1
    else:
        m = re.search(r'([+−-])$', text)                    # C7H11O6−: a formula's digits are not the charge
        if not m:
            return 0
        n = 1
    return n if m.group(m.lastindex) == '+' else -n


def isotopes(rows, tol, unit='ppm'):
    """Peaks left unmatched that are the M+1, M+2 ... isotope peaks of a matched ion (their m/z within tol of the
    calculated isotope peak, their height within a factor of 3 of the expected share) become matched as
    ('isotope of <ion> (<where>)', '<ion> M+k', formula, m/z calc., atom ids). Returns (rows, {monoisotopic m/z:
    'M+1 9.8 % (calc. 10.2), M+2 ...'}) - how well each matched ion's measured isotope peaks fit its formula."""
    out = list(rows)
    index = {r[0]: k for k, r in enumerate(out)}
    fits = {}
    for mz, y, match, _err in rows:
        if not match or match[0].startswith('isotope of'):
            continue
        counts = parse_formula(match[2])
        z = abs(ion_charge(match[1])) or 1
        if not counts:
            continue
        pattern = chem.isotope_pattern(counts, z, min_percent=0.5)
        if len(pattern) < 2:
            continue
        _mono_mz, mono_pct = pattern[0]
        bits = []
        for k, (iso_mz, pct) in enumerate(pattern[1:], 1):
            want = pct / mono_pct
            if want < 0.005:
                continue
            best = None
            for pm, py, pmatch, _e in out:
                if pm == mz:
                    continue
                e = ppm(pm, iso_mz) if unit == 'ppm' else pm - iso_mz
                if abs(e) <= max(tol, 3.0) if unit == 'ppm' else abs(e) <= tol:
                    if best is None or abs(e) < abs(best[1]):
                        best = (pm, e, py, pmatch)
            if best is None:
                if want >= 0.3:
                    bits.append(f'M+{k} missing (calc. {100 * want:.0f} %)')
                continue
            pm, e, py, pmatch = best
            got = py / y if y else 0.0
            bits.append(f'M+{k} {100 * got:.1f} % (calc. {100 * want:.1f})')
            if pmatch is None and want / 3 <= got <= want * 3:
                out[index[pm]] = (pm, py, (f'isotope of {match[1]} ({match[0]})', f'{match[1]} M+{k}', match[2],
                                           iso_mz, match[4]), e)
        if bits:
            fits[mz] = ', '.join(bits)
    return out, fits


# ------------------------------------------------------------------------------------------- ranking structures
def rank_structures(doc, units, peaks, positive, tol, unit='ppm', cleave=False):
    """How well each drawn structure explains a measured spectrum, best first: [{'name', 'ids', 'explained' (share
    of the measured intensity matched by its ions, fragments and their isotope peaks), 'matched', 'peaks',
    'similarity' (cosine of measured and predicted sticks, a predicted stick's height being the fragment's
    likelihood), 'score'}]. Peaks: [(m/z, intensity)]."""
    import math
    import fragments
    total = sum(y for _m, y in peaks) or 1.0
    out = []
    for name, ids in units:
        cands = candidate_ions(doc, [(name, ids)], positive, cleave)
        rows, _fits = isotopes(annotate(peaks, cands, tol, unit), tol, unit)
        explained = sum(y for _m, y, match, _e in rows if match)
        matched = sum(1 for r in rows if r[2])
        predicted = []                                      # (m/z, weight) of the ions it should give
        info = chem.analyse(doc, ids)
        for ion in info.ions:
            if (ion.z > 0) == positive:
                predicted.append((ion.mz, 0.5))
        for mode in (('ESI+', 'EI') if positive else ('ESI-',)):
            try:
                _pre, frs = fragments.predict(doc, set(ids), mode, 2)
            except Exception:
                frs = []
            predicted += [(f.mz, f.score) for f in frs]
        a, b, used = [], [], set()
        for mz, w in predicted:
            hit = None
            for k, (pm, _py) in enumerate(peaks):
                e = ppm(pm, mz) if unit == 'ppm' else pm - mz
                if abs(e) <= tol:
                    hit = k
                    break
            if hit is None and w < 0.4:                     # an unlikely ion that is absent says little
                continue
            a.append(w)
            b.append(math.sqrt(peaks[hit][1] / total) if hit is not None else 0.0)
            if hit is not None:
                used.add(hit)
        for k, (_pm, py) in enumerate(peaks):
            if k not in used and py / total > 0.02:
                a.append(0.0)
                b.append(math.sqrt(py / total))
        dot = sum(x * y for x, y in zip(a, b))
        na, nb = math.sqrt(sum(x * x for x in a)), math.sqrt(sum(y * y for y in b))
        cos = dot / (na * nb) if na and nb else 0.0
        share = explained / total
        out.append({'name': name, 'ids': ids, 'explained': share, 'matched': matched, 'peaks': len(peaks),
                    'similarity': cos, 'score': 0.6 * share + 0.4 * cos})
    out.sort(key=lambda r: -r['score'])
    return out


# ------------------------------------------------------------------------------------------- LC-MS runs
ADDUCT_MODE = {'[M+H]+': 'ESI+', '[M+Na]+': 'ESI+Na', '[M+NH4]+': 'ESI+NH4', '[M-H]-': 'ESI-'}


def search_run(doc, ids, spectra, positive, tol, unit='ppm', min_rel=0.01):
    """Where a drawn structure shows up in an LC-MS run (a list of msfile.Spectrum): the MS/MS scans whose
    precursor is one of its ions, scored by the share of their intensity its predicted fragments explain, and the
    MS1 scan where each of its ions is tallest. [{'kind': 'MS2' | 'MS1', 'index', 'rt', 'precursor', 'ion',
    'explained', 'matched', 'peaks', 'intensity'}], MS/MS best first, then the MS1 apexes."""
    import fragments
    import msfile
    info = chem.analyse(doc, ids)
    ions = [(ion.label, ion.mz) for ion in info.ions if (ion.z > 0) == positive and abs(ion.z) == 1]
    if positive and not info.charge:
        ions.append(('[M]•+', info.exact - chem.ELECTRON))
    cache = {}

    def frags_for(label):
        mode = ADDUCT_MODE.get(label.replace('−', '-'), 'EI' if '•' in label else None)
        if mode is None:
            return []
        if mode not in cache:
            try:
                pre, frs = fragments.predict(doc, set(ids), mode, 3)
            except Exception:
                pre, frs = None, []
            cache[mode] = [f.mz for f in frs] + ([pre.mz] if pre else [])
        return cache[mode]

    def close(a, b):
        return abs(ppm(a, b)) <= tol if unit == 'ppm' else abs(a - b) <= tol
    out, apex = [], {}
    for k, sp in enumerate(spectra):
        if sp.level >= 2 and sp.precursor:
            hit = next(((lab, mz) for lab, mz in ions if abs(sp.precursor - mz) <= max(0.02, mz * 20e-6)), None)
            if hit is None:
                continue
            peaks = msfile.centroid(sp, min_rel)
            total = sum(y for _m, y in peaks) or 1.0
            cands = frags_for(hit[0])
            matched, explained = 0, 0.0
            for pm, py in peaks:
                if any(close(pm, c) for c in cands):
                    matched += 1
                    explained += py
            out.append({'kind': 'MS2', 'index': k, 'rt': sp.rt, 'precursor': sp.precursor, 'ion': hit[0],
                        'explained': explained / total, 'matched': matched, 'peaks': len(peaks), 'intensity': total})
        elif sp.level == 1 and len(sp.mz):
            for lab, mz in ions:
                lo, hi = (mz * (1 - tol * 1e-6), mz * (1 + tol * 1e-6)) if unit == 'ppm' else (mz - tol, mz + tol)
                sel = (sp.mz >= lo) & (sp.mz <= hi)
                if sel.any():
                    y = float(sp.intensity[sel].max())
                    if y > apex.get(lab, (0.0, None))[0]:
                        apex[lab] = (y, k)
    for lab, (y, k) in apex.items():
        out.append({'kind': 'MS1', 'index': k, 'rt': spectra[k].rt, 'precursor': None, 'ion': lab, 'explained': None,
                    'matched': None, 'peaks': None, 'intensity': y})
    out.sort(key=lambda r: (r['kind'] != 'MS2', -(r['explained'] or 0), -(r['intensity'] or 0)))
    return out
