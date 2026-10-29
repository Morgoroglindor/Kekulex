"""
pubchem.py - names through PubChem (online; only when asked from the Structure menu).

    pubchem.smiles_for_name('nitrobenzene')        # 'C1=CC=C(C=C1)[N+](=O)[O-]'
    pubchem.name_for_inchikey('LQNUZADURLCDLV-...')# ('nitrobenzene', 'nitrobenzene') = (common name, IUPAC name)

Only the name typed, or the structure's InChIKey, is sent - nothing else from the drawing. PubChem knows the compounds
in its database; a new structure has no name there.
"""
import json
import urllib.error
import urllib.parse
import urllib.request

BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'
TIMEOUT = 15


class PubChemError(Exception):
    pass


def _get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'StructureEditor'})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            return json.loads(r.read().decode('utf-8'))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise PubChemError(f'PubChem answered {e.code}') from e
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise PubChemError(f'PubChem could not be reached ({getattr(e, "reason", e)})') from e


def smiles_for_name(name):
    """The SMILES PubChem has for a name (None: not found)."""
    url = f'{BASE}/name/{urllib.parse.quote(name.strip())}/property/SMILES,IsomericSMILES/JSON'
    data = _get(url)
    if not data:
        return None
    props = data.get('PropertyTable', {}).get('Properties', [])
    if not props:
        return None
    return props[0].get('SMILES') or props[0].get('IsomericSMILES') or props[0].get('CanonicalSMILES')


def name_for_inchikey(key):
    """(title, IUPAC name) PubChem has for an InChIKey, or None."""
    data = _get(f'{BASE}/inchikey/{urllib.parse.quote(key)}/property/Title,IUPACName/JSON')
    if not data:
        return None
    props = data.get('PropertyTable', {}).get('Properties', [])
    if not props:
        return None
    return props[0].get('Title', ''), props[0].get('IUPACName', '')
