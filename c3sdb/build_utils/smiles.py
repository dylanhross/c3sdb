"""
    c3sdb/build_utils/smiles.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with utilities for generating SMILES structures
"""


from typing import Optional, Dict, Tuple
import json
import sqlite3
import os
import time

import requests

from c3sdb.build_utils._parsing import parse_lipid, parse_peptide
from c3sdb.build_utils._remote import (
    pubchem_cid_fetch_smiles, pubchem_search_by_name, lmaps_fetch_smiles
)


# path to SMILES search cache file (embedded in the package)
_SMILES_SEARCH_CACHE: str = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                                         "_include/smiles_search_cache.json")


# SMILES strings for each amino acid (excludes C-terminal OH)
_AA_BASE: str = 'N[C@@H]({})C(=O)'
_AA_TO_SMI: Dict[str, str] = {
    'A': _AA_BASE.format('C'),
    'C': _AA_BASE.format('CS'),
    'D': _AA_BASE.format('CC(=O)O'),
    'E': _AA_BASE.format('CCC(=O)O'),
    'F': _AA_BASE.format('CC{0}=CC=CC=C{0}'),
    'G': 'NCC(=O)',
    'H': _AA_BASE.format('CC{0}=CNC=N{0}'),
    'I': _AA_BASE.format('C(C)CC'),
    'K': _AA_BASE.format('CCCCN'),
    'L': _AA_BASE.format('CC(C)C'),
    'M': _AA_BASE.format('CCSC'),
    'N': _AA_BASE.format('CC(=O)N'),
    'P': 'N{0}C[C@@H](CCC{0})C(=O)',
    'Q': _AA_BASE.format('CCC(=O)N'),
    'R': _AA_BASE.format('CCCNC(=N)N'),
    'S': _AA_BASE.format('CO'),
    'T': _AA_BASE.format('C[C@@H](C)O'),
    'V': _AA_BASE.format('C(C)C'),
    'W': _AA_BASE.format('CC{0}=CNC{1}=C(C=CC=C{1}){0}'),
    'Y': _AA_BASE.format('CC{0}=CC=C(O)C=C{0}')
}


def _peptide_seq_to_smiles(seq: str
                           ) -> str :
    """
    generates a SMILES structure for a linear peptide sequence
    
    Parameters
    ----------
    seq : ``str``
        peptide sequence
    
    Returns
    -------
    smi : ``str``
        peptide SMILES
    """ 
    smi = ''
    rings = 0
    for aa in seq:
        aa_smi = _AA_TO_SMI[aa]
        if '{1}' in aa_smi:
            rings += 1
            r0 = rings
            rings += 1
            r1 = rings
            smi += aa_smi.format(r0, r1)
        elif '{0}' in aa_smi:
            rings += 1
            r0 = rings
            smi += aa_smi.format(r0)
        else:
            smi += aa_smi
    smi += 'O'
    return smi


def _carbon_chain(c: int, 
                  u: int
                  ) -> str :
    """
    Returns a SMILES representation of a carbon chain with c carbons and u unsaturations
    
    Parameters
    ----------
    c : ``int``
        number of carbons
    u : ``int``
        number of unsaturations

    Returns
    -------
    smi : ``str`` 
        SMILES representation of carbon chain
    """
    smi = "C"
    c -= 1
    while u > 0:
        smi += "C=CC"
        c -= 3
        u -= 1
    while c > 0:
        smi += "C"
        c -= 1
    return smi


def _sphingo_smiles(c: str, 
                    n_carbon: int, 
                    n_unsat: int, 
                    fa_mod: str
                    ) -> Optional[str] :
    """
    Returns a SMILES structure (str) for a sphingolipid using a base SMILES structure that 
    includes head group

    Parameters
    ----------
    c : ``str``
        lipid class (abbreviated)
    n_carbon : ``int``
        number of FA carbons (sum composition)
    n_unsat : ``int``
        number of FA unsaturations (sum composition)
    fa_mod : ``str``
        FA modifier (should be "d" for sphingolipids)

    Returns
    -------
    smi : ``str`` or ``None`` 
        SMILES structure, None if anything goes wrong
    """
    # check fatty acid modifier (should only be 'd' for this class)
    if fa_mod and fa_mod != 'd':
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'SM': '[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(NC({})=O)[C@]([H])(O){}',
        'Cer': '[C@](CO)([H])(NC({})=O)[C@]([H])(O){}',
        'GlcCer': '[C@](CO[C@@H]1O[C@H](CO)[C@@H](O)C(O)C1O)([H])(NC({})=O)[C@]([H])(O){}',
        'HexCer': 'C(COC1OC(CO)C(O)C(O)C1O)([H])(NC({})=O)[C@]([H])(O){}'
    }
    if c not in base_smi_c:
        # the specifiec class is not implemented for this sphingolipid smiles generator
        # no error, just return None
        return None
    base_smi = base_smi_c[c]
    if n_unsat == 0:
        return base_smi.format(_carbon_chain(n_carbon - 19, 0), _carbon_chain(15, 0))
    elif n_unsat == 1:
        return base_smi.format(_carbon_chain(n_carbon - 19, 0), _carbon_chain(15, 1))
    else:
        return base_smi.format(_carbon_chain(n_carbon - 19, n_unsat - 1), _carbon_chain(15, 1))


def _phospho_smiles(c: str, 
                    n_carbon: int, 
                    n_unsat: int, 
                    fa_mod: Optional[str]
                    ) -> Optional[str] :
    """
    Returns a SMILES structure (str) for a phospholipid using a base SMILES structure which
    includes head group

    Parameters
    ----------
    c : ``str``
        lipid class
    n_carbon : ``int``
        fatty acid carbons
    n_unsat : ``int``
        fatty acid unsaturation
    fa_mod : ``str``
        fatty acid modifier 
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES structure, None if anything goes wrong
    """
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'PC': 'C[N+](C)(C)CCOP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PE': 'NCCOP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PS': 'C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PA': 'OP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O',
        'PG': 'OCC(O)COP(OC[C@]([H])(OC({})=O)COC({})=O)(=O)O'
    }
    if c not in base_smi_c:
        # this specific class not implemented
        # no error, return None
        return None
    base_smi = base_smi_c[c]
    # split the fatty acid carbons and unsaturations evenly between the two fatty acids
    nc, nu = n_carbon - 2, n_unsat
    nc_a, nc_b = nc / 2, nc / 2 + nc % 2
    nu_a, nu_b = nu / 2, nu / 2 + nu % 2 
    return base_smi.format(_carbon_chain(nc_a, nu_a), _carbon_chain(nc_b, nu_b))


def _lysophospho_smiles(c: str, 
                        n_carbon: int, 
                        n_unsat: int, 
                        fa_mod: Optional[str]
                        ) -> Optional[str] :
    """
    Returns a SMILES structure (str) for a lysophospholipid using a base SMILES structure which
    includes head group

    .. note::
        the generated SMILES structure corresponds to the 2-lysophospholipid
    
    Parameters
    ----------
    c : ``str``
        lipid class
    n_carbon : ``int``
        fatty acid carbons
    n_unsat : ``int``
        fatty acid unsaturation
    fa_mod : ``str``
        fatty acid modifier 
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES structure, None if anything goes wrong
    """
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'LPC': 'C[N+](C)(C)CCOP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPE': 'NCCOP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPS': 'C(O)(=O)[C@@]([H])(N)COP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPA': 'OP(OC[C@]([H])(O)COC({})=O)(=O)O',
        'LPG': 'OCC(O)COP(OC[C@]([H])(O)COC({})=O)(=O)O'
    }
    if c not in base_smi_c:
        # the specifiec class is not implemented for this smiles generator
        # no error, return None
        return None
    base_smi = base_smi_c[c]
    # only one fatty acid position
    nc, nu = n_carbon - 1, n_unsat
    return base_smi.format(_carbon_chain(nc, nu))


def _glycero_smiles(c: str, 
                    n_carbon: int, 
                    n_unsat: int, 
                    fa_mod: Optional[str]
                    ) -> Optional[str] :
    """
    Returns a SMILES structure (str) for a glycerolipid using a base SMILES structure which
    includes head group

    Parameters
    ----------
    c : ``str``
        lipid class
    n_carbon : ``int``
        fatty acid carbons
    n_unsat : ``int``
        fatty acid unsaturation
    fa_mod : ``str``
        fatty acid modifier 
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES structure, None if anything goes wrong
    """
    # check fatty acid modifier 
    if fa_mod:
        # return None for unimplemented fa_mod
        return None 
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi_c = {
        'DG': 'C(O)[C@]([H])(OC({})=O)COC({})=O',
        'DGDG': 'C(O[C@@H]1O[C@H](CO[C@H]2O[C@H](CO)[C@H](O)C(O)C2O)[C@H](O)C(O)C1O)[C@]([H])(OC({})=O)COC({})=O'
    }
    if c not in base_smi_c:
        # the specifiec class is not implemented for this smiles generator
        # no error, return None
        return None
    base_smi = base_smi_c[c]
    # split the fatty acid carbons and unsaturations evenly between the two fatty acids
    nc, nu = n_carbon - 2, n_unsat
    nc_a, nc_b = nc / 2, nc / 2 + nc % 2
    nu_a, nu_b = nu / 2, nu / 2 + nu % 2 
    return base_smi.format(_carbon_chain(nc_a, nu_a), _carbon_chain(nc_b, nu_b))


def _tg_smiles(n_carbon: int, 
               n_unsat: int, 
               fa_mod: str
               ) -> Optional[str] :
    """
    Returns a SMILES structure (str) for a TG using a base SMILES structure which 
    includes head group

    Parameters
    ----------
    c : ``str``
        lipid class
    n_carbon : ``int``
        fatty acid carbons
    n_unsat : ``int``
        fatty acid unsaturation
    fa_mod : ``str``
        fatty acid modifier
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES structure, None if anything goes wrong
    """
    # check fatty acid modifier
    if fa_mod:
        # return None for unimplemented fa_mod
        return None
    # base SMILES structure (includes head group and positions for fatty acids)
    base_smi = '{}C(=O)OCC(COC(=O){})OC(=O){}'
    # split the fatty acid carbons and unsaturations evenly between the three fatty acids
    nc, nu = n_carbon - 3, n_unsat
    # NC needs to be at least 3 (so at least 6 in input)
    if nc < 3:
        return None
    nc_a = nc // 3 + nc % 3
    nc_b = nc // 3
    nc_c = nc // 3
    nu_a = nu // 3 + nu % 3
    nu_b = nu // 3
    nu_c = nu // 3
    return base_smi.format(_carbon_chain(nc_a, nu_a), _carbon_chain(nc_b, nu_b), _carbon_chain(nc_c, nu_c))


def _generate_lipid_smiles(lipid_cls: str, 
                           n_carbon: int, 
                           n_unsat: int, 
                           fa_mod: Optional[str] = None
                           ) -> Optional[str]:
    """
    Generates a SMILES structure given a lipid class, sum fatty acid composition, and
    optional fatty acid modifier. Returns None if there are any errors.
    
    Parameters
    ----------
    lipid_cls : ``str``
        lipid class
    n_carbon : ``int``
        fatty acid carbons
    n_unsat : ``int``
        fatty acid unsaturation
    fa_mod : ``str``, default=None
        fatty acid modifier [optional, default=None]
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES structures (or None if any errors)
    """
    sphingo_cls = ['SM', 'Cer', 'GlcCer', 'HexCer']
    phospho_cls = ['PC', 'PE', 'PS', 'PA', 'PG']
    lysophospho_cls = ['LPC', 'LPE', 'LPS', 'LPA', 'LPG']
    glycero_cls = ['DG', 'DGDG']
    smi = None
    # sphingolipids
    if lipid_cls in sphingo_cls:
        smi = _sphingo_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)
    # phospholipids
    elif lipid_cls in phospho_cls:
        smi = _phospho_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)
    # lysophospholipids
    elif lipid_cls in lysophospho_cls:
        smi = _lysophospho_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)
    # glycerolipids
    elif lipid_cls in glycero_cls:
        smi = _glycero_smiles(lipid_cls, n_carbon, n_unsat, fa_mod)
    # TG
    elif lipid_cls == 'TG':
        smi = _tg_smiles(n_carbon, n_unsat, fa_mod)
    return smi


def save_smiles_search_cache(smiles_search_cache: Dict[str, str], 
                             cache_file_name: str
                             ) -> None :
    """
    save the search cache to file in .json format. Overwrites any previous existing cache
    
    Parameters
    ----------
    smiles_search_cache : ``dict(str:str)``
        search cache mapping compound names to SMILES structures
    cache_file_name : ``str``
        path to save the search cache file into
    """
    with open(cache_file_name, "w") as j:
        json.dump(smiles_search_cache, j, indent=4)


def load_smiles_search_cache(cache_file_name: Optional[str] = None
                             ) -> Dict[str, str] :
    """
    open the search cache if one exists, return it as a dictionary
    
    Parameters
    ----------
    cache_file : ``str`` or ``None``, default=None
        path to load the search cache from, if None load the 
        search cache file embedded in this package
    
    Returns
    -------
    smiles_search_cache : ``dict(str:str)``
        search cache mapping compound names to SMILES structures, return empty dict if
        search cache file not found
    """
    cache_file_name = _SMILES_SEARCH_CACHE if cache_file_name is None else cache_file_name
    if os.path.isfile(cache_file_name):
        with open(cache_file_name, "r") as jf:
            return json.load(jf)
    else:
        return {}


def _smi_from_cid(session: requests.Session, 
                  cid: int
                  ) -> Tuple[Optional[str], int]:
    """
    a wrapper around pubchem_cid_fetch_smiles that tries to get the isomeric SMILES
    first, and failing that it tries to get the canonical SMILES. Returns None if both of
    those fail. 

    Parameters
    ----------
    session : ``requests.Session``
        requests session for doing web requests
    cid : ``str``
        PubChem CID
    
    Returns
    -------
    smi : ``str`` or ``None``
        SMILES string or None if unsuccessful
    """
    n_requests = 0
    smiles, n_req = pubchem_cid_fetch_smiles(session, cid, canonical=False)
    n_requests += n_req
    if not smiles:
        smiles, n_req = pubchem_cid_fetch_smiles(session, cid)
        n_requests += n_req
    return smiles, n_requests


def add_smiles_to_db(cursor: sqlite3.Cursor, 
                     session: requests.Session, 
                     smiles_search_cache: Dict[str, str], 
                     gen_lipid_smi: bool = True
                     ) -> Tuple[int, int] :
    """
    Fetches SMILES structures for the entries from the C3S.db using compound names

    workflow:
    - check the search cache and see if there is an entry for the compound name
    - if the name can be parsed as a lipid 
        - try to scrape LIPID MAPS for the SMILES structure 
        - or else just use the lipid SMILES generator 
    - try to search PubChem by compound name to get a CID then use that to retrieve a SMILES 

    Parameters
    ----------
    cursor : ``sqlite3.cursor``
        cursor for running queries against the drugs.db database
    session : ``requests.Session`` 
        requests session for doing web requests
    search_cache:  ``dict(str:str)``
        search cache mapping compound names to SMILES structures
    gen_lipid_smi : ``bool``, default=True
        if a lipid name is able to be parsed but searching LIPID MAPS does not yield a
        SMILES structure, then use the lipid SMILES generator to generate a generic
        SMILES structure matching the lipid class and fatty acid composition 
    
    Returns
    -------
    n_smiles : ``int``
        number of SMILES structures added to the database
    n_requests : ``int``
        number of web requests that were sent
    """
    # track the number of web requests that were sent
    n_requests = 0
    # map global identifiers to SMILES structures
    gid_to_smi = {}
    qry_sel = "SELECT g_id, name FROM master WHERE smi IS NULL"
    n = 0
    for g_id, name in cursor.execute(qry_sel).fetchall():
        n += 1
        print(f"\r\t({n:6d}) {g_id} {name:<100s}", end="")
        if name in smiles_search_cache:
            # first check the search cache for the compound
            gid_to_smi[g_id] = smiles_search_cache[name]
        elif p_lipid := parse_lipid(name):
            #  search lipid maps
            smi, n_req = lmaps_fetch_smiles(session, p_lipid)
            n_requests += n_req
            if smi:
                gid_to_smi[g_id] = smi
                # add entry to search cache
                smiles_search_cache[name] = smi
            elif gen_lipid_smi:
                # try to generate a lipid SMILES structure
                if "fa_mod" not in p_lipid:
                    p_lipid["fa_mod"] = None
                lc, nc, nu, fm = p_lipid["lipid_class"], p_lipid["n_carbon"], p_lipid["n_unsat"], p_lipid["fa_mod"]
                smi = _generate_lipid_smiles(lc, nc, nu, fa_mod=fm)
                if smi:
                    gid_to_smi[g_id] = smi
        elif parse_peptide(name):
            # if parsed the name as a peptide, generate a SMILES structure from that
            smi = _peptide_seq_to_smiles(name)
            gid_to_smi[g_id] = smi
        else:
            # finally search pubchem by name (try to find CID first, then grab the SMILES)  
            cids, n_req = pubchem_search_by_name(session, name)
            n_requests += n_req
            # do one retry for failed queries
            if cids is None:
                cids, n_req = pubchem_search_by_name(session, name)
                n_requests += n_req
            if cids is not None:
                # if multiple CIDs were found just use the first
                cid = cids[0]
                smi, n_req = _smi_from_cid(session, cid)
                n_requests += n_req
                if smi:
                    gid_to_smi[g_id] = smi
                    # add entry to search cache
                    smiles_search_cache[name] = smi
    # add the SMILES structures to the master table by g_id for any compounds that were matched
    print("\n\tadding SMILES structures to database ...", end=" ")
    qry_updt = "UPDATE master SET smi=? WHERE g_id=?"
    for g_id in gid_to_smi:
        qdata = (gid_to_smi[g_id], g_id)
        cursor.execute(qry_updt, qdata)
    print("done")
    # return number of SMILES structures added
    return len(gid_to_smi), n_requests
