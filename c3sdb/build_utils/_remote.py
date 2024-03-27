"""
    c3sdb/build_utils/_remote.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Module with utilities for accessing remote information via web APIs

    - PubChem REST API to search for compounds and obtain SMILES structures 
        (see: http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial)
    - Lipid MAPS REST API to search the LMSD for lipid SMILES structures
"""


from typing import List, Optional, Dict, Any
from time import sleep

import requests


# delay (in seconds) before sending any request 
_REQUEST_DELAY: float = 0.25


def pubchem_search_by_name(session: requests.Session, 
                           name: str
                           ) -> Optional[List[int]] :
    """
    Searches for a PubChem CID using a compound name, returning a list of results,
    returns None if any errors

    Parameters
    ----------
    session : ``requests.Session``
        requests session
    name : ``str``
        name of the compound to search
    
    Returns
    -------
    pubchem_cids
        parameters:
            session (requests.Session) -- requests session
            name (str) -- name of the compound to search
        returns:
            (list(int) or None) -- list of PubChem CID(s) matching the search name or None if unsuccessful
    """
    # construct the request URL per REST API documentation
    url_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url_input = "compound/name/"
    url_command = "/cids/"
    url_output = "TXT"
    url = url_prolog + url_input + name + url_command + url_output
    try:
        sleep(_REQUEST_DELAY)
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        # split response on whitespace
        return [int(_) for _ in resp.split()]
    except:
        print("Failed to retrieve PubChem CID for compound: {}".format(name), end=" ")
        return None


def pubchem_cid_fetch_smiles(session: requests.Session, 
                             cid: int, 
                             canonical: bool = True
                             ) -> Optional[str] :
    """
    Fetche the SMILES structure from the record with the specified PubChem CID

    Either the Canonical SMILES or the Isomeric SMILES may be fetched, cannonical
    is default since that is more likely to be present but isomeric is preferred
    when available since the overall goal will be to produce 3D structures

    Parameters
    ----------
    session : ``requests.Session``
        requests session
    cid : ``int``
        PubChem CID
    cannonical : ``bool``, default=True
        Whether to fetch the canonical SMILES (True) or the isomeric SMILES (False)

    Returns
    -------
    smiles : ``str`` or None
        SMILES structure or None if no results
    """
    # construct the request URL per REST API documentation
    url_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url_input = "compound/cid/"
    # option to search for Isomeric SMILES
    stype = "Canonical"
    if not canonical:
        stype = "Isomeric"
    url_command = "/property/{}SMILES/".format(stype)
    url_output = "TXT"
    url = url_prolog + url_input + str(cid) + url_command + url_output
    try:
        sleep(_REQUEST_DELAY)
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        return resp.strip()
    except:
        print("Failed to retrieve {} SMILES for PubChem CID: {}".format(stype, cid), end=" ")
        return None


def _str_from_lipid_dict(lipid: Dict[Any], 
                         ign_fa_comp: bool
                         ) -> str :
    """
    Converts a lipid definition in the form of a dictionary reflecting the structure of the lipids.db
    schema into a single string with standardized formatting such that it may be recognized by the
    Lipid MAPS structure database REST API.
    
    Parameters
    ----------
    lipid : ``dict(...)``
        dict defining target lipid with structure reflecting that of the lipids.db schema
        TODO (Dylan Ross): What is this referring to?
    ign_fa_comp : ``bool``
        flag indicating whether to ignore individual fatty acid composition and just construct
        the name using the sum composition

    Returns
    -------
    lipid_name : ``str``
        lipid name in standardized format
    """
    lipid_str = "{}(".format(lipid["lipid_class"])
    if "fa_mod" in lipid:
        # add in fatty acid modifier 
        lipid_str += "{}-".format(lipid["fa_mod"].capitalize())
    if not ign_fa_comp and "fa1_n_carbon" in lipid: 
        # deal with individual chains
        lipid_str += "{}:{}_{}:{}".format(lipid["fa1_n_carbon"], lipid["fa1_n_unsat"], 
                                          lipid["fa2_n_carbon"], lipid["fa2_n_unsat"])
        if "fa3_n_carbon" in lipid:
            lipid_str += "_{}:{}".format(lipid["fa3_n_carbon"], lipid["fa3_n_unsat"])
    else:
        # total fatty acid composition
        lipid_str += "{}:{}".format(lipid["n_carbon"], lipid["n_unsat"])
    lipid_str += ")"
    return lipid_str
    

def lmaps_fetch_smiles(session: requests.Session, 
                       lipid: Dict[Any]
                       ) -> Optional[str] :
    """
    Uses the Lipid MAPS REST API to search the LMSD for lipid SMILES structures. 
    
    Takes as input a lipid in the form of a dictionary with a structure reflecting 
    that of lipids.db schema. If individual fatty acid composition is specified, then 
    try searching using the 'abbrev_chains' input item first, then the 'abbrev' input 
    item if that fails.

    Parameters
    ----------
    session : ``requests.Session``
        requests session
    lipid : ``dict(...)``
        dict defining target lipid with structure reflecting that of the lipids.db schema
        TODO (Dylan Ross): What is this referring to?
    
    Returns
    -------
    smiles : ``str`` or ``None``
        SMILES structure or None in case of error
    """
    url = "https://lipidmaps.org/rest/compound/abbrev{}/{}/smiles/"
    lipid_name = _str_from_lipid_dict(lipid, False)
    result = None
    if "_" in lipid_name or "/" in lipid_name:
        # includes individual fatty acid composition
        try:
            sleep(_REQUEST_DELAY)
            result = session.get(url.format("_chains", lipid_name)).json()
        except:
            msg = "fetch_lipid_smiles: search with 'abbrev_chains' for lipid {} failed".format(lipid_name)
            print(msg, end=" ")
        # try again with 
        if not result:
            try:
                sleep(_REQUEST_DELAY)
                result = session.get(url.format("", lipid_name)).json()
            except Exception as e:
                print(e)
                msg = "fetch_lipid_smiles: search with 'abbrev' for lipid {} failed".format(lipid_name)
                print(msg, end=" ")
        # regenerate the name but with total fatty acid composition in case the above did not work
        lipid_name = _str_from_lipid_dict(lipid, True)
    if not result:
        try:
            sleep(_REQUEST_DELAY)
            result = session.get(url.format("", lipid_name)).json()
        except:
            msg = "fetch_lipid_smiles: search with 'abbrev' for lipid {} failed".format(lipid_name)
            print(msg, end=" ")
    # if there was a result, just return the first one
    smi = None
    if result:
        if 'Row1' in result:
            smi = result['Row1']['smiles']
        elif 'smiles' in result:
            smi = result['smiles']
    return smi

