"""
    pubchem_scrape.py
    Dylan H. Ross
    2018/09/17

        Tool for using PubChem's REST API to search for compounds and obtain SMILES structures.

        more info here:
            http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
"""


def search_by_name(session, name):
    """
search_by_name
    description:
        Searches for a PubChem CID using a compound name, returning a list of results.
        Returns None if any errors.
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
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        # split response on whitespace
        return [int(_) for _ in resp.split()]
    except:
        print("Failed to retrieve PubChem CID for compound: {}".format(name), end=" ")
        return None


def cid_fetch_smiles(session, cid, canonical=True):
    """
cid_fetch_smiles
    description:
        Fetches the SMILES structure from the record with the specified PubChem
        CID, returning it as a string. Either the Canonical SMILES or the Isomeric
        SMILES may be fetched (see canonical kwarg). The Default is Canonical
        since that is more likely to be present, however, the Isomeric is preferred
        when available since the overall goal will be to produce 3D structures.
    parameters:
        session (requests.Session) -- requests session
        cid (int) -- PubChem CID
        [canonical (bool)] -- Whether to fetch the canonical SMILES (True) or the
                                isomeric SMILES (False) [optional, default=True]
    returns:
        (str) -- SMILES structure or empty string if no results
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
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        return resp.strip()
    except:
        print("Failed to retrieve {} SMILES for PubChem CID: {}".format(stype, cid), end=" ")
        return ""
