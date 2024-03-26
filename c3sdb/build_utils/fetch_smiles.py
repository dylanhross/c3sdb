#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
    fetch_smiles.py
    Dylan H. Ross
    2018/09/17
    
        
"""


import os
from json import dump as jdump, load as jload
from time import sleep


from c3sdb.build_utils.parsing import parse_carbohydrate, parse_lipid, parse_peptide


# change how we import other tools/build scripts depending upon whether this is being called directly
# or being used in the context of the db_setup module (i.e. to generate documentation)
if __name__ == '__main__':
    from lipid_parser import parse_lipid
    from peptide_parser import parse_peptide
    from lmaps_scrape import search_lmaps_smiles
    from generate_lipid_smiles import generate_lipid_smiles
    from generate_peptide_smiles import seq_to_smiles
    from pubchem_scrape import search_by_name, cid_fetch_smiles
else:
    from db_setup.lipid_parser import parse_lipid
    from db_setup.peptide_parser import parse_peptide
    from db_setup.lmaps_scrape import search_lmaps_smiles
    from db_setup.generate_lipid_smiles import generate_lipid_smiles
    from db_setup.generate_peptide_smiles import seq_to_smiles
    from db_setup.pubchem_scrape import search_by_name, cid_fetch_smiles


def save_search_cache(search_cache, cache_name="search_cache.json"):
    """
save_search_cache
    description:
        save the search cache to file in .json format. Overwrites any previous existing cache
    parameters:
        search_cache (dict(str:str)) -- search cache mapping compound names to SMILES structures
        [cache_name (str)] -- path to the search cache file location [optional, default="search_cache.json"]
"""
    with open(cache_name, "w") as j:
        jdump(search_cache, j, indent=4)


def load_search_cache(cache_name="search_cache.json"):
    """
load_search_cache
    description:
        open the search cache if one exists, return it as a dictionary
    parameters:
        [cache_name (str)] -- path to the search cache file location [optional, default="search_cache.json"]
    returns:
        (dict(str:str)) -- the saved search cache mapping compound names to SMILES structures or an empty dict if the
                            cache file could not be found
"""
    # try to open the cache, if the file does not exist then just return an empty dict
    try:
        with open(cache_name, "r") as j:
            cache = jload(j)
    except:
        cache = {}
    return cache


def smi_from_cid(sess, cid, delay=0.2):
    """
smi_from_cid
    description:
        a wrapper around pubchem_scrape.cid_fetch_smiles that tries to get the isomeric SMILES
        first, and failing that it tries to get the canonical SMILES. Returns None if both of
        those fail. Also adds time delays between requests so as not to hammer the PubChem
        servers unduly
        per guidelines here:
            https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access$_RequestVolumeLimitations
            no more than 5 requests/s or 400 requests/min -> default delay of 0.2 s is fine
    parameters:
        session (requests.Session) -- requests session for doing web requests
        cid (str) -- PubChem CID
        [delay (float)] -- time delay (in seconds) between successive requests [optional,
                            default=0.25]
    returns:
        (str or None) -- SMILES string or None if unsuccessful
"""
    sleep(delay)
    smiles = cid_fetch_smiles(sess, cid, canonical=False)
    if not smiles:
        print(" no isomeric SMILES, trying canonical ", end="")
        # double the delay on back-to-back requests
        sleep(delay * 2)
        smiles = cid_fetch_smiles(sess, cid)
    if smiles:
        return smiles
    else:
        # I want None instead of the empty string that normally
        # gets returned by cid_fetch_smiles on error
        return None


def cid_from_name(sess, name, delay=0.25):
    """
smi_from_cid
    description:
        a wrapper around pubchem_scrape.search_by_name that tries to get the PubChem CID for a
        compound by name
        Also adds time delays between requests so as not to hammer the PubChem
        servers unduly
        per guidelines here:
            https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access$_RequestVolumeLimitations
            no more than 5 requests/s or 400 requests/min -> default delay of 0.25 s is fine
    parameters:
        session (requests.Session) -- requests session for doing web requests
        name (str) -- compound name to search
        [delay (float)] -- time delay (in seconds) between successive requests [optional, default=0.25]
    returns:
        (int or None) -- PubChem CID or None if unsuccessful
"""
    sleep(delay)
    cids = search_by_name(sess, name)
    if not cids:
        return None
    # only return the first CID (if multiple were found...)
    return cids[0]


def fetch_smi_byname(cursor, session, search_cache, delay=0.25, gen_lipid_smi=True):
    """
fetch_smi_byname
    description:
        Fetches SMILES structures for the entries from the C3S.db using compound names.
        First:
            * check the search cache and see if there is an entry for the compound name
        Next, if the name can be parsed as a lipid (lipid_parser.py):
            * try to scrape LIPID MAPS for the SMILES structure (lmaps_scrape.py)
            * or else just use the lipid SMILES generator (generate_lipid_smiles.py)
        Finally:
            * try to search PubChem by compound name to get a CID then use that to retrieve a SMILES (pubchem_scrape.py)
    parameters:
        cursor (sqlite3.cursor) -- cursor for running queries against the drugs.db database
        session (requests.Session) -- requests session for doing web requests
        search_cache (dict(str:str)) -- search cache mapping compound names to SMILES structures
        [delay (float)] -- time delay (in seconds) between successive requests [optional, default=0.25]
        [gen_lipid_smi (bool)] -- if a lipid name is able to be parsed but searching LIPID MAPS does not yield a
                                    SMILES structure, then use the lipid SMILES generator to generate a generic
                                    SMILES structure matching the lipid class and fatty acid composition [optional,
                                    default=True]
"""
    # map global identifiers to SMILES structures
    gid_to_smi = {}

    qry = "SELECT g_id, name FROM master WHERE smi IS NULL"
    for g_id, name in cursor.execute(qry).fetchall():
        
        if name in search_cache:
            # first check the search cache for the compound
            print("found SMILES in search cache for compound: {} (g_id: {}) ... ".format(name, g_id), end="")
            gid_to_smi[g_id] = search_cache[name]
            print("ok")

        elif parse_lipid(name):
            # next try to parse the name as a lipid then search lipid maps
            p_lipid = parse_lipid(name)
            print("parsed lipid: {} (g_id: {}) ... searching LIPID MAPS ... ".format(name, g_id), end="")
            smi = search_lmaps_smiles(sess, p_lipid)
            if smi:
                gid_to_smi[g_id] = smi
                # add entry to search cache
                search_cache[name] = smi
                print("ok")
            elif gen_lipid_smi:
                # try to generate a generic SMILES structure
                print("generating generic SMILES ... ", end="")
                if "fa_mod" not in p_lipid:
                    p_lipid["fa_mod"] = None
                lc, nc, nu, fm = p_lipid["lipid_class"], p_lipid["n_carbon"], p_lipid["n_unsat"], p_lipid["fa_mod"]
                smi = generate_lipid_smiles(lc, nc, nu, fa_mod=fm)
                if smi:
                    gid_to_smi[g_id] = smi
                    # we DO NOT add entry to search cache for this since it is generated NOT scraped
                    #search_cache[name] = smi
                    print("ok")
                else:
                    print("FAILED TO GENERATE LIPID SMILES")
            else:
                print("NO SMILES FROM LIPID MAPS")
        
        elif parse_peptide(name):
            # try to parse the name as a peptide and generate a SMILES structure from that
            print("matched peptide regex ... ", end="")
            try:
                print("generating peptide SMILES ... ", end="")
                smi = seq_to_smiles(name)
                gid_to_smi[g_id] = smi
                # we DO NOT add entry to search cache for this since it is generated NOT scraped
                #search_cache[name] = smi
                print("ok")
            except Exception as e:
                print("FAILED TO GENERATE PEPTIDE SMILES")
        
        else:
            # finally search pubchem by name (try to find CID first, then grab the SMILES)  
            print("searching PubChem by name: '{}' (g_id: {}) ... ".format(name, g_id), end="")
            cid = cid_from_name(session, name, delay=delay)
            # wait a delay, then try again if the first time didn't work
            if not cid:
                sleep(delay)
                cid = cid_from_name(session, name, delay=delay)
                print("retrying ... ", end="")
            if cid:
                print("matched PubChem CID: {} -- fetching SMILES ... ".format(cid), end="")
                smi = smi_from_cid(session, cid, delay=delay)
                if smi:
                    gid_to_smi[g_id] = smi
                    # add entry to search cache
                    search_cache[name] = smi
                    print("ok")
                else:
                    print("NO SMILES FROM CID")
            else:
                print("NO CID")

    # add the SMILES structures to the master table by g_id for any compounds that were matched
    qry = "UPDATE master SET smi=? WHERE g_id=?"
    for g_id in gid_to_smi:
        qdata = (gid_to_smi[g_id], g_id)
        cursor.execute(qry, qdata)


if __name__ == '__main__':

    from sqlite3 import connect
    from requests import Session

    # connect to database
    con = connect("C3S.db")
    cur = con.cursor()

    # load name->SMILES search cache
    name_to_smi = load_search_cache()

    # requests Session
    sess = Session()

    print("fetching SMILES structures...")
    try:
        # fetch n_smi for each dataset (that does not already have it)
        fetch_smi_byname(cur, sess, name_to_smi)
        # save changes to the database
        con.commit()
    except KeyboardInterrupt:
        # save the search cache but do not write changes to the database
        print("\nINTERRUPTED!\nsaving search cache and exiting ...")
    except Exception as e:
        print()
        print(e)
        print("ERROR!\nsaving search cache and exiting ...")
    
    # close database connection
    con.close()

    # save the search cache to file
    save_search_cache(name_to_smi)

