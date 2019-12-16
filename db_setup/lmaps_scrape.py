"""
    fetch_lmaps_smiles.py
    Dylan Ross
    2018/10/15

        Uses the Lipid MAPS REST API to search the LMSD for lipid SMILES structures.
"""


from time import sleep


def str_from_lipid_dict(lipid, ign_fa_comp=False):
    """
str_from_lipid_dict
    description:
        Converts a lipid definition in the form of a dictionary reflecting the structure of the lipids.db
        schema into a single string with standardized formatting such that it may be recognized by the
        Lipid MAPS structure database REST API.
    parameters:
        lipid (dict(...)) -- dictionary defining target lipid with a structure reflecting that of the 
                                lipids.db schema
        [ign_fa_comp (bool)] -- ignore individual fatty acid composition and just construct the name using
                                the total fatty acid composition [optional, default=False]
    returns:
        (str) -- string version of the lipid dict
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
    


def search_lmaps_smiles(sess, lipid, delay=0.1):
    """
search_lmaps_smiles
    description:
        uses the Lipid MAPS REST API to search the LMSD for lipid SMILES structures. Takes as input
        a lipid in the form of a dictionary with a structure reflecting that of lipids.db schema.
        If individual fatty acid composition is specified, then try searching using the 'abbrev_chains'
        input item first, then the 'abbrev' input item if that fails.
    parameters:
        sess (requests.Session) -- web request session
        lipid (dict(...)) -- dictionary defining target lipid with a structure reflecting that of the 
                                lipids.db schema
        delay (float) -- delay in seconds between consecutive web requests (so as not to hammer the
                            Lipid MAPS servers) [optional, default=0.1]
    returns:
        (str) -- SMILES structure or empty string in case of error
"""
    url = "https://lipidmaps.org/rest/compound/abbrev{}/{}/smiles/"
    lipid_name = str_from_lipid_dict(lipid)

    result = None
    if "_" in lipid_name:
        # includes individual fatty acid composition
        try:
            sleep(delay)
            result = sess.get(url.format("_chains", lipid_name)).json()
        except:
            msg = "fetch_lipid_smiles: search with 'abbrev_chains' for lipid {} failed".format(lipid_name)
            print(msg, end=" ")
        # try again with 
        if not result:
            try:
                sleep(delay)
                result = sess.get(url.format("", lipid_name)).json()
            except Exception as e:
                print(e)
                msg = "fetch_lipid_smiles: search with 'abbrev' for lipid {} failed".format(lipid_name)
                print(msg, end=" ")
        # regenerate the name but with total fatty acid composition in case the above did not work
        lipid_name = str_from_lipid_dict(lipid, ign_fa_comp=True)
    if not result:
        try:
            sleep(delay)
            result = sess.get(url.format("", lipid_name)).json()
        except:
            msg = "fetch_lipid_smiles: search with 'abbrev' for lipid {} failed".format(lipid_name)
            print(msg, end=" ")

    # if there was a result, just return the first one
    smi = ""
    if result:
        if 'Row1' in result:
            smi = result['Row1']['smiles']
        elif 'smiles' in result:
            smi = result['smiles']
    return smi

