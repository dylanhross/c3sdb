"""
    c3sdb/build_utils/src_data.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Module for adding source datasets to the database
"""


import hashlib
import re


def _gen_id(name: str, 
            adduct: str, 
            ccs: float, 
            ccs_type: str, 
            src_tag: str
            ) -> str :
    """ 
    computes a unique string identifier for an entry by hashing on name+adduct+ccs+ccs_type+src_tag
    """
    s = '{}{}{}{}{}'.format(name, adduct, ccs, ccs_type, src_tag)
    h = hashlib.sha1(s.encode()).hexdigest()[-10:].upper()
    return 'CCSBASE_' + h


def add_dataset(cursor, src_tag, metadata):
    """
    Adds values from a source dataset to the database, specified by a source tag
    
    Parameters
    ----------
        cursor (sqlite3.cursor) -- cursor for running queries against the drugs.db database
        src_tag (str) -- source tag 
        metadata (dict(...)) -- CCS metadata: CCS type and method
        [gid_start (int)] -- starting number for g_id integer identifier [optional, default=0]
    
    Returns
    -------
    (int) -- the next available g_id value
"""
    # regex for identifying multiple charges
    multi_z = re.compile('.*\]([0-9])[+-]')

    with open("cleaned_data/{}.json".format(src_tag), "r") as j:
        jdata = jload(j)
    # query string
    # g_id, name, adduct, mz, ccs, smi, src_tag
    qry = "INSERT INTO master VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"

    # keep track of identifiers, skip the compound if there is a collision
    g_ids = []
    
    for cmpd in jdata:

        # strip whitespace off of the name
        name = cmpd["name"].strip()

        # fix messed up adducts on the fly...
        adduct = cmpd["adduct"]
        fixed_adducts = {
            "[M+]+": "[M]+", 
            "M+NH4]+": "[M+NH4]+",
            "[M+H]+*": "[M+H]+",
            "[M+Na]+*": "[M+Na]+",
            "[M+H20-H]-": "[M+H2O-H]-",
        }
        adduct = fixed_adducts[adduct] if adduct in fixed_adducts else adduct

        # check for multiple charges
        mz = float(cmpd["mz"])
        is_multi = multi_z.match(adduct)
        z = 1
        if is_multi:
            z = int(is_multi.group(1))

        # calculate the mass
        mass = mz * z

        # use smi if included
        smi = None
        if "smi" in cmpd:
            smi = cmpd["smi"]
        
        # make sure CCS is a float
        ccs = float(cmpd["ccs"])

        # CCS metadata
        ccs_type, ccs_method = metadata["type"], metadata["method"]

        # unique identifier
        g_id = gen_id(name, adduct, ccs, ccs_type, src_tag)

        if g_id not in g_ids:
            qdata = (
                g_id, name, adduct, mass, z, mz, ccs, smi, None, src_tag, ccs_type, ccs_method
            )
            cursor.execute(qry, qdata)
            g_ids.append(g_id)
        else:
            print('\t\tID: {} already present ({}, {}, {}, {}, {})'.format(g_id, name, adduct, ccs, ccs_type, src_tag))

