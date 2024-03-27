"""
    c3sdb/build_utils/src_data.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Module for adding source datasets to the database
"""


import hashlib
import re
import json
import sqlite3


def _gen_id(name: str, 
            adduct: str, 
            ccs: float, 
            ccs_type: str, 
            src_tag: str
            ) -> str :
    """ 
    computes a unique string identifier for an entry by hashing on name+adduct+ccs+ccs_type+src_tag
    """
    s = f"{name}{adduct}{ccs}{ccs_type}{src_tag}"
    h = hashlib.sha1(s.encode()).hexdigest()[-10:].upper()
    return 'CCSBASE_' + h


def add_dataset(cursor: sqlite3.Cursor, 
                src_dset_file: str
                ) -> None :
    """
    Adds values from a source dataset to the database, specified by a source tag
    
    Parameters
    ----------
    cursor : ``sqlite3.cursor``
        cursor for running queries against the drugs.db database
    src_dset : ``str``
        source dataset (JSON file)
    """
    with open(src_dset_file, "r") as j:
        jdata = json.load(j)
    # regex for identifying multiple charges
    multi_z = re.compile('.*\]([0-9])[+-]')
    # query string
    # g_id, name, adduct, mz, ccs, smi, src_tag
    qry = "INSERT INTO master VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
    for cmpd in jdata["data"]:
        # strip whitespace off of the name
        name = cmpd["name"].strip()
        # fix messed up adducts on the fly
        adduct = cmpd["adduct"]
        adduct = {
            "[M+]+": "[M]+", 
            "M+NH4]+": "[M+NH4]+",
            "[M+H]+*": "[M+H]+",
            "[M+Na]+*": "[M+Na]+",
            "[M+H20-H]-": "[M+H2O-H]-",
        }.get(adduct, adduct)
        # check for multiple charges
        mz = float(cmpd["mz"])
        is_multi = multi_z.match(adduct)
        z = 1
        if is_multi:
            z = int(is_multi.group(1))
        # calculate the mass
        mass = mz * z
        # use smi if included
        smi = cmpd.get("smi")
        # make sure CCS is a float
        ccs = float(cmpd["ccs"])
        # metadata
        ccs_type, ccs_method = jdata["metadata"]["type"], jdata["metadata"]["method"]
        src_tag = jdata["metadata"]["src_tag"]
        # unique identifier
        qdata = (
            _gen_id(name, adduct, ccs, ccs_type, src_tag), 
            name, adduct, mass, z, mz, ccs, smi, None, src_tag, ccs_type, ccs_method
        )
        cursor.execute(qry, qdata)

