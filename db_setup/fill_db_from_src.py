#!/usr/bin/python3
"""
    fill_db_from_src.py
    Dylan H. Ross
    2018/10/04

        Initializes the C3S.db Sqlite3 database, filling it with values from the the cleaned datasets.
        After this script runs, the database will contain only RELEVANT information contained in the 
        original source datasets (name, adduct, mz, ccs, and smi if provided).
"""


from json import load as jload
import re


def add_src_dataset(cursor, src_tag, metadata, gid_start=0):
    """
add_src_dataset
    description:
        Adds values from a source dataset to the database, specified by a source tag
    parameters:
        cursor (sqlite3.cursor) -- cursor for running queries against the drugs.db database
        src_tag (str) -- source tag 
        metadata (dict(...)) -- CCS metadata: CCS type and method
        [gid_start (int)] -- starting number for g_id integer identifier [optional, default=0]
    returns:
        (int) -- the next available g_id value
"""
    # regex for identifying multiple charges
    multi_z = re.compile('.*\]([0-9])[+-]')

    with open("cleaned_data/{}.json".format(src_tag), "r") as j:
        jdata = jload(j)
    # query string
    # g_id, name, adduct, mz, ccs, smi, src_tag
    qry = "INSERT INTO master VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
    # s_id starts at 0 and goes up from there
    g_id = gid_start
    for cmpd in jdata:

        # fix messed up adducts on the fly...
        adduct = cmpd["adduct"]
        fixed_adducts = {
            "[M+]+": "[M]+", 
            "M+NH4]+": "[M+NH4]+",
            "[M+H]+*": "[M+H]+",
            "[M+Na]+*": "[M+Na]+",
            "[M+H20-H]-": "[M+H2O-H]-"
        }
        if adduct in fixed_adducts:
            adduct = fixed_adducts[adduct]

        # check for multiple charges
        mz = cmpd["mz"]
        adduct = cmpd["adduct"]
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
        
        # CCS metadata
        ccs_type, ccs_method = metadata["type"], metadata["method"]

        qdata = (
            g_id, cmpd["name"].strip(), adduct, mass, z, mz, cmpd["ccs"], smi, None, src_tag, ccs_type, ccs_method
        )
        cursor.execute(qry, qdata)
        g_id += 1

    return g_id


if __name__ == '__main__':

    from sqlite3 import connect


    # connect to database
    con = connect("C3S.db")
    cur = con.cursor()

    # source datasets
    dsets = [
        "zhou1016",
        "zhou0817",
        "zhen0917",
        "pagl0314",
        "righ0218",
        "nich1118",
        "may_0114",
        "moll0218",
        "hine1217",
        "hine0217",
        "hine0817",
        "groe0815",
        "bijl0517",
        "stow0817",
        "hine0119",
        "leap0219"
    ]

    # CCS metadata by source
    metadata = {
        "zhou1016": {"type": "DT", "method": "single field, calibrated with Agilent tune mix (Agilent)"},
        "zhou0817": {"type": "DT", "method": "single field, calibrated with Agilent tune mix (Agilent)"},
        "zhen0917": {"type": "DT", "method": "stepped-field"},
        "pagl0314": {"type": "TW", "method": "calibrated with polyalanine"},
        "righ0218": {"type": "TW", "method": "calibrated with polyalanine"},
        "nich1118": {"type": "DT", "method": "single field, calibrated with ESI Low Concentration Tuning Mix (Agilent)"},
        "may_0114": {"type": "DT", "method": "stepped-field"},
        "moll0218": {"type": "TW", "method": "Major Mix IMS/Tof Calibration Kit (Waters)"},
        "hine1217": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)"},
        "hine0217": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-)"},
        "hine0817": {"type": "TW", "method": "calibrated with polyalanine and drug standards"},
        "groe0815": {"type": "DT", "method": "stepped-field"},
        "bijl0517": {"type": "TW", "method": "?"},
        "stow0817": {"type": "DT", "method": "stepped-field"},
        "hine0119": {"type": "TW", "method": "calibrated with phosphatidylcholines (ESI+) and phosphatidylethanolamines (ESI-), doubly charged cardiolipins calibrated with poly-DL-alanine"},
        "leap0219": {"type": "DT", "method": "stepped-field"}
    }

    # add each src dataset
    print("adding cleaned datasets into C3S.db")
    gid_next = 0
    for dset in dsets:
        print("\tadding dataset: {} ...".format(dset), end=" ")
        gid_next = add_src_dataset(cur, dset, metadata[dset], gid_start=gid_next)
        print("ok")
    print()

    # save changes to the database
    con.commit()
    con.close()
