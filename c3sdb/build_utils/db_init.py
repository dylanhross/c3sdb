"""
    c3sdb/build_utils/db_init.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Module for inititializing the database file
"""


import os
import sqlite3


# store path to _include directory in this package
_INCLUDE_PATH: str = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "_include/")


def create_db(f: str
              ) -> None :
    """
    initialize the C3S.db from SQLite3 scripts (included in this package in the c3sdb/_include/ directory)
    
    .. note:: 

        overwrites the database file if it already exists

    Parameters
    ----------
    f : ``str``
        filename/path of the database
    """
    # see if the file exists
    if os.path.exists(f):
        os.remove(f)
    # initial connection creates the DB
    con = sqlite3.connect(f)  
    cur = con.cursor()
    # execute SQL scripts to set up the database
    sql_scripts = [
        os.path.join(_INCLUDE_PATH, "C3SDB_schema.sqlite3"),
        os.path.join(_INCLUDE_PATH, "mqn_schema.sqlite3"),
        os.path.join(_INCLUDE_PATH, "pred_CCS_schema.sqlite3")
    ]
    for sql_script in sql_scripts:
        with open(sql_script, "r") as sql_f:
            cur.executescript(sql_f.read())
    # save and close the database
    con.commit()
    con.close()

