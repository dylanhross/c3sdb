"""
    c3sdb/build_utils/__init__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    sub-package with utilities for building the database
"""

# build script
"""
#!/bin/bash
# create a new sqlite3 database using a SQL schema file then run the python initialization scripts
# to fill it with data from files in the cleaned_data directory and pull any missing information
# from other sources


# get rid of the database file (if it exists)
rm -f C3S.db

# initialize the database using the SQL schema file
cat C3SDB_schema.sql | sqlite3 C3S.db

# fill the database with initial values from the reference files
./fill_db_from_src.py

# get missing neutral SMILES structures
./fetch_smiles.py

# add the MQN table to the database then add MQNs for all available SMILES structures
cat mqn_schema.sql | sqlite3 C3S.db
./generate_mqns.py

# label each compound with a (rough) chemical class, then fix specific mislabeled compounds
./label_chem_class.py
#./fix_mislabeled_chem_class.py

# add the predicted CCS table to the database
cat pred_CCS_schema.sql | sqlite3 C3S.db
"""

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
    "leap0219",
    "blaz0818", 
    #"vasi0120",
    "tsug0220",
    "lian0118",
    "teja0918",
    "pola0620",
    "dodd0220",
    "celm1120",
    "belo0321"
]