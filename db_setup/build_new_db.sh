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
