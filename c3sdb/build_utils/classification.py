"""
    c3sdb/build_utils/classification.py

    Dylan Ross (dylan.ross@pnnl.gov)

    functions for performing chemcial classifications
"""


import sqlite3

from c3sdb.build_utils._parsing import parse_carbohydrate, parse_lipid, parse_peptide


def label_class_byname(cursor: sqlite3.Cursor
                       ) -> None :
    """
    adds rough chemical class labels to the master table of C3S.db based on the name
    
    Currently, if the name can be parsed as a lipid it will be assigned the 'lipid' 
    class, likewise for the 'peptide' and 'carbohydrate' classes, and failing all 
    of those it is simply assigned 'small molecule'.
    
    Parameters
    ----------
    cursor : ``sqlite3.Cursor``
        cursor for C3S.db database
    """
    # map global identifiers to class labels
    gid_to_lab = {}
    qry = "SELECT g_id, name FROM master WHERE chem_class_label IS NULL"
    for g_id, name in cursor.execute(qry).fetchall():
        #print('(g_id: {}) {:40s} looks like a'.format(g_id, name), end=' ')
        if parse_lipid(name):
            label = 'lipid'
        elif parse_peptide(name):
            label = 'peptide'
        elif parse_carbohydrate(name):
            label = 'carbohydrate'
        else:
            label = 'small molecule'
        #print(label)
        gid_to_lab[g_id] = label
    # add the chem class label to the master table by g_id for any compounds that were matched
    qry = "UPDATE master SET chem_class_label=? WHERE g_id=?"
    for g_id in gid_to_lab:
        qdata = (gid_to_lab[g_id], g_id)
        cursor.execute(qry, qdata)
