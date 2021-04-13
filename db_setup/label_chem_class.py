#!/usr/bin/python3
"""
    label_chem_class.py
    Dylan H. Ross
    2019/07/17
    
        TODO
"""


# change how we import other tools/build scripts depending upon whether this is being called directly
# or being used in the context of the db_setup module (i.e. to generate documentation)
if __name__ == '__main__':
    from lipid_parser import parse_lipid
    from peptide_parser import parse_peptide
    from carbohydrate_parser import parse_carbohydrate

else:
    from db_setup.lipid_parser import parse_lipid
    from db_setup.peptide_parser import parse_peptide
    from db_setup.carbohydrate_parser import parse_carbohydrate


def label_class_byname(cursor):
    """
label_class_byname
    description:
        adds rough chemical class labels to the master table of C3S.db based on the name. Currently, if the name can
        be parsed as a lipid it will be assigned the 'lipid' class, likewise for the 'peptide' and 'carbohydrate' 
        classes, and failing all of those it is simply assigned 'small molecule'.
    parameters:
        cursor (sqlite3.cursor) -- cursor for running queries against the C3S.db database
"""
    # map global identifiers to class labels
    gid_to_lab = {}

    qry = "SELECT g_id, name FROM master WHERE chem_class_label IS NULL"
    for g_id, name in cursor.execute(qry).fetchall():
        
        print('(g_id: {}) {:40s} looks like a'.format(g_id, name), end=' ')

        if parse_lipid(name):
            label = 'lipid'
        elif parse_peptide(name):
            label = 'peptide'
        elif parse_carbohydrate(name):
            label = 'carbohydrate'
        else:
            label = 'small molecule'
  
        print(label)

        gid_to_lab[g_id] = label


    # add the chem class label to the master table by g_id for any compounds that were matched
    qry = "UPDATE master SET chem_class_label=? WHERE g_id=?"
    for g_id in gid_to_lab:
        qdata = (gid_to_lab[g_id], g_id)
        cursor.execute(qry, qdata)



if __name__ == '__main__':

    from sqlite3 import connect

    # connect to database
    con = connect("C3S.db")
    cur = con.cursor()

    print("adding (rough) chemical class labels ...")
    label_class_byname(cur)
    # save changes to the database
    con.commit()
    
    # close database connection
    con.close()



