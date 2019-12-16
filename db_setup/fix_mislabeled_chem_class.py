#!/usr/bin/python3
"""
    fix_mislabeled_chem_class.py
    Dylan H. Ross
    2019/08/31
    
        TODO
"""




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

    print("fixing mislabeled compounds ...")
    #do_the_thing(cur)
    # save changes to the database
    #con.commit()
    
    # close database connection
    con.close()

