"""    
    C3SData/db_interface.py
    Dylan H. Ross
    2019/04/09
    
    description:
        Functions for getting data from the C3S.db SQLite3 database by source dataset

"""


from sqlite3 import connect
from numpy import array


def all_dset(db_path):
    """
all_dset
    description:
        Returns a list of all src_tags in the C3S.db
    parameters:
        db_path (str) -- path to C3S.db database file
    returns:
        (list(str)) -- list of src_tags in C3S.db
"""
    con = connect(db_path)
    cur = con.cursor()
    qry = 'SELECT DISTINCT src_tag FROM master'
    src_tags = [_[0] for _ in cur.execute(qry).fetchall()]
    con.close()
    return src_tags


def fetch_single_dset(db_path, src_tag):
    """
fetch_single_dset
    description:
        fetch name, mz, adduct, ccs, src_tag, smi, and chem_class_label for a dataset based on the provided src_tag
    parameters:
        db_path (str) -- path to C3S.db database file
        src_tag (str) -- specify the source dataset
    returns:
        (numpy.array, ...) -- arrays of the desired attributes
"""
    con = connect(db_path)
    cur = con.cursor()
    qryA = 'SELECT g_id, name, mz, adduct, ccs, src_tag, smi, chem_class_label FROM master WHERE src_tag="{}"'.format(src_tag) + \
          ' AND smi IS NOT NULL AND g_id IN (SELECT g_id FROM mqns)'
    qryB = 'SELECT mqns.* FROM mqns INNER JOIN master ON mqns.g_id=master.g_id WHERE master.src_tag="{}"'.format(src_tag)
    names, mzs, adducts, ccss, srcs, smis, cls_labs = [], [], [], [], [], [], []
    for _, name, mz, adduct, ccs, src, smi, cls_lab in cur.execute(qryA).fetchall():
        names.append(name)
        mzs.append(mz)
        adducts.append(adduct)
        ccss.append(ccs)
        srcs.append(src)
        smis.append(smi)
        cls_labs.append(cls_lab)
    # fetch mqns separately from the rest of the data
    mqns = []
    for _, *mqn in cur.execute(qryB).fetchall():
        mqns.append(mqn)
    con.close()
    # self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_
    return array(names), array(mzs), array(adducts), array(ccss), array(srcs), array(smis), array(mqns), array(cls_labs)


def fetch_multi_dset(db_path, src_tags):
    """
fetch_multi_dset
    description:
        fetch name, mz, adduct, ccs, src_tag, smi and chem_class_label for a dataset based on the provided list of 
        src_tags. 
        * Internally, this function does not use multiple calls to fetch_single_dset, rather, it structures a 
        single SQL query to fetch all of the data at once for sake of efficiency. *
    parameters:
        db_path (str) -- path to C3S.db database file
        src_tags (list(str)) -- list of src_tags to specify source datasets
    returns:
        (numpy.array, ...) -- arrays of the desired attributes
"""
    con = connect(db_path)
    cur = con.cursor()
    qryA = 'SELECT g_id, name, mz, adduct, ccs, src_tag, smi, chem_class_label FROM master WHERE src_tag IN ('
    qryB = 'SELECT mqns.* FROM mqns INNER JOIN master ON mqns.g_id=master.g_id WHERE master.src_tag IN ('
    for tag in src_tags:
        qryA += '"{}",'.format(tag)
        qryB += '"{}",'.format(tag)
    qryA = qryA.rstrip(',') + ') AND smi IS NOT NULL AND g_id IN (SELECT g_id FROM mqns)'
    qryB = qryB.rstrip(',') + ')'
    names, mzs, adducts, ccss, srcs, smis, cls_labs = [], [], [], [], [], [], []
    for _, name, mz, adduct, ccs, src, smi, cls_lab in cur.execute(qryA).fetchall():
        names.append(name)
        mzs.append(mz)
        adducts.append(adduct)
        ccss.append(ccs)
        srcs.append(src)
        smis.append(smi)
        cls_labs.append(cls_lab)
    # fetch mqns separately from the rest of the data
    mqns = []
    for _, *mqn in cur.execute(qryB).fetchall():
        mqns.append(mqn)
    con.close()
    # self.cmpd_, self.mz_, self.adduct_, self.ccs_, self.src_, self.smi_, self.mqn_
    return array(names), array(mzs), array(adducts), array(ccss), array(srcs), array(smis), array(mqns), array(cls_labs)






