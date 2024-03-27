"""
    c3sdb/build_utils/mqns.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module for computing and adding MQNs to the database
"""


import sqlite3
from typing import Optional, List

from rdkit import Chem
from rdkit.Chem import Descriptors


def _get_mqns(smi: str
              ) -> Optional[List[float]] :
    """
    Computes the complete set of 42 MQNs as described in:
    Nguyen et al. ChemMedChem 4:1803-5 (2009)
    
    Parameters
    ----------
    smi : ``str``
        SMILES string of structure
    
    Returns
    -------
    mqns : ``list(float)`` or ``None`` 
        array of MQN features, None if anything goes wrong
    """
    try:
        features = Descriptors.rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smi))
        return features
    except Exception as e:
        return None


def add_mqns_to_db(cur: sqlite3.Cursor
                   ) -> None :
    """
    
    Parameters
    ----------
    cur : ``sqlite3.Cursor``
        C3S.db database cursor
    """
    # generate the rdk features and MQNs
    qry = "SELECT g_id, smi FROM master WHERE smi IS NOT NULL"
    gid_to_mqn = {}
    for g_id, n_smi in cur.execute(qry).fetchall():
        mqn = _get_mqns(n_smi)
        if mqn:
            gid_to_mqn[g_id] = mqn
    # update the database with the generated MQNs
    qry = f"INSERT INTO mqns VALUES ({','.join('?' * 43)})" 
    for g_id in gid_to_mqn:
        qdata = (g_id, *gid_to_mqn[g_id])
        cur.execute(qry, qdata)

