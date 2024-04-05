"""
    c3sdb/test/build_utils/smiles.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Unit tests for the c3sdb.build_utils.smiles module
"""


import unittest
import os

from c3sdb.build_utils.smiles import _SMILES_SEARCH_CACHE


class Test_SmilesSearchCache(unittest.TestCase):
    """ tests for the _SMILES_SEARCH_CACHE constant """

    def test_smiles_search_cache_exists(self):
        """ ensure the smiles search cache file is found under _SMILES_SEARCH_CACHE """
        self.assertTrue(os.path.isfile(_SMILES_SEARCH_CACHE))
