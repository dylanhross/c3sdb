"""
    c3sdb/test/build_utils/db_init.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Unit tests for the c3sdb.build_utils.db_init module
"""


import unittest
import os

from c3sdb.build_utils.db_init import _INCLUDE_PATH


class Test_IncludePath(unittest.TestCase):
    """ tests for the _INCLUDE_PATH constant """

    def test_include_path_dir_exists(self):
        """ ensure the directory pointed to by _INCLUDE_PATH exists """
        self.assertTrue(os.path.isdir(_INCLUDE_PATH))

    def test_include_path_has_schemas(self):
        """ ensure that the expected SQLite3 schema files are found under _INCLUDE_PATH """
        expected_schema_files = [
            "C3SDB_schema.sqlite3", 
            "mqn_schema.sqlite3", 
            "pred_CCS_schema.sqlite3"
        ]
        for schema_file in expected_schema_files:
            sf_path = os.path.join(_INCLUDE_PATH, schema_file)
            self.assertTrue(os.path.isfile(sf_path))

