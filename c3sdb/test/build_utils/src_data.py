"""
    c3sdb/test/build_utils/src_data.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Unit tests for the c3sdb.build_utils.src_data module
"""


import unittest
import os

from c3sdb.build_utils.src_data import _SRC_DATA_PATH, ALL_SRC_TAGS


class Test_SrcDataPath(unittest.TestCase):
    """ tests for the _SRC_DATA_PATH constant """

    def test_srd_data_path_dir_exists(self):
        """ ensure the directory pointed to by _SRC_DATA_PATH exists """
        self.assertTrue(os.path.isdir(_SRC_DATA_PATH))

    def test_src_data_path_has_source_datasets(self):
        """ ensure all of the expected source datasets are present under _SRC_DATA_PATH """
        for src_tag in ALL_SRC_TAGS:
            src_data_path = os.path.join(_SRC_DATA_PATH, f"{src_tag}.json")
            self.assertTrue(os.path.isfile(src_data_path))
