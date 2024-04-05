"""
    c3sdb/test/__main__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Run all defined unit tests in the test sub-package
"""


import unittest

# import TestCase subclasses and unittest.main will run them all
from c3sdb.test.build_utils.db_init import (
    Test_IncludePath
)
from c3sdb.test.build_utils.smiles import (
    Test_SmilesSearchCache
)
from c3sdb.test.build_utils.src_data import (
    Test_SrcDataPath
)



if __name__ == "__main__":
    unittest.main(verbosity=2)

