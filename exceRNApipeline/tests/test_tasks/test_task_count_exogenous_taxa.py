import exceRNApipeline.tasks.task_count_exogenous_taxa as task
from unittest import TestCase
import os
import shutil
import sys


_FILE_PATH = os.path.dirname(__file__)
_TEST_PATH = os.path.dirname(_FILE_PATH)
_DATA_DIR = os.path.join(os.path.dirname(_TEST_PATH), 'includes/taxonomy/tests')
_TEMP_DIR = os.path.join(_FILE_PATH, '_test_temp')

_BASE_ARGS = [
    "",
    "-i", os.path.join(_DATA_DIR, "sample_readCounts.txt"),
    "-o", os.path.join(_DATA_DIR, "sample_taxaCount_"),
    "-m", os.path.join(_DATA_DIR, "names_test2.dmp"),
    "-d", os.path.join(_DATA_DIR, "nodes_test2.dmp")
]

class TestTaskCountExogenousTaxa(TestCase):

    def setUp(self):
        super().setUp()
        if os.path.exists(_TEMP_DIR):
            shutil.rmtree(_TEMP_DIR)
        os.mkdir(_TEMP_DIR)
        self.sys_argv = sys.argv
    
    def tearDown(self):
        super().tearDown()
        shutil.rmtree(_TEMP_DIR)
        sys.argv = self.sys_argv
    
    def test_count_exogenous_taxa(self):
        sys.argv = _BASE_ARGS
        task.main()