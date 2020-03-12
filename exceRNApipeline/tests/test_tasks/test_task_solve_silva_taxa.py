import exceRNApipeline.tasks.task_solve_silva_taxa as task
from unittest import TestCase
import os
import sys
import shutil


_FILE_PATH = os.path.dirname(__file__)
_TEST_TEMP = os.path.join(_FILE_PATH, "_test_temp")
_DATA_DIR = os.path.join(_FILE_PATH, "data")

_BASE_ARGS = [
    '',
    '-f', os.path.join(_DATA_DIR, "silva.fasta.gz"),
    '-t', os.path.join(_DATA_DIR, "taxmap.txt"),
    '-s', 'lsu',
    '-o', os.path.join(_TEST_TEMP, 'silva.fasta.gz')
]

class TestTaskSolveSilvaTaxa(TestCase):

    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        if os.path.exists(_TEST_TEMP):
            shutil.rmtree(_TEST_TEMP)
        os.mkdir(_TEST_TEMP)
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
        shutil.rmtree(_TEST_TEMP)
    
    def test_solva_silva_taxa(self):
        sys.argv = _BASE_ARGS
        task.main()
