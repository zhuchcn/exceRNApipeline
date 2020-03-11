import src.tasks.task_split_fasta as task
import unittest
import os
import shutil
import sys


_FILE_DIR = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_FILE_DIR, 'data')
_INPUT_FA = os.path.join(_DATA_DIR, 'test.fa.gz')
_TEMP_DIR = os.path.join(_FILE_DIR, '_test_temp')
_OUTPUT_PREFIX = os.path.join(_TEMP_DIR, 'split-output')
_BASE_ARGS = [
    '',
    '-i', _INPUT_FA,
    '-o', _OUTPUT_PREFIX
]

class TestTaskSplitFasta(unittest.TestCase):

    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        if os.path.exists(_TEMP_DIR):
            shutil.rmtree(_TEMP_DIR)
        os.mkdir(_TEMP_DIR)
    
    def tearDown(self):
        sys.argv = self.sys_argv
        shutil.rmtree(_TEMP_DIR)
        super().tearDown()
    
    def test_split_fasta_r(self):
        sys.argv = _BASE_ARGS + ['-r', '5']
        task.main()
        self.assertEqual(len(os.listdir(_TEMP_DIR)), 2)
    
    def test_split_fasta_b(self):
        sys.argv = _BASE_ARGS + ['-b', '2']
        task.main()
        self.assertEqual(len(os.listdir(_TEMP_DIR)), 2)