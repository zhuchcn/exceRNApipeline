import exceRNApipeline.tasks.task_extract_fastx as task
from unittest import TestCase
import os
import shutil
import sys


_FILE_PATH = os.path.dirname(__file__)
_DATA_DIR = os.path.join(_FILE_PATH, 'data')
_TEMP_DIR = os.path.join(_FILE_PATH, '_test_temp')

_BASE_ARGS = [
    '',
    '-i', os.path.join(_DATA_DIR, 'test.fa.gz'),
    '-o', os.path.join(_TEMP_DIR, 'test_output_fa.gz'),
    '-l', os.path.join(_DATA_DIR, 'fasta_namelist.txt'),
    '-f', 'fasta'
]

class TestTaskExtractFasta(TestCase):

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
    
    def test_extract_fasta(self):
        sys.argv = _BASE_ARGS
        task.main()