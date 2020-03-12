import exceRNApipeline.tasks.task_anno_to_fasta as task
from unittest import TestCase
import os
import shutil
import sys
from Bio import SeqIO


_TEMP_DIR = os.path.join(os.path.dirname(__file__), "_test_temp")
_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
_ANNO_FILE = os.path.join(_DATA_DIR, "anno_rm_hg38.txt")
_GENOME_FILE = os.path.join(_DATA_DIR, "hg38_chr1_trunc.fa")
_OUTPUT_FILE = os.path.join(_TEMP_DIR, "output.fa")


_BASE_ARGS = [
    "",
    "-a", _ANNO_FILE,
    "-g", _GENOME_FILE,
    "-o", _OUTPUT_FILE,
    "-c", "5",
    "-s", "6",
    "-e", "7",
    "-k", "3",
    "-n", "10", "11"
]

class TestAnnoToFasta(TestCase):

    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        if os.path.exists(_TEMP_DIR):
            shutil.rmtree(_TEMP_DIR)
        os.mkdir(_TEMP_DIR)
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
        shutil.rmtree(_TEMP_DIR)
    
    def test_anno_to_fasta(self):
        '''Test for anno_to_fasta.py'''
        sys.argv = _BASE_ARGS
        task.main()
        seqs = list(SeqIO.parse(_OUTPUT_FILE, 'fasta'))
        self.assertEqual(len(seqs), 10)
