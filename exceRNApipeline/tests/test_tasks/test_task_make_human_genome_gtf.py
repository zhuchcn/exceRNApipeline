import exceRNApipeline.tasks.task_make_human_genome_gtf as task
import unittest
from unittest.mock import patch
from io import StringIO
import sys
import os
import shutil


_FILE_PATH = os.path.dirname(__file__)
_OUTPUT_DIR = os.path.join(_FILE_PATH, '_test_temp')
_GENCODE_GTF = os.path.join(_OUTPUT_DIR, "gencode.gtf")
_TRNA_GTF = os.path.join(_OUTPUT_DIR, "tRNA.gtf")
_PIRNA_GTF = os.path.join(_OUTPUT_DIR, "piRNA.gtf")

_BASE_ARGS = [
    '',
    '-c', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz',
    '-r', 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.tRNAs.gtf.gz',
    '-i', 'https://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.v1_7_6.hg38.gtf.gz',
    '-g', _GENCODE_GTF,
    '-t', _TRNA_GTF,
    '-p', _PIRNA_GTF,
    '-o', _OUTPUT_DIR
]

@unittest.skip
class TestTaskMakeHumanGenomeGTF(unittest.TestCase):

    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        os.mkdir(_OUTPUT_DIR)
        if os.path.exists(_OUTPUT_DIR):
            shutil.rmtree(_OUTPUT_DIR)
    
    @staticmethod
    def run_task():
        with patch("sys.stdout", new=StringIO()) as stdout:
            task.main()
            return stdout.getvalue()
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
        shutil.rmtree(_OUTPUT_DIR)
    
    def test_make_human_genome_gtf(self):
        sys.argv = _BASE_ARGS
        self.run_task()