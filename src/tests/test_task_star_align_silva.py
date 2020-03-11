import src.tasks.task_star_align_silva as task
import src.tasks.slurm_job as slurm_job
from unittest import TestCase
from unittest.mock import patch
import unittest
from io import StringIO
import sys
import os


# overwrite the shell function because we don't want to really execute the 
# commands.
task.shell = lambda *args, **kwargs : None
slurm_job.shell = lambda *args, **kwargs : "".encode('utf-8')
task.SlurmJob = slurm_job.SlurmJob

_SAMPLE_NAME = "sample1"
_SILVA_IND = "01"
_INPUT_FQ = f"output/05-RE/{_SAMPLE_NAME}/Unmapped.fastq.gz"
_GENOME_INDEX = f"genomes/silva_tax/star_index_silva_{_SILVA_IND}"
_OUTPUT_TXT = f"output/06-SILVA/{_SAMPLE_NAME}/Aligned_{_SILVA_IND}.txt.gz"
_OUTPUT_PREFIX = f"06-SILVA/{_SAMPLE_NAME}/Aligned_{_SILVA_IND}/"
_NTHREADS = "16"
_BASE_ARGS = [
    "",
    "-i", _INPUT_FQ,
    "-n", _SAMPLE_NAME,
    "-g", _GENOME_INDEX,
    "-o", _OUTPUT_TXT,
    "-p", _OUTPUT_PREFIX,
    "-t", _NTHREADS
]

class TestTaskStarAlignSilva(TestCase):
    
    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
    
    @staticmethod
    def run_task():
        with patch('sys.stdout', new=StringIO()) as stdout:
            task.main()
            return stdout.getvalue()
        
    def test_star_align_silva_base(self):
        """test for star align for base case"""
        sys.argv = _BASE_ARGS
        out = self.run_task()
        # print(out)
        # pattern = f"STAR\s+\W+\s+--runMode alignReads\s+\W+\s+" +\
        #           f"--readFilesIn {_INPUT_FQ}\s+\W+\s+" +\
        #           f"--genomeDir {_GENOME_INDEX}\s+\W+\s+" + \
        #           f"--outFileNamePrefix {_OUTPUT_PREFIX}\s+\W+\s+" +\
        #           f"--runThreadN {_NTHREADS}.+"
        #self.assertRegex(out, pattern)
    
    def test_star_align_silva_scratch(self):
        """test for star align with scratch"""
        sys.argv = _BASE_ARGS + ["-s", "/scratch"]
        out = self.run_task()
        # print(out)
        # pattern = f"STAR\s+\W+\s+--runMode alignReads\s+\W+\s+" +\
        #           f"--readFilesIn {_INPUT_FQ}\s+\W+\s+" +\
        #           f"--genomeDir {_GENOME_INDEX}\s+\W+\s+" + \
        #           f"--outFileNamePrefix /scratch/_/{_SAMPLE_NAME}\s+\W+\s+" +\
        #           f"--runThreadN {_NTHREADS}.+"
        # self.assertRegex(out, pattern)
        