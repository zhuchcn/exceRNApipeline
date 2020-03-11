import src.tasks.task_preprocess as task
import src.tasks.slurm_job as slurm_job
from unittest import TestCase
from unittest.mock import patch
from io import StringIO
import sys
import re
import os


# overwrite the shell function because we don't want to really execute the 
# commands.
task.shell = lambda *args, **kwargs : None
slurm_job.shell = lambda *args, **kwargs : "".encode('utf-8')
task.SlurmJob = slurm_job.SlurmJob

_BASE_ARGS = [
    "",
    "-i", "input/some.fastq.gz",
    "-o", "output/preprocessed/some_SE.fastq.gz",
    "-n", "some",
    "-a", "TGGAATTCTCGGGTGCCAAGGAACTC",
    "-l", "output/preprocessed/some_SE.log",
    "-p", "output/preprocessed/some"
]

class TestTaskPreProcess(TestCase):
    def setUp(self):
        self.sys_argv = sys.argv
        super().setUp()
    
    def tearDown(self):
        sys.argv = self.sys_argv
        super().tearDown()
    
    @staticmethod
    def run_task():
        with patch('sys.stdout', new=StringIO()) as stdout:
            task.main()
            return stdout.getvalue()

    def test_preprocess_base(self):
        """Test for all flags except -s"""
        sys.argv = _BASE_ARGS 
        out = self.run_task()
        occurs = re.findall("hts_Stats", out)
        self.assertEqual(len(occurs), 2)
        self.assertTrue(bool(re.search("hts_AdapterTrimmer", out)))
        self.assertTrue(bool(re.search("hts_QWindowTrim", out)))
        self.assertTrue(bool(re.search("hts_NTrimmer", out)))
    
    def test_preprocess_scratch(self):
        """Test for -s"""
        sys.argv = _BASE_ARGS + ['-s', '/scratch']
        out = self.run_task()
        self.assertTrue(bool(re.search("mkdir /scratch/", out)))
        self.assertTrue(bool(re.search("mv /scratch/.+fastq\.gz", out)))
        self.assertTrue(bool(re.search("mv /scratch/.+_SE\.log", out)))