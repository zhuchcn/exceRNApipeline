import exceRNApipeline.tasks.task_star_index as task
import exceRNApipeline.tasks.slurm_job as slurm_job
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
    "-i", "genome/some.fasta",
    "-o", "genome/some_star_index",
    "-t", "16"
]

class TestTaskStarIndex(TestCase):

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
    
    def test_star_index_base(self):
        """Test for -i, -o, and -t"""
        sys.argv = _BASE_ARGS
        out = self.run_task()

        self.assertTrue(bool(re.search("cd genome", out)))
        self.assertFalse(bool(re.search("STAR None", out)))
        self.assertTrue(bool(re.search("--runMode genomeGenerate", out)))
        self.assertTrue(bool(
            re.search(f"--genomeFastaFiles {os.getcwd()}/genome/some.fasta",
            out
        )))
        self.assertTrue(bool(re.search("--genomeDir some_star_index", out)))
        self.assertTrue(bool(
            re.search(f"--limitGenomeGenerateRAM {16*2*1024**3}", out
        )))

    def test_star_index_mem_gb(self):
        """Test for -m"""
        sys.argv = _BASE_ARGS + ['-m', '60']
        out = self.run_task()
        
        self.assertTrue(bool(
            re.search(f"--limitGenomeGenerateRAM {60*1024**3}", out)
        ))
    
    def test_star_index_extra_args(self):
        """Test for -a"""
        extra_args = "--genomeSAindexNbases 8 --genomeChrBinNbits 16"
        sys.argv = _BASE_ARGS + ["-a", extra_args]
        out = self.run_task()

        self.assertTrue(bool(re.search(extra_args, out)))
    
    def test_star_index_scratch(self):
        """Test for -s"""
        sys.argv = _BASE_ARGS + ['-s', '/scratch']
        out = self.run_task()
        self.assertTrue(bool(re.search("mkdir /scratch/", out)))
        self.assertTrue(bool(re.search(
            "mv /scratch/_/some_star_index genome/some_star_index", out
        )))
