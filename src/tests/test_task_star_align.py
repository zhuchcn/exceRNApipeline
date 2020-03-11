import src.tasks.task_star_align as task
import src.tasks.slurm_job as slurm_job
from unittest import TestCase
from unittest.mock import patch
from io import StringIO
import sys
import os


# overwrite the shell function because we don't want to really execute the 
# commands.
task.shell = lambda *args, **kwargs : None
slurm_job.shell = lambda *args, **kwargs : "".encode('utf-8')
task.SlurmJob = slurm_job.SlurmJob

_INPUT_FQ = "output/preprocess/sample1_SE.fastq.gz"
_SAMPLE_NAME = "sample1"
_GENOME_INDEX = "genome/genome_star_index"
_OUTPUT_BAM = "output/align/sample1/sample1.bam"
_OUTPUT_UNMAPPED = "output/align/sample1/Unmapped.fastq.gz"
_OUTPUT_PREFIX = "output/align/sample1/"
_NTHREADS = "16"
_BASE_ARGS = [
    "",
    "-i", _INPUT_FQ,
    "-n", _SAMPLE_NAME,
    "-g", _GENOME_INDEX,
    "-b", _OUTPUT_BAM,
    "-u", _OUTPUT_UNMAPPED,
    "-p", _OUTPUT_PREFIX,
    "-t", _NTHREADS
]

class TestTaskStarAlign(TestCase):
    
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
        
    def test_star_align_base(self):
        """test for star align for base case"""
        sys.argv = _BASE_ARGS
        out = self.run_task()
        # print(out)
        pattern = f"STAR\s+\W+\s+--runMode alignReads\s+\W+\s+" +\
                  f"--readFilesIn {_INPUT_FQ}\s+\W+\s+" +\
                  f"--genomeDir {_GENOME_INDEX}\s+\W+\s+" + \
                  f"--outFileNamePrefix {_OUTPUT_PREFIX}\s+\W+\s+" +\
                  f"--runThreadN {_NTHREADS}.+"
        self.assertRegex(out, pattern)
    
    def test_star_align_scratch(self):
        """test for star align with scratch"""
        sys.argv = _BASE_ARGS + ["-s", "/scratch"]
        out = self.run_task()
        # print(out)
        pattern = f"STAR\s+\W+\s+--runMode alignReads\s+\W+\s+" +\
                  f"--readFilesIn {_INPUT_FQ}\s+\W+\s+" +\
                  f"--genomeDir {_GENOME_INDEX}\s+\W+\s+" + \
                  f"--outFileNamePrefix /scratch/_/{_SAMPLE_NAME}\s+\W+\s+" +\
                  f"--runThreadN {_NTHREADS}.+"
        self.assertRegex(out, pattern)
        