import exceRNApipeline.tasks.task_silva_extract_unmapped as task
import exceRNApipeline.tasks.slurm_job as slurm_job
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
task.extract_fastx = lambda *args, **kargs : None

_SAMPLE = "sample1"
_ALIGNED_TXT = f"output/06-SILVA/{_SAMPLE}/Aligned_combine.txt.gz"
_UNMAPPED_FQ = f"output/05-RE/{_SAMPLE}/Unmapped.fastq.gz"
_OUTPUT_FQ = f"output/06-SILVA/{_SAMPLE}/Unmapped.fastq.gz"
_NAMELIST_TXT = f"output/06-SILVA/{_SAMPLE}/aligned_reads.txt"

_BASE_ARGS = [
    '',
    '-a', _ALIGNED_TXT,
    '-u', _UNMAPPED_FQ,
    '-o', _OUTPUT_FQ,
    '-l', _NAMELIST_TXT
]


class TestTaskSilvaExtractUnmapped(TestCase):

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
    
    def test_silva_extract_unmapped(self):
        sys.argv = _BASE_ARGS
        out = self.run_task()
        print(out)
