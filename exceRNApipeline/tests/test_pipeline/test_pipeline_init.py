import exceRNApipeline.pipeline.__main__ as pipeline
from unittest import TestCase
import sys
import os
import shutil


_FILE_DIR = os.path.dirname(__file__)
_TEMP_DIR = os.path.join(_FILE_DIR, "_test_tmp")

class TestCasePipelineInit(TestCase):
    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        self.cwd = os.getcwd()
        if os.path.exists(_TEMP_DIR):
            shutil.rmtree(_TEMP_DIR)
        os.mkdir(_TEMP_DIR)
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
        shutil.rmtree(_TEMP_DIR)
        os.chdir(self.cwd)
    
    def test_pipeline_init(self):
        os.chdir(_TEMP_DIR)
        sys.argv = ["", "init"]
        pipeline.main()
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "input")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "genomes")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "output")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "slurmout")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "slurm_job_status")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "pipeline_config.yml")))
        self.assertTrue(os.path.exists(os.path.join(_TEMP_DIR, "slurm_config.yml")))