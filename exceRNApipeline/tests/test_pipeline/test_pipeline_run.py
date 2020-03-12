import exceRNApipeline.pipeline.__main__ as pipeline
import exceRNApipeline.pipeline.run as run
from unittest import TestCase
from unittest.mock import patch
from io import StringIO
import sys
import os
import shutil


_FILE_DIR = os.path.dirname(__file__)
_TEMP_DIR = os.path.join(_FILE_DIR, "_test_tmp")

run.sp.run = lambda *args, **kwargs: None
pipeline.parse_args_run = run.parse_args


class TestCasePipelineRun(TestCase):
    """Test pipeline run"""

    def setUp(self):
        super().setUp()
        self.sys_argv = sys.argv
        self.cwd = os.getcwd()
        if os.path.exists(_TEMP_DIR):
            shutil.rmtree(_TEMP_DIR)
        os.mkdir(_TEMP_DIR)
        self.pipeline_init()
    
    def tearDown(self):
        super().tearDown()
        sys.argv = self.sys_argv
        shutil.rmtree(_TEMP_DIR)
        os.chdir(self.cwd)
    
    def pipeline_init(self):
        os.chdir(_TEMP_DIR)
        sys.argv = ["", "init"]
        pipeline.main()
        os.chdir(self.cwd)
    
    @staticmethod
    def run_pipeline():
        with patch('sys.stdout', new=StringIO()) as stdout:
            pipeline.main()
            return stdout.getvalue()
    
    def test_pipeline_run(self):
        """Test for pipeline run -np"""
        os.chdir(_TEMP_DIR)
        sys.argv = [
            "", "run",
            "--configfile", "pipeline_config.yml",
            "-np"
        ]
        out = self.run_pipeline()
        self.assertRegex(out, "^snakemake --snakefile " +\
            ".+/pipeline/smk/Snakefile --configfile pipeline_config.yml.*")
    
    def test_pipeline_run_slurm(self):
        """Test for pipeline run slurm"""
        os.chdir(_TEMP_DIR)
        sys.argv = [
            "", "run",
            "--configfile", "pipeline_config.yml",
            "--slurm-config", "slurm_config.yml"
        ]
        out = self.run_pipeline()
        self.assertTrue(out, "srun.+snakemake --snakefile.+--cluster-config"+\
            " slurm_config.yml --cluster \"sbatch.+\"")
    
    def test_pipeline_run_singularity(self):
        """Test for pipeline run with singularity"""
        os.chdir(_TEMP_DIR)
        sys.argv = [
            "", "run",
            "--configfile", "pipeline_config.yml",
            "--slurm-config", "slurm_config.yml",
            "--use-singularity"
        ]
        out = self.run_pipeline()
        print(out)
        self.assertTrue(out, "srun.+snakemake --snakefile.+ "+\
            "--use-singularity" +\
            " --cluster-config slurm_config.yml --cluster \"sbatch.+\"")
    