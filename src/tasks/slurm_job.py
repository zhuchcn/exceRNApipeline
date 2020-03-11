from snakemake.shell import shell
from datetime import datetime
from src.utils import logger


class SlurmJob():
    def __init__(self, scratch):
        print('Slurm init')
        self.hostname = shell('hostname', read=True).decode('utf-8').rstrip()
        self.slurm_jobid = shell("echo ${{SLURM_JOBID}}", read=True)\
            .decode('utf-8').rstrip()
        self.user = shell('echo ${{USER}}', read=True)\
            .decode('utf-8').rstrip()
        self.scratch_path = scratch
    
    @property
    def scratch(self):
        return f'{self.scratch_path}/{self.user}_{self.slurm_jobid}'
    
    @property
    def logfile(self):
        return f'slurm_job_status/{self.user}_{self.slurm_jobid}'
    
    def __enter__(self):
        self.setUp()
        return self
    
    def __exit__(self, type, value, traceback):
        self.tearDown()

    def setUp(self):
        cmd = f"""
        mkdir {self.scratch}
        echo {self.hostname} > {self.logfile}
        """
        logger(cmd)
        shell(cmd)

    def tearDown(self):
        cmd = f"""
        rm -rf {self.scratch}
        rm {self.logfile}
        """
        logger(cmd)
        shell(cmd)
