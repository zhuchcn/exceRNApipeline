import os
import shutil
from snakemake.shell import shell
from datetime import datetime


class SlurmJob():
    def __init__(self):
        print('Slurm init')
        self.hostname = shell('hostname', read=True).decode('utf-8').rstrip()
        self.slurm_jobid = shell("echo ${{SLURM_JOBID}}", read=True)\
            .decode('utf-8').rstrip()
        self.user = shell('echo ${{USER}}', read=True)\
            .decode('utf-8').rstrip()
    
    @property
    def scratch(self):
        return f'/scratch/{self.user}_{self.slurm_jobid}'
    
    @property
    def logfile(self):
        return f'slurm_job_status/{self.user}_{self.slurm_jobid}'
    
    def setUp(self):
        os.mkdir(self.scratch)
        with open(self.logfile, "w") as fh:
            fh.write(self.hostname)
    
    def tearDown(self):
        shutil.rmtree(self.scratch)
        os.remove(self.logfile)


def log(msg):
    print(
        '[ ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' ]' + msg,
        flush=True
    )