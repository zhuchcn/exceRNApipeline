import os
import shutil


def init(args):
    smk_dir = os.path.join(os.path.dirname(__file__), "smk")
    shutil.copy2(
        os.path.join(smk_dir, "pipeline_config.yml"),
        os.path.join(os.getcwd(), "pipeline_config.yml")
    )
    shutil.copy2(
        os.path.join(smk_dir, "slurm_config.yml"),
        os.path.join(os.getcwd(), "slurm_config.yml")
    )
    os.mkdir(os.path.join(os.getcwd(), "input"))
    os.mkdir(os.path.join(os.getcwd(), "genomes"))
    os.mkdir(os.path.join(os.getcwd(), "output"))
    os.mkdir(os.path.join(os.getcwd(), "slurmout"))
    os.mkdir(os.path.join(os.getcwd(), "slurm_job_status"))

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="init",
        description="Initialize pipeline",
        help="Initialize pipeline"
    )
    parser.set_defaults(func=init)