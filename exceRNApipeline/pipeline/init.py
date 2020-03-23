import os
import shutil


def init(args):
    smk_dir = os.path.join(os.path.dirname(__file__), "smk")
    copy(
        os.path.join(smk_dir, "pipeline_config.yml"),
        os.path.join(os.getcwd(), "pipeline_config.yml")
    )
    copy(
        os.path.join(smk_dir, "slurm_config.yml"),
        os.path.join(os.getcwd(), "slurm_config.yml")
    )
    
    mkdir(os.path.join(os.getcwd(), "input"))
    mkdir(os.path.join(os.getcwd(), "genomes"))
    mkdir(os.path.join(os.getcwd(), "output"))
    mkdir(os.path.join(os.getcwd(), "slurmout"))
    mkdir(os.path.join(os.getcwd(), "slurm_job_status"))

def mkdir(path):
    try:
        os.mkdir(path)
        print("directory created at: " + path)
    except FileExistsError:
        print(f"directory exists: '{path}'")

def copy(src, dst):
    if os.path.exists(dst):
        print(f"file already exists: '{dst}'")
        return
    shutil.copy2(src, dst)
    print("config file initialized at: " + dst)

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="init",
        description="Initialize pipeline",
        help="Initialize pipeline"
    )
    parser.set_defaults(func=init)