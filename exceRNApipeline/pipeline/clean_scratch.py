import subprocess as sp
import argparse
import os
import sys


def clean_scratch(args):
    if not os.path.exists(args.slurm_status) or \
        len(os.listdir(args.slurm_status)) == 0:
        print("no record found for scratch usage.")
        sys.exit(0)
    print("cleaning scratches..")
    for job in os.listdir(args.slurm_status):
        job_file = os.path.join(args.slurm_status, job)
        with open(job_file, "rt") as fh:
            hostname = fh.readlines()[0].rstrip()
        print(hostname + " " + job)
        sp.run(f"ssh {hostname} rm -rf {args.scratch_dir}/{job}".split())
        os.remove(job_file)

def parse_args(subparsers):
    parser = subparsers.add_parser(
        name="clean-scratch",
        description="clean scratch temporary files in each node",
        help="clean scratch temporary files in each node"
    )
    parser.set_defaults(func=clean_scratch)
    parser.add_argument(
        "--slurm-status", type=str, default="slurm_job_status",
        help="the directory which stores the slurm scratch usage records" +\
             " (default: slurm_job_status)"
    )
    parser.add_argument(
        "--scratch-dir", type=str, default="/scratch",
        help="the scratch path (default: /scratch)"
    )