from exceRNApipeline.includes.ensembl import getBacteriaGenomes
from .slurm_job import SlurmJob
from exceRNApipeline.includes.utils import logger
from snakemake.shell import shell
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir', type=str)
    parser.add_argument('-v', '--version', type=str)
    parser.add_argument('-s', '--scratch-dir', type=str)
    parser.add_argument('-t', '--nthreads', type=int)
    return parser.parse_args()

def main():
    args = parse_args()

    if args.scratch_dir == "None":
        args.scratch_dir = None

    if args.scratch_dir:
        with SlurmJob(args.scratch_dir) as slurm:
            output_dir = os.path.join(slurm.scratch, "bacteria_genomes")
            getBacteriaGenomes(output_dir, args.version, args.nthreads, 
                               verbose=True)
            cmd = f"mv {output_dir}/* {args.output_dir}/"
            logger(cmd)
            shell(cmd)
    else:
        getBacteriaGenomes(args.output_dir, args.verbose, args.nthreads,
                           verbose=True)
