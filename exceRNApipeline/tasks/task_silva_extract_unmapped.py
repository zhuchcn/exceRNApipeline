from .task_extract_fastx import main as extract_fastx
from exceRNApipeline.includes.utils import logger
from .slurm_job import SlurmJob
from snakemake.shell import shell
import os
import argparse
import sys


def extract_silva_unmapped(aligned_txt, unmapped_fq, output_fq, namelist, verbose):
    cmd = f"""
    zcat {aligned_txt} |\\
        cut -f 1 | sort | uniq >> {namelist}
    """
    logger(cmd)
    shell(cmd)

    _sys_argv = sys.argv
    sys.argv = [
        "",
        "-i", unmapped_fq,
        "-o", output_fq,
        "-l", namelist,
        "-f", 'fastq'
    ]
    if verbose:
        sys.argv.append("-v")
    extract_fastx()
    sys.argv = _sys_argv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aligned-txt', type=str)
    parser.add_argument('-u', '--unmapped-fq', type=str)
    parser.add_argument('-o', '--output-fq', type=str)
    parser.add_argument('-l', '--namelist-txt', type=str)
    parser.add_argument('-s', '--scratch-dir', type=str)
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser.parse_args()

def main():
    args = parse_args()
    
    if args.scratch_dir == "None":
        args.scratch_dir = None
    
    namelist = args.namelist_txt
    output = args.output_fq

    if args.scratch_dir:
        with SlurmJob(args.scratch_dir) as slurm:
            namelist = os.path.basename(namelist)
            namelist = slurm.scratch + '/' + namelist
            output = os.path.basename(output)
            output = slurm.scratch + '/' + output
            extract_silva_unmapped(args.aligned_txt, args.unmapped_fq, output, 
                                   namelist, args.verbose)
            cmd = f"mv {output} {args.output_fq}"
            logger(cmd)
            shell(cmd)
    else:
        extract_silva_unmapped(args.aligned_txt, args.unmapped_fq, output, 
                               namelist, args.verbose)
        cmd = f"rm {namelist}"
        logger(cmd)
        shell(cmd)


if __name__ == '__main__':
    main()