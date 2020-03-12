from exceRNApipeline.includes.utils import logger
from .slurm_job import SlurmJob
import os
import argparse
from snakemake.shell import shell


def index_genome(input_fas, wd, genome_dir, nthreads, mem, extra_args):
    cmd = f"""
    cd {wd}
    mkdir {genome_dir}
    STAR {extra_args} \\
        --runMode genomeGenerate \\
        --genomeFastaFiles {input_fas} \\
        --genomeDir {genome_dir} \\
        --runThreadN {nthreads} \\
        --limitGenomeGenerateRAM {mem}
    """
    logger(cmd)
    shell(cmd)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-fastas', type=str)
    # parser.add_argument('-g', '--genome-dir', type=str)
    # parser.add_argument('-w', '--working-dir', type=str)
    parser.add_argument('-o', '--output-dir', type=str)
    parser.add_argument('-t', '--nthreads', type=int)
    parser.add_argument('-m', '--mem-gb', type=int)
    parser.add_argument('-a', '--extra-args', type=str)
    parser.add_argument("-s", "--scratch-dir", type=str,
                        help="Path to the scratch diractory.")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    if args.scratch_dir == "None":
        args.scratch_dir = None
    
    nthreads = args.nthreads
    genome_dir = os.path.basename(args.output_dir)
    input_fastas = args.input_fastas if args.input_fastas.startswith("/") \
                else os.path.join(os.getcwd(), args.input_fastas)
    
    mem = args.mem_gb or nthreads * 2
    mem *= 1024 ** 3

    extra_args = args.extra_args 
    if extra_args == "None" or extra_args is None:
        extra_args = ""
    
    if args.scratch_dir:
        with SlurmJob(args.scratch_dir) as slurm:
            wd = slurm.scratch
            index_genome(input_fastas, wd, genome_dir, nthreads, mem, extra_args)
            # move genome index from scratch to the output dir
            cmd = f"mv {slurm.scratch}/{genome_dir} {args.output_dir}"
            logger(cmd)
            shell(cmd)
    else:
        wd = os.path.dirname(args.output_dir)
        index_genome(input_fastas, wd, genome_dir, nthreads, mem, extra_args)


if __name__ == "__main__":
    main()