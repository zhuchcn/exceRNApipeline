from .slurm_job import SlurmJob
from src.utils import logger
import os
import argparse
from snakemake.shell import shell


def star_align(input_fq, genome_index, output_prefix, nthreads, extra_args):
    cmd = f"""
    STAR \\
        --runMode alignReads \\
        --readFilesIn {input_fq} \\
        --genomeDir {genome_index} \\
        --outFileNamePrefix {output_prefix} \\
        --runThreadN {nthreads} \\
        --outSAMtype BAM Unsorted \\
        --readFilesCommand  zcat \\
        --outReadsUnmapped Fastx \\
        --outFilterMatchNminOverLread 0.9 \\
        --outFilterMatchNmin  18 \\
        --outFilterMismatchNmax 1  \\
        --outFilterMultimapNmax  1000000 \\
        --alignIntronMax  1 \\
        --alignIntronMin  2  {extra_args}
    """
    logger(cmd)
    shell(cmd)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-fq', type=str)
    parser.add_argument('-n', '--sample-name', type=str)
    parser.add_argument('-g', '--genome-index', type=str)
    parser.add_argument('-b', '--output-bam', type=str)
    parser.add_argument('-u', '--output-unmapped', type=str)
    parser.add_argument('-p', '--output-prefix', type=str)
    parser.add_argument('-t', '--nthreads', type=int)
    parser.add_argument('-s', '--scratch-dir', type=str)
    parser.add_argument('-a', '--extra-args', type=str)
    return parser.parse_args()
    

def main():
    args = parse_args()
    if args.scratch_dir == "None":
        args.scratch_dir = None
    
    extra_args = args.extra_args
    if extra_args == "None" or extra_args is None:
        extra_args = ""

    if args.scratch_dir:
        name = args.sample_name
        with SlurmJob(args.scratch_dir) as slurm:
            output_prefix = f'{slurm.scratch}/{name}'
            # create a temp dir in scratch
            cmd = f"mkdir {output_prefix}"
            logger(cmd)
            shell(cmd)

            star_align(args.input_fq, args.genome_index, output_prefix,
                       args.nthreads, extra_args)
            cmd = f"""
    gzip -c {output_prefix}Unmapped.out.mate1 \\
            > {output_prefix}Unmapped.fastq.gz
    mv {output_prefix}Aligned.out.bam {args.output_bam}
    mv {output_prefix}Unmapped.fastq.gz {args.output_unmapped}
    mv {output_prefix}Log.final.out {args.output_prefix}Log.final.out
    """
            logger(cmd)
            shell(cmd)
    else:
        star_align(args.input_fq, args.genome_index, args.output_prefix,
                   args.nthreads, extra_args)
        cmd = f"""
    gzip -c {args.output_prefix}Unmapped.out.mate1 \\
            > {args.output_prefix}Unmapped.fastq.gz
    rm {args.output_prefix}Unmapped.out.mate1
    """
        logger(cmd)
        shell(cmd)

if __name__ == "__main__":
    main()