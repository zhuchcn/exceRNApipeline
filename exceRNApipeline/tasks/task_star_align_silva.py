from .slurm_job import SlurmJob
from exceRNApipeline.includes.utils import logger
import os
import argparse
from snakemake.shell import shell


def mkdir(path):
    if os.path.isdir(path):
        return
    parent = os.path.dirname(path)
    if parent != "":
        mkdir(parent)
    os.mkdir(path)

def star_align(input_fq, genome_index, output_prefix, output_txt, nthreads, extra_args):
    mkdir(os.path.dirname(output_prefix))
    cmd = f"""
    STAR \\
        --runMode alignReads \\
        --readFilesIn {input_fq} \\
        --genomeDir {genome_index} \\
        --outFileNamePrefix {output_prefix} \\
        --runThreadN {nthreads} \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes  Standard \\
        --readFilesCommand  zcat \\
        --alignEndsType EndToEnd \\
        --outFilterMatchNminOverLread 1.0 \\
        --outFilterMatchNmin  18 \\
        --outFilterMismatchNmax 0  \\
        --outFilterMismatchNoverLmax 0.3 \\
        --outFilterMultimapNmax  10000 \\
        --alignIntronMax  1 \\
        --alignIntronMin  2  {extra_args}
    """
    logger(cmd)
    shell(cmd)

    cmd = f"""
    samtools view \\
        {output_prefix}Aligned.out.bam |\\
    awk -F '\\t' '{{{{split($3, tax, ":"); print $1\"\\t\"tax[2]}}}}' | uniq | gzip \\
    > {output_txt}
    """
    logger(cmd)
    shell(cmd)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-fq', type=str)
    parser.add_argument('-n', '--sample-name', type=str)
    parser.add_argument('-c', '--silva-index', type=str)
    parser.add_argument('-g', '--genome-index', type=str)
    parser.add_argument('-o', '--output-txt', type=str)
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
        index = args.silva_index
        with SlurmJob(args.scratch_dir) as slurm:
            output_prefix = f'{slurm.scratch}/{name}/silva_{index}_'
            # create a temp dir in scratch
            cmd = f"mkdir {slurm.scratch}/{name}"
            logger(cmd)
            shell(cmd)

            star_align(args.input_fq, args.genome_index, output_prefix,
                       args.output_txt, args.nthreads, extra_args)

    else:
        star_align(args.input_fq, args.genome_index, args.output_prefix,
                   args.output_txt, args.nthreads, extra_args)
        cmd = f"rm -rf {args.output_prefix}*"
        print(cmd, flush=True)
        shell(cmd)

if __name__ == "__main__":
    main()