import argparse
from snakemake.shell import shell
from .slurm_job import SlurmJob
from exceRNApipeline.includes.utils import logger


def pre_process(input_fq, adapter, log_file, prefix):
    cmd = f"""
    hts_Stats -L {log_file} -U {input_fq} | \\
    hts_AdapterTrimmer -A -L {log_file} -a {adapter} | \\
    hts_QWindowTrim -n -A -L {log_file} | \\
    hts_NTrimmer -n -A -L {log_file}  | \\
    hts_Stats -A -L {log_file} -f {prefix}
    """
    logger(cmd)
    shell(cmd)

def parse_args():
    parser = argparse.ArgumentParser(
        description="[exRNA-pipeline] pre-processing"
    )
    parser.add_argument("-i", "--input-fq", type=str,
                        help="Path to the input fastq files.")
    parser.add_argument("-o", "--output-fq", type=str,
                        help="Path to t he output fastq files.")
    parser.add_argument("-n", "--sample-name", type=str,
                        help="Sample name")
    parser.add_argument("-a", "--adapter", type=str,
                        help="Adapter sequence.")
    parser.add_argument("-l", "--log-file", type=str,
                        help="Path to the log file.")
    parser.add_argument("-p", "--prefix", type=str,
                        help="Output prefix")
    parser.add_argument("-s", "--scratch-dir", type=str,
                        help="Path to the scratch diractory.")
    args = parser.parse_args()
    if args.scratch_dir == "None":
        args.scratch_dir = None
    return args

def main():
    args = parse_args()
    if args.scratch_dir:
        with SlurmJob(args.scratch_dir) as slurm:
            pre_process(
                args.input_fq, args.adapter,
                f"{slurm.scratch}/{args.sample_name}.htsStats.log",
                f"{slurm.scratch}/{args.sample_name}"
            )
            cmd = f"""
            mv {slurm.scratch}/{args.sample_name}_SE.fastq.gz {args.output_fq}
            mv {slurm.scratch}/{args.sample_name}.htsStats.log {args.log_file}
            """
            logger(cmd)
            shell(cmd)
    else:
        pre_process(args.input_fq, args.adapter,
                   args.log_file, args. prefix)

if __name__ == "__main__":
    main()
    