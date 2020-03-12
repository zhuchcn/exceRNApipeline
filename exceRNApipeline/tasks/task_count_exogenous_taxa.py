from _exceRNApipeline_taxaCounter import taxaCounter
from exceRNApipeline.includes.utils import logger
import os
import argparse
from snakemake.shell import shell


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-file", type=str)
    parser.add_argument("-o", "--output-prefix", type=str)
    parser.add_argument("-d", "--nodes-dmp", type=str)
    parser.add_argument("-m", "--names-dmp", type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    if(not os.path.exists(args.input_file)):
        raise argparse.ArgumentError(f"--input-file {args.input_file} does not exist.")
    if(not os.path.exists(args.nodes_dmp)):
        raise argparse.ArgumentError(f"--nodes-dmp {args.nodes_dmp} does not exist.")
    if(not os.path.exists(args.names_dmp)):
        raise argparse.ArgumentError(f"--names-dmp {args.names_dmp} does not exist.")
    if(not os.path.isdir(os.path.dirname(args.output_prefix))):
        msg = f"--output-prefix: you parsed the argument as {args.output_prefix}, " \
            + f"but the directory {os.path.dirname(args.output_prefix)} does not exist."
        raise argparse.ArgumentError(msg)
    taxaCounter(args.input_file, args.output_prefix, args.names_dmp, args.nodes_dmp)

if __name__ == '__main__':
    main()