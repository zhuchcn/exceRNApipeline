from _exRNApipeline_taxaCounter import taxaxCounter
from src.utils import logger
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
    taxaCounter(args.input_file, args.output_prefix, args.names_dmp, args.nodes_dmp)

if __names__ == '__main__':
    main()