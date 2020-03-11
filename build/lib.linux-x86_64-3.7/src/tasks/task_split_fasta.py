from Bio import SeqIO
import os
import gzip
from src.utils import logger
import math
import argparse


def open_file(file):
    if os.path.splitext(file)[1] == '.gz':
        fh = gzip.open(file, 'rt')
    else:
        fh = open(file, 'rt')
    return fh

def split_fa(input_file, output_prefix, n_record=None, n_batch=None):
    if n_record:
        if n_batch:
            logger("n_batch is ignored")        
    elif n_batch:
        n = 0
        fh = open_file(input_file)
        for line in fh:
            if line.startswith(">"):
                n += 1
        fh.close()
        n_record = math.ceil( n / n_batch )
    else:
        raise ValueError('At least one of n_record or n_batch must be given.')
    print(n_record)

    fh = open_file(input_file)
    
    split_fa_n_record(fh, output_prefix, n_record)
    fh.close()

def split_fa_n_record(ih, output_prefix, n_record):
    i = 0
    j = 1
    seqs = []
    def write():
        logger(f"Writing {len(seqs)} records to {output_prefix}{j}.fasta")
        with open(f"{output_prefix}{j}.fasta", 'w') as oh:
            SeqIO.write(seqs, oh, 'fasta')
    for record in SeqIO.parse(ih, 'fasta'):
        seqs.append(record)
        i += 1
        if i >= n_record:
            write()
            i = 0
            j += 1
            seqs = []
    if i != 0:
        write()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', type=str)
    parser.add_argument('-o', '--output-prefix', type=str)
    parser.add_argument('-r', '--n-record', type=int)
    parser.add_argument('-b', '--n-batch', type=int)
    return parser.parse_args()

def main():
    args = parse_args()
    split_fa(args.input_file, args.output_prefix, args.n_record, args.n_batch)


if __name__ == '__main__':
    main()