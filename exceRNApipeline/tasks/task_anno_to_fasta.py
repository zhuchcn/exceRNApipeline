import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from exceRNApipeline.includes.utils import logger
import argparse


def anno2fasta(anno, genome, output, col_chr, col_start, col_end, col_name, 
               skip=0, one_based=False):
    try:
        col_name = [int(i) for i in col_name]
    except ValueError:
        raise ValueError('chol_name must all be int')
    # Parse  the genome. This is going to use a lot of RAM
    logger('start reading genome')
    if os.path.splitext(genome)[1] == '.gz':
        with gzip.open(genome, 'rt') as fh:
            genome = SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))
    else:
        genome = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))

    if os.path.splitext(anno)[1] == '.gz':
        ih = gzip.open(anno, 'rt')
    else:
        ih = open(anno, 'rt')

    logger("start extracting sequences..")
    oh = open(output, 'a')
    i = 0
    for l in ih:
        if i < skip:
            i += 1
            continue
        l = l.rstrip().rsplit()
        chrom = l[col_chr - 1]
        start = int(l[col_start - 1])
        end = int(l[col_end - 1])
        name = ':'.join([l[i - 1] for i in col_name])
        if not one_based:
            start = start - 1
        if chrom not in genome:
            continue
        record = SeqRecord(
            genome[chrom].seq[start:end],
            id=name,
            name='',
            description=''
        )
        SeqIO.write(record, oh, 'fasta')

    oh.close()
    ih.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-a', '--anno-file', type=str, default=None,
        help='File path to the annotation file'
    )
    parser.add_argument(
        '-g', '--genome-file', type=str, default=None,
        help='File path to the genome file'
    )
    parser.add_argument(
        '-o', '--output-file', type=str, default=None,
        help='File path to the output fasta'
    )
    parser.add_argument(
        '-c', '--column-chr', type=int, default=None,
        help='Column index for chromosome in the annotation file.'
    )
    parser.add_argument(
        '-s', '--column-start', type=int, default=None,
        help='Column index for the start position in the annotation file.'
    )
    parser.add_argument(
        '-e', '--column-end', type=int, default=None,
        help='Column index for the end position in the annotation file.'
    )
    parser.add_argument(
        '-n', '--column-name', nargs="+", type=int,
        help='Column index for the name of the gene in the annotation file.'
    )
    parser.add_argument(
        '-k', '--skip-lines', type=int, default=0,
        help="Number of lines to skip."
    )
    parser.add_argument(
        '-b', '--one-based', action='store_true',
        help='Whether the annotation file postion starts are one based.'
    )
    return parser.parse_args()

def main():
    args = parse_args()
    anno2fasta(args.anno_file, args.genome_file, args.output_file, 
               args.column_chr, args.column_start, args.column_end,
               args.column_name, args.skip_lines, args.one_based)

if __name__ == '__main__':
    main()