import gzip
from Bio import SeqIO
import argparse


def dump_taxmap(filepath):
    taxmap = {}
    with open(filepath, 'rt') as fh:
        for l in fh:
            l = l.rstrip().split('\t')
            silva_id = l[0]
            ncbi_tax = l[4]
            if silva_id not in taxmap:
                taxmap[silva_id] = ncbi_tax
    return taxmap


def solve_taxa(inpath, outpath, su, taxmap):
    with gzip.open(inpath, 'rt') as ih, gzip.open(outpath, 'wt') as oh:
        for record in SeqIO.parse(ih, 'fasta'):
            silva_id = record.id.split('.')[0]
            silva_id_long = record.id.split(' ')[0]
            ncbi_tax = taxmap[silva_id] + ':' + su.upper()
            record_id = silva_id_long + ":" + ncbi_tax
            record.id = record_id
            record.name = record_id
            record.description = record_id
            SeqIO.write(record, oh, 'fasta')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--fasta', type=str)
    parser.add_argument('-t', '--taxmap', type=str)
    parser.add_argument('-s', '--subunit', type=str)
    parser.add_argument('-o', '--output', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    taxmap = dump_taxmap(args.taxmap)
    solve_taxa(args.fasta, args.output, args.subunit, taxmap)

if __name__ == "__main__":
    main()
