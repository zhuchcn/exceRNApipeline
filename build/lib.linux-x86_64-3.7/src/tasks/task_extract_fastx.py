from src.utils import logger
from Bio import SeqIO
import gzip
import argparse
    

def extract_fastx(input_path, output_path, seq_format, verbose, namelist_file):
    ih = gzip.open(input_path, 'rt')
    
    names = set()
    with open(namelist_file, 'rt') as fh:
        for l in fh:
            name = l.rstrip()
            if not name in names:
                names.add(name)
    if verbose:
        logger('Reading namelist file finished.')

    if verbose:
        logger('Start writing sequences.')
    
    oh = gzip.open(output_path, 'wt')

    for record in SeqIO.parse(ih, seq_format):
        record_name = record.description
        save = record_name in names
        if save:
            if verbose:
                logger(f'Saving {record_name}')
            SeqIO.write(record, oh, seq_format)

    oh.close()
    ih.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', type=str)
    parser.add_argument('-o', '--output-file', type=str)
    parser.add_argument('-l', '--namelist-file', type=str)
    parser.add_argument('-f', '--seq-format', type=str)
    parser.add_argument('-v', '--verbose', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    extract_fastx(args.input_file, args.output_file, args.seq_format,
                  args.verbose, args.namelist_file)


if __name__ == '__main__':
    main()