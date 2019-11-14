import gzip
from Bio import SeqIO

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


def solva_taxa(inpath, outpath, su, taxmap):
    with gzip.open(inpath, 'rt') as ih, gzip.open(outpath, 'wt') as oh:
        for record in SeqIO.parse(ih, 'fasta'):
            silva_id = record.id.split('.')[0]
            ncbi_tax = taxmap[silva_id] + ':' + su.upper()
            record.id = ncbi_tax
            record.name = ncbi_tax
            record.description = ncbi_tax
            SeqIO.write(record, oh, 'fasta')


if __name__ == "__main__":
    taxmap = dump_taxmap(snakemake.input.taxmap)
    solva_taxa(
        snakemake.input.fasta,
        snakemake.output[0],
        snakemake.wildcards.su,
        taxmap
    )

