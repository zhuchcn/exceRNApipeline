from src.utils import logger
import os
import argparse
from snakemake.shell import shell


def bed2gtf(input_file, output_file, source, feature, gene_type, keep_input):
    cmd = f"""
    {"zcat" if input.endswith("gz") else "cat"} {input_file} |\\
        awk -F '\t' '{{ \\
            start = $2 + 1 \\
            stop = $3 + 1 \\
            print $1"\\t{source}\\t{feature}\\t"start"\\t"stop"\\t"$5"\\t"$6"\\t.\\tgene_id \\""$4"\\"; transcript_id \\""$4"\\"; gene_type=\\"{gene_type}\\"; gene_name \\""$4"\\"; transcript_type \\""$4"\\"; transcript_name \\""$4"\\";" \\
        }}'  \\
        > {output_file}
    """
    if not keep_input:
        cmd += f"rm {input_file}"
    logger(cmd)
    shell(cmd)

def make_gencode_gtf(url, output_dir, output_file):
    # gencode annotation
    # The comprehensive annotation for chromosomes and scaffolds are used.
    logger("Download and preprocess the gencode annotation")
    cmd = f"""
    wget -q \\
        {url} \\
        -O {output_dir}/gencode_genome.gtf.gz
    zcat {output_dir}/gencode_genome.gtf.gz |\\
        grep -v '^#' \\
        >> {output_file}
    rm {output_dir}/gencode_genome.gtf.gz
    """
    logger(cmd)
    shell(cmd)

def make_tRNA_gtf(url, output_dir, output_file):
    # Download the tRNA annotation from gencode, and process it.
    logger("Download and preprocess the tRNA annotation")
    cmd = f"""
    # tRNA gtf from gencode
    wget -q {url} -O {output_dir}/tRNA.gtf.gz
    gunzip -c {output_dir}/tRNA.gtf.gz \\
        > {output_dir}/tRNA.tmp.gtf
    rm {output_dir}/tRNA.gtf.gz
    """
    logger(cmd)
    shell(cmd)

    with open(f'{output_dir}/tRNA.tmp.gtf', 'rt') as ih, \
            open(output_file, 'w') as oh:
        for l in ih:
            if l.startswith('#'):
                continue
            l = l.rstrip().split('\t')
            l[1] = 'GtRNAdb'
            l[2] = 'exon'
            attrs = {}
            for attr in l[-1].split('; '):
                key, val = attr.replace(';', '').split(' ')
                val = val.replace('\"', '')
                attrs[key] = val
            attrs['gene_id'] = f"tRNA-{attrs['gene_type']}-{attrs['gene_id']}"
            attrs['transcript_id'] = f"tRNA-{attrs['transcript_name']}-{attrs['transcript_id']}"
            attrs['gene_type'] = 'tRNA'
            attrs['gene_name'] = attrs['gene_type']
            l[-1] = ' '.join([key + ' \"' + val + '\";' for key, val in attrs.items()])
            oh.write('\t'.join(l) + '\n')
    os.remove(f'{output_dir}/tRNA.tmp.gtf')

def make_piRNA_gtf(url, output_dir, output_file):
    # Download piRNA annotation from piRBase and convert to gtf
    logger("Download and preprocess piRNA annotation")
    tmp_name = url.split("/")[-1]
    cmd = f"""
    wget -q {url} -O {output_dir}/{tmp_name}
    """
    print(cmd, flush=True)
    shell(cmd)
    if tmp_name.endswith("bed") or tmp_name.endswith("bed.gz"):
        bed2gtf(input_file=f"{output_dir}/{tmp_name}",
                output_file=output_file,
                source="piRNA",
                feature="exon",
                gene_type="piRNA",
                keep_input=False)
    elif tmp_name.endswith("gz"):
        cmd = f"gunzip -c {output_dir}/{tmp_name} > {output_file}"
        print(cmd, flush=True)
        shell(cmd)
    else:
        cmd = f"mv {output_dir}/{tmp_name} {output_file}"
        logger(cmd)
        shell(cmd)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--genecode-url', type=str)
    parser.add_argument('-r', '--tRNA-url', type=str)
    parser.add_argument('-i', '--piRNA-url', type=str)
    parser.add_argument('-g', '--genecode-gtf', type=str)
    parser.add_argument('-t', '--tRNA-gtf', type=str)
    parser.add_argument('-p', '--piRNA-gtf', type=str)
    parser.add_argument('-o', '--output-dir', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    make_gencode_gtf(args.genecode_url, args.output_dir, args.genecode_gtf)
    make_tRNA_gtf(args.tRNA_url, args.output_dir, args.tRNA_gtf)
    make_piRNA_gtf(args.piRNA_url, args.output_dir, args.piRNA_gtf)

if __name__ == "__main__":
    main()
