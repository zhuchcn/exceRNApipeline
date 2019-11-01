from snakemake.shell import shell
import os

# Download the tRNA annotation from gencode, and process it.
cmd = """
# tRNA gtf from gencode
wget {snakemake.params.tRNA} -O {snakemake.params.output_dir}/tRNA.gtf.gz
gunzip -c {snakemake.params.output_dir}/tRNA.gtf.gz \
    > {snakemake.params.output_dir}/tRNA.tmp.gtf
rm {snakemake.params.output_dir}/tRNA.gtf.gz
"""
shell(cmd)

with open(f'{snakemake.params.output_dir}/tRNA.tmp.gtf', 'rt') as ih, \
        open(snakemake.output.tRNA, 'w') as oh:
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
os.remove(f'{snakemake.params.output_dir}/tRNA.tmp.gtf')

# Download piRNA annotation from piRBase and convert to gtf
cmd = """
wget {snakemake.params.piRNA} -O {snakemake.params.output_dir}/piRNAs.bed.gz
zcat {snakemake.params.output_dir}/piRNAs.bed.gz |\
    awk -F '\t' '{{
        start = $2 + 1
        stop = $3 + 1
        print $1"\tpiRbase\texon\t"start"\t"stop"\t"$5"\t"$6"\t.\tgene_id \"piRNA-"$4"\"; transcript_id \"piRNA-"$4"\"; gene_type=\"piRNA\"; gene_name \""$4"\"; transcript_type \""$4"\"; transcript_name \""$4"\";"
    }}'  \
    > {snakemake.output.piRNA}
rm {snakemake.params.output_dir}/piRNAs.bed.gz
"""
shell(cmd)

# gencode annotation
# The comprehensive annotation for chromosomes and scaffolds are used.
cmd = """
wget {snakemake.params.gencode} -O {snakemake.params.output_dir}/gencode_genome.gtf.gz
zcat {snakemake.params.output_dir}/gencode_genome.gtf.gz |\
    grep -v '^#' \
    >> {snakemake.output.gencode}
rm {snakemake.params.output_dir}/gencode_genome.gtf.gz
"""
shell(cmd)