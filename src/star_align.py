import os
from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.scratch:
    slurm = SlurmJob(snakemake.params.scratch)
    slurm.setUp()
    sample = snakemake.wildcards.sample
    output_prefix = f'{slurm.scratch}/{sample}/'
    os.mkdir(os.path.join(slurm.scratch, sample))
else:
    output_prefix = snakemake.params.prefix

extra_args = snakemake.params.get('extra')
if extra_args is None:
    extra_args = ""

shell(f"""
STAR \\
    --runMode alignReads \\
    --readFilesIn {snakemake.input.sample} \\
    --genomeDir {snakemake.input.genome} \\
    --outFileNamePrefix {output_prefix} \\
    --runThreadN {snakemake.threads} \\
    --outSAMtype BAM Unsorted \\
    --readFilesCommand  zcat \\
    --outReadsUnmapped Fastx \\
    --outFilterMatchNminOverLread 0.9 \\
    --outFilterMatchNmin  18 \\
    --outFilterMismatchNmax 1  \\
    --outFilterMultimapNmax  1000000 \\
    --alignIntronMax  1 \\
    --alignIntronMin  2  {extra_args}
""")

if snakemake.params.scratch:
    shell(f"""
    gzip -c {output_prefix}Unmapped.out.mate1 \\
            > {output_prefix}Unmapped.fastq.gz
    mv {output_prefix}Aligned.out.bam {snakemake.output.bam}
    mv {output_prefix}Unmapped.fastq.gz {snakemake.output.unmapped}
    mv {output_prefix}Log.final.out {snakemake.params.prefix}Log.final.out
    """)
    slurm.tearDown()
else:
    shell(f"""
    gzip -c {snakemake.params.prefix}Unmapped.out.mate1 \\
            > {snakemake.params.prefix}Unmapped.fastq.gz
    rm {snakemake.params.prefix}Unmapped.out.mate1
    """)