import os
from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.scratch:
    slurm = SlurmJob(snakemake.params.scratch)
    slurm.setUp()
    sample = snakemake.wildcards.sample
    index = snakemake.wildcards.bacteria_ind
    output_prefix = f'{slurm.scratch}/{sample}/bacteria_{index}_'
    os.mkdir(os.path.join(slurm.scratch, sample))
else:
    output_prefix = snakemake.params.prefix

extra_args = snakemake.params.get('extra')
if extra_args is None:
    extra_args = ""

print(extra_args)

shell(f"""
STAR \\
    --runMode alignReads \\
    --readFilesIn {snakemake.input.sample} \\
    --genomeDir {snakemake.input.genome} \\
    --outFileNamePrefix {output_prefix} \\
    --runThreadN {snakemake.threads} \\
    --outSAMtype BAM Unsorted \\
    --outSAMattributes  Standard \\
    --readFilesCommand  zcat \\
    --alignEndsType EndToEnd \\
    --outFilterMatchNminOverLread 1.0 \\
    --outFilterMatchNmin  18 \\
    --outFilterMismatchNmax 0  \\
    --outFilterMismatchNoverLmax 0.3 \\
    --outFilterMultimapNmax  10000 \\
    --alignIntronMax  1 \\
    --alignIntronMin  2  {extra_args}
""")

shell(f"""
samtools view \\
    {output_prefix}Aligned.out.bam |\\
    awk -F '\\t' '{{{{print($1\"\\t\"$3)}}}}' | uniq | gzip \\
    > {snakemake.output}
""")

if snakemake.params.scratch:
    slurm.tearDown()
else:
    shell(f"rm -rf {snakemake.params.prefix}*")