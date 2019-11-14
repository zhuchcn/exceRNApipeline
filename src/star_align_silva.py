import os
from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.get('useScratch'):
    slurm = SlurmJob()
    slurm.setUp()
    sample = snakemake.wildcards.sample
    index = snakemake.wildcards.silva_ind
    output_prefix = f'{slurm.scratch}/{sample}/silva_{index}_'
    os.mkdir(os.path.join(slurm.scratch, sample))
else:
    output_prefix = snakemake.prams.prefix

extra_args = snakemake.params.get('extra')
if extra_args is None:
    extra_args = ""

print(extra_args)

shell(f"""
{snakemake.params.STAR} \\
    --runMode alignReads \\
    --readFilesIn {snakemake.input.sample} \\
    --genomeDir {snakemake.input.genome} \\
    --outFileNamePrefix {output_prefix} \\
    --runThreadN {snakemake.threads} \\
    --limitBAMsortRAM {snakemake.params.ram} \\
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
{snakemake.params.samtools} view \\
    {output_prefix}Aligned.out.bam |\\
    awk -F '\\t' '{{{{print($1\"\\t\"$3)}}}}' | uniq | gzip \\
    > {snakemake.output}
""")

if snakemake.params.get('useScratch'):
    slurm.tearDown()
else:
    shell(f"rm -rf {snakemake.params.prefix}*")