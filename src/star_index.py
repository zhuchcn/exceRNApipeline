import os
from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.get('use_scratch'):
    slurm = SlurmJob()
    slurm.setUp()
    wd = slurm.scratch
else:
    wd = os.path.dirname(snakemake.output[0])

genome_dir = os.path.basename(snakemake.output[0])
if snakemake.input.genome.startswith('/'):
    input_genome = snakemake.input.genome
else:
    input_genome = os.path.join(os.getcwd(), snakemake.input.genome)
cmd = f"""
    cd {wd}
    mkdir {genome_dir}
    STAR {snakemake.params.extra} \\
        --runMode genomeGenerate \\
        --genomeFastaFiles {input_genome} \\
        --genomeDir {genome_dir} \\
        --runThreadN {snakemake.threads} 
"""
shell(cmd)

if snakemake.params.get('use_scratch'):
    shell(f"mv {slurm.scratch}/{genome_dir} {snakemake.output}")
    slurm.tearDown()
