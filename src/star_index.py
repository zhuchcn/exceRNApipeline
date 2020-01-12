import os
from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.scratch:
    slurm = SlurmJob(snakemake.params.scratch)
    slurm.setUp()
    wd = slurm.scratch
else:
    wd = os.path.dirname(snakemake.output[0])

genome_dir = os.path.basename(snakemake.output[0])
if snakemake.input.genome.startswith('/'):
    input_genome = snakemake.input.genome
else:
    input_genome = os.path.join(os.getcwd(), snakemake.input.genome)

mem = snakemake.params.get("mem") or snakemake.threads * 2
mem = mem * 1024 ** 3

cmd = f"""
    cd {wd}
    mkdir {genome_dir}
    STAR {snakemake.params.extra} \\
        --runMode genomeGenerate \\
        --genomeFastaFiles {input_genome} \\
        --genomeDir {genome_dir} \\
        --runThreadN {snakemake.threads} \\
        --limitGenomeGenerateRAM {mem}
"""
print(cmd, flush=True)
shell(cmd)

if snakemake.params.scratch:
    cmd = f"mv {slurm.scratch}/{genome_dir} {snakemake.output}"
    print(cmd, flush=True)
    shell(cmd)
    slurm.tearDown()
