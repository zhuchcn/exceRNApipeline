import os
from snakemake.shell import shell
from ensembl import getBacteriaGenomes
from utils import SlurmJob


if snakemake.params.get('use_scratch'):
    slurm = SlurmJob()
    slurm.setUp()
    output_dir = os.join(slurm.scratch, "bacteria_genomes")
else:
    output_dir = snakemake.params.output_dir

getBacteriaGenomes(output_dir, snakemake.params.version, snakemake.threads)

if snakemake.params.get('use_scratch'):
    shell(f"mv {output_dir} {sankemake.params.output_dir}")