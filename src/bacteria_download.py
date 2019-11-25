import os
from snakemake.shell import shell
from ensembl import getBacteriaGenomes
from utils import SlurmJob


if snakemake.params.scratch:
    slurm = SlurmJob()
    slurm.setUp()
    output_dir = os.path.join(slurm.scratch, "bacteria_genomes")
else:
    output_dir = snakemake.params.output_dir

getBacteriaGenomes(output_dir, snakemake.params.version,
                   snakemake.threads, verbose=True)

if snakemake.params.scratch:
    shell(f"mv {output_dir}/* {snakemake.params.output_dir}/")