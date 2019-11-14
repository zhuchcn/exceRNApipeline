from snakemake.shell import shell
from utils import SlurmJob
import os


if snakemake.params.get('useScratch'):
    slurm = SlurmJob()
    slurm.setUp()
    sample = snakemake.wildcards.sample
    namelist = f"{slurm.scratch}/aligned_reads.txt"
    output = f"{slurm.scratch}/{os.path.basename(snakemake.output.unmapped)}"
else:
    namelist = snakemake.params.namelist_path
    output = os.path.basenaem(snakemake.output.unmapped)

cmd = f"""
    zcat {snakemake.input.aligned} | cut -f 1 | sort | uniq >> {namelist}
    {snakemake.params.hts_fastx} extract-fastx \\
        --input-file {snakemake.aligned} \\
        --output-file {output} \\
        --namelist-file $namelist_file \\
        -v -u
"""

if snakemake.params.get('useScratch'):
    shell(f"mv {output} {snakemake.output.unmapped}")
    slurm.tearDown()
else:
    shell(f"rm {namelist}")