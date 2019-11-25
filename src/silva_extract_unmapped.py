from snakemake.shell import shell
from utils import SlurmJob
import os


if snakemake.params.scratch:
    slurm = SlurmJob(snakemake.params.scratch)
    slurm.setUp()
    sample = snakemake.wildcards.sample
    namelist = f"{slurm.scratch}/aligned_reads.txt"
    output = f"{slurm.scratch}/{os.path.basename(snakemake.output.unmapped)}"
else:
    namelist = snakemake.params.namelist_path
    output = os.path.basenaem(snakemake.output.unmapped)

shell(f"""
    zcat {snakemake.input.aligned} |\\
        cut -f 1 | sort | uniq >> {namelist}
    hts_fastx extract-fastx \\
        --input-file {snakemake.input.fastq} \\
        --output-file {output} \\
        --namelist-file {namelist} \\
        -v -u
""")

if snakemake.params.scratch:
    shell(f"mv {output} {snakemake.output.unmapped}")
    slurm.tearDown()
else:
    shell(f"rm {namelist}")