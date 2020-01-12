from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.scratch:
    slurm = SlurmJob(snakemake.params.scratch)
    slurm.setUp()
    sample = snakemake.wildcards.sample
    prefix = f'{slurm.scratch}/{sample}'
    log = f'{slurm.scratch}/{sample}.htsStats.log'
else:
    log = snakemake.output.log
    prefix = snakemake.params.prefix

cmd = f"""
    hts_Stats -L {log} -U {snakemake.input} | \\
    hts_AdapterTrimmer -A -L {log} -a {snakemake.params.adapter} | \\
    hts_QWindowTrim -n -A -L {log} | \\
    hts_NTrimmer -n -A -L {log}  | \\
    hts_Stats -A -L {log} -f {prefix}
"""
print(cmd, flush=True)
shell(cmd)

if snakemake.params.scratch:
    cmd = f"""
        mv {slurm.scratch}/{sample}_SE.fastq.gz  {snakemake.output.fastq}
        mv {slurm.scratch}/{sample}.htsStats.log  {snakemake.output.log}
    """
    print(cmd, flush=True)
    shell(cmd)
    slurm.tearDown()
