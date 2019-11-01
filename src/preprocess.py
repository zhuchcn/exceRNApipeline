from snakemake.shell import shell
from utils import SlurmJob


if snakemake.params.get('useScratch'):
    slurm = SlurmJob()
    slurm.setUp()
    sample = snakemake.wildcards.sample
    prefix = f'{slurm.scratch}/{sample}'
    log = f'{slurm.scratch}/{sample}.htsStats.log'
else:
    log = snakemake.log
    prefix = snakemake.params.prefix


cmd = f"""
    {snakemake.params.hts_Stats} -L {log} -U {snakemake.input} | \
    {snakemake.params.hts_AdapterTrimmer} -A -L {log} -a {snakemake.params.adapter} | \
    {snakemake.params.hts_QWindowTrim} -n -A -L {log} | \
    {snakemake.params.hts_NTrimmer} -n -A -L {log}  | \
    {snakemake.params.hts_Stats} -A -L {log} -f {prefix}
"""
shell(cmd)

if snakemake.params.get('useScratch'):
    cmd = f"""
    mv {slurm.scratch}/{sample}_SE.fastq.gz  {snakemake.output}
    mv {slurm.scratch}/{sample}.htsStats.log  {snakemake.log}
    """
    shell(cmd)
    slurm.tearDown()
