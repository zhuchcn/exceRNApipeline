rule preprocess:
    input:
        sample=lambda wildcards: config['samples'][wildcards.sample]
    output: 
        fastq = temp("output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"),
        log = "output/01-Preprocess/{sample}/{sample}.htsStats.log"
    params:
        adapter=config['adapter'],
        prefix="output/01-Preprocess/{sample}/{sample}",
        scratch = config.get('scratch')
    threads: 5
    script: '../src/preprocess.py'

rule preprocess_stats:
    input: 
        logs = expand("output/01-Preprocess/{sample}/{sample}.htsStats.log",
                      sample=config["samples"]),
    output: 
        path = "output/results/qc/preprocess.tsv",
        plot = report("output/results/qc/preprocess_barplot.png", 
                      caption="../reports/preprocess_barplot.rst",
                      category="Preprocess")
    params:
        samples = list(config["samples"].keys())
    threads: 1
    script: '../src/preprocess_stats.py'

rule fastqc:
    input:
        fastq = "output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    output:
        report("output/results/qc/{sample}_SE_fastqc/Images/per_base_quality.png",
               caption="../reports/qc_per_base_quality.rst",
               category="QC")
    params:
        outdir = "output/results/qc"
    threads: 1
    shell: """
    fastqc -o {params.outdir} {input.fastq}
    unzip {params.outdir}/{wildcards.sample}_SE_fastqc.zip \\
        -d {params.outdir}/
    rm -rf {params.outdir}/{wildcards.sample}_SE_fastqc.zip
    """
    
    