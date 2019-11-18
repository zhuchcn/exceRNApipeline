rule preprocess:
    input:
        sample=lambda wildcards: config['samples'][wildcards.sample]
    output: "output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    params:
        adapter=config['adapter'],
        prefix="output/01-Preprocess/{sample}/{sample}",
        use_scratch = config['use_scratch']
    log: "output/01-Preprocess/{sample}/{sample}.htsStats.log"
    threads: 5
    script: '../src/preprocess.py'