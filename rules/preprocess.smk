rule preprocess:
    input:
        sample=lambda wildcards: config['samples'][wildcards.sample]
    output: "output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    params:
        hts_Stats=config["softwares"]["hts_Stats"],
        hts_AdapterTrimmer=config["softwares"]["hts_AdapterTrimmer"],
        hts_QWindowTrim=config["softwares"]["hts_QWindowTrim"],
        hts_NTrimmer=config["softwares"]["hts_NTrimmer"],
        adapter=config['adapter'],
        prefix="output/01-Preprocess/{sample}/{sample}",
        useScratch = config['useScratch']
    log: "output/01-Preprocess/{sample}/{sample}.htsStats.log"
    threads: 5
    script: '../src/preprocess.py'