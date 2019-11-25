rule univec_get_sequence:
    output: "genomes/UniVec_Core.fasta"
    params: url = config['urls']['univec']
    threads: 2
    shell: """
    wget {params.url} -O {output}
    """

rule univec_index:
    input:
        genome = "genomes/UniVec_Core.fasta"
    output: 
        temp(directory("genomes/star_index_univec"))
    params:
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        scratch = config.get("scratch"),
        mem = config.get("univec_index_mem_gb")
    threads: config["univec_index_cpu"]
    script: '../src/star_index.py'

rule univec_mapping:
    input:
        genome="genomes/star_index_univec",
        sample="output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    output:
        bam=temp("output/02-UniVec/{sample}/Aligned.out.bam"),
        unmapped=temp("output/02-UniVec/{sample}/Unmapped.fastq.gz")
    params:
        prefix="output/02-UniVec/{sample}/",
        scratch = config.get("scratch")
    threads: config["univec_align_cpu"]
    script: '../src/star_align.py'