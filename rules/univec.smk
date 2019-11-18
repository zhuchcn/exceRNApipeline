rule univec_get_sequence:
    output: "genomes/UniVec_Core.fasta"
    params: url = config['urls']['univec']
    threads: 2
    shell: """
    wget {params.url} -O {output}
    """

rule univec_index:
    input:
        genome = config['genomes'].get('univec') or \
            "genomes/UniVec_Core.fasta"
    output: 
        directory("genomes/star_index_univec")
    params:
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        use_scratch = config['use_scratch']
    threads: 8
    script: '../src/star_index.py'

rule univec_mapping:
    input:
        genome="genomes/star_index_univec",
        sample="output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    output:
        bam="output/02-UniVec/{sample}/Aligned.out.bam",
        unmapped="output/02-UniVec/{sample}/Unmapped.fastq.gz"
    params:
        prefix="output/02-UniVec/{sample}/",
        use_scratch = config['use_scratch']
    threads: 12
    script: '../src/star_align.py'