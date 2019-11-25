rule download_mirbase:
    output: "genomes/miRBase_hairpin_v22.fa"
    params: url = config['urls']['miRBase']
    threads: 1
    shell: """
    wget {params.url} -O {output}.gz
    zcat {output}.gz > {output}
    """

rule mirbase_index:
    input: "genomes/miRBase_hairpin_v22.fa"
    output: 
        directory("genomes/star_index_mirbase")
    params:
        extra = "--genomeSAindexNbases 10 --genomeChrBinNbits 7",
        scratch = config.get('scratch'),
        mem = config.get("mirbase_index_mem_gb")
    threads: config["mirbase_index_cpu"]
    script: '../src/star_index.py'

rule mirbase_mapping:
    input:
        genome="genomes/star_index_mirbase",
        sample="output/05-RE/{sample}/Unmapped.fastq.gz"
    output:
        bam="output/06-miRBase/{sample}/Aligned.out.bam",
        unmapped="output/06-miRBase/{sample}/Unmapped.fastq.gz"
    params:
        prefix="output/06-miRBase/{sample}/",
        scratch = config.get('scratch'),
        extra="""
        --outFilterMismatchNoverLmax 0.3 \
        --outFilterMatchNminOverLread 1.0 \
        --outFilterMismatchNmax 0 \
        --outFilterMultimapNmax  10000 \
        --alignEndsType  EndToEnd
        """
    threads: config["mirbase_align_cpu"]
    script: '../src/star_align.py'