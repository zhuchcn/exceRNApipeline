rule rrna_make_sequence:
    input: config['rRNA']
    output: "genomes/RiboRNA_homo_sapiens.fasta"
    params: 
        rna5s1=config['urls']['rna5s1']
    threads: 2
    shell: """
    curl \
        -o genomes/rna5s1.fasta \
        -H 'Content-type:text/x-fasta' \
        {params.rna5s1}
    cat genomes/rna5s1.fasta {input} \
        > {output}
    """

rule rrna_index:
    input:
        genome = "genomes/RiboRNA_homo_sapiens.fasta"
    output: 
        temp(directory("genomes/star_index_rrna_human"))
    params:
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        scratch = config.get('scratch'),
        mem = config.get("rrna_index_mem_gb")
    threads: config["rrna_index_cpu"]
    script: '../src/star_index.py'

rule rrna_mapping:
    input:
        genome = "genomes/star_index_rrna_human",
        sample = "output/02-UniVec/{sample}/Unmapped.fastq.gz"
    output:
        bam=temp("output/03-rRNA/{sample}/Aligned.out.bam"),
        unmapped=temp("output/03-rRNA/{sample}/Unmapped.fastq.gz")
    params:
        prefix="output/03-rRNA/{sample}/",
        scratch = config.get('scratch')
    threads: config["rrna_align_cpu"]
    script: '../src/star_align.py'