rule rrna_make_sequence:
    input: U13369=config['genomes']['U13369']
    output: "genomes/RiboRNA_homo_sapiens.fasta"
    params: 
        rna5s1=config['urls']['rna5s1']
    threads: 2
    shell: """
    curl \
        -o genomes/rna5s1.fasta \
        -H 'Content-type:text/x-fasta' \
        {params.rna5s1}
    cat genomes/rna5s1.fasta {input.U13369} \
        > {output}
    """

rule rrna_index:
    input:
        genome = config['genomes'].get('rRNA') or \
                 "genomes/RiboRNA_homo_sapiens.fasta"
    output: 
        directory("genomes/star_index_rrna_human")
    params:
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        useScratch = config['useScratch']
    threads: 8
    script: '../src/star_index.py'

rule rrna_mapping:
    input:
        genome = "genomes/star_index_rrna_human",
        sample = "output/02-UniVec/{sample}/Unmapped.fastq.gz"
    output:
        bam="output/03-rRNA/{sample}/Aligned.out.bam",
        unmapped="output/03-rRNA/{sample}/Unmapped.fastq.gz"
    params:
        prefix="output/03-rRNA/{sample}/",
        useScratch = config['useScratch']
    threads: 12
    script: '../src/star_align.py'