rule download_repeat_masker:
    output: temp("genomes/repeat_masker.fa.out")
    params: url = config['urls']['repeat_masker']
    threads: 1
    shell: """
    wget {params.url} -O {output}.gz
    zcat {output}.gz > {output}
    """

rule make_re_fasta:
    input: 
        annotation = "genomes/repeat_masker.fa.out",
        genome = "genomes/GRCh38.primary_assembly.genome.fa"
    output: "genomes/RepeatMasker.fa"
    threads: 24
    shell: """
    hts_fastx anno2fasta \\
        -a {input.annotation} -g {input.genome} -o {output} \\
        -c 5 -s 6 -e 7 -n 10 -n 11 -n 5 -n 6 -n 7 -k 3
    """

rule re_index:
    input: genome = "genomes/RepeatMasker.fa"
    output: 
        temp(directory("genomes/star_index_re"))
    params:
        extra = "--genomeSAindexNbases 14 --genomeChrBinNbits 8",
        use_scratch = config['use_scratch']
    threads: 24
    script: '../src/star_index.py'

rule re_mapping:
    input:
        genome="genomes/star_index_re",
        sample="output/04-Genome/{sample}/Unmapped.fastq.gz"
    output:
        bam=temp("output/05-RE/{sample}/Aligned.out.bam"),
        unmapped=temp("output/05-RE/{sample}/Unmapped.fastq.gz")
    params:
        prefix="output/05-RE/{sample}/",
        use_scratch = config['use_scratch']
    threads: 24
    script: '../src/star_align.py'