rule genome_download:
    output: "genomes/human_genome_primary_assembly.fa"
    params: url = config['urls']['human']
    threads: 2
    shell: """
    wget {params.url} -O {output}.gz
    gunzip -c {output}.gz > {output}
    """

rule genome_index:
    input:
        genome = "genomes/human_genome_primary_assembly.fa"
    output: 
        temp(directory("genomes/star_index_human_genome"))
    params:
        extra = "--genomeSAindexNbases 14 --genomeChrBinNbits 18",
        use_scratch = config['use_scratch']
    threads: 12
    script: '../src/star_index.py'

rule genome_mapping:
    input:
        genome = "genomes/star_index_human_genome",
        sample = "output/03-rRNA/{sample}/Unmapped.fastq.gz"
    output:
        bam="output/04-Genome/{sample}/Aligned.out.bam",
        unmapped=temp("output/04-Genome/{sample}/Unmapped.fastq.gz")
    params:
        prefix="output/04-Genome/{sample}/",
        use_scratch = config['use_scratch'],
        ram = config["endogenous_align_ram"]
    threads: 24
    script: '../src/star_align.py'

rule make_genome_annotation:
    output: 
        gencode = 'genomes/gencode_human_annotation.gtf',
        tRNA = 'genomes/tRNA.gtf',
        piRNA = 'genomes/piRNA.gtf'
    params: 
        output_dir = 'genomes',
        tRNA = config['urls']['tRNA'],
        piRNA = config['urls']['piRNA'],
        gencode = config['urls']['gencode']
    threads: 2
    script: "../src/make_human_genome_gtf.py"