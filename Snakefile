shell.prefix("set +o pipefail; ")

configfile: "pipeline_config.yml"

# TODO: figure this out
# report: "report/workflow.rst"

rule all:
    input:
        gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
        tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
        piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
        summary = "output/04-Genome/ReadsPerGene_summary.txt"

# Step 1: preprocessing
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
    script: 'src/preprocess.py'

# Step 2: map to the UniVec library
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
        STAR = config['softwares']['STAR'],
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        useScratch = config['useScratch']
    threads: 8
    script: 'src/star_index.py'

rule univec_mapping:
    input:
        genome="genomes/star_index_univec",
        sample="output/01-Preprocess/{sample}/{sample}_SE.fastq.gz"
    output:
        bam="output/02-UniVec/{sample}/Aligned.out.bam",
        unmapped="output/02-UniVec/{sample}/Unmapped.fastq.gz"
    params:
        prefix="output/02-UniVec/{sample}/",
        STAR=config['softwares']['STAR'],
        useScratch = config['useScratch']
    threads: 12
    script: 'src/star_align.py'

# Step 3: map to the rRNA sequences.
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
        STAR=config['softwares']['STAR'],
        extra = "--genomeSAindexNbases 8 --genomeChrBinNbits 16",
        useScratch = config['useScratch']
    threads: 8
    script: 'src/star_index.py'

rule rrna_mapping:
    input:
        genome = "genomes/star_index_rrna_human",
        sample = "output/02-UniVec/{sample}/Unmapped.fastq.gz"
    output:
        bam="output/03-rRNA/{sample}/Aligned.out.bam",
        unmapped="output/03-rRNA/{sample}/Unmapped.fastq.gz"
    params:
        prefix="output/03-rRNA/{sample}/",
        STAR=config['softwares']['STAR'],
        useScratch = config['useScratch']
    threads: 12
    script: 'src/star_align.py'

# Step 4: map to the human genome. 
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
        genome=config["genomes"].get('human') or \
            "genomes/human_genome_primary_assembly.fa"
    output: 
        directory("genomes/star_index_human_genome")
    params:
        STAR=config['softwares']['STAR'],
        extra = "--genomeSAindexNbases 14 --genomeChrBinNbits 18",
        useScratch = config['useScratch']
    threads: 12
    script: 'src/star_index.py'

rule genome_mapping:
    input:
        genome = "genomes/star_index_human_genome",
        sample = "output/03-rRNA/{sample}/Unmapped.fastq.gz"
    output:
        bam="output/04-Genome/{sample}/Aligned.out.bam",
        unmapped="output/04-Genome/{sample}/Unmapped.fastq.gz"
    params:
        STAR=config['softwares']['STAR'],
        prefix="output/04-Genome/{sample}/",
        useScratch = config['useScratch']
    threads: 24
    script: 'src/star_align.py'

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
    script: "src/make_human_genome_gtf.py"
    
rule genome_count:
    input:
        gencode = config['annotations'].get('gencode') \
                or 'genomes/gencode_human_annotation.gtf',
        tRNA = config['annotations'].get('tRNA') \
                or 'genomes/tRNA.gtf',
        piRNA = config['annotations'].get('piRNA') \
                or 'genomes/piRNA.gtf',
        bam="output/04-Genome/{sample}/Aligned.out.bam"
    output: 
        gencode = "output/04-Genome/{sample}/ReadsPerGene_gencode.txt",
        tRNA = "output/04-Genome/{sample}/ReadsPerGene_tRNA.txt",
        piRNA = "output/04-Genome/{sample}/ReadsPerGene_piRNA.txt"
    params:
        htseq_count = config['softwares']['htseq_count']
    threads: 4
    shell: """
    {params.htseq_count} -f bam -s no -i gene_id \
        --additional-attr gene_name --additional-attr gene_type \
        {input.bam} {input.gencode} > {output.gencode}
    {params.htseq_count} -f bam -s no -i gene_id \
        --additional-attr gene_name --additional-attr gene_type \
        {input.bam} {input.tRNA} > {output.tRNA};
    {params.htseq_count} -f bam -s no -i gene_id \
        --additional-attr gene_name --additional-attr gene_type \
        {input.bam} {input.piRNA} > {output.piRNA};
    """

rule summarize_counts:
    input:
        gencode = expand("output/04-Genome/{sample}/ReadsPerGene_gencode.txt", sample=config["samples"]),
        tRNA = expand("output/04-Genome/{sample}/ReadsPerGene_tRNA.txt", sample=config["samples"]),
        piRNA = expand("output/04-Genome/{sample}/ReadsPerGene_piRNA.txt", sample=config["samples"])
    output:
        gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
        tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
        piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
        summary = "output/04-Genome/ReadsPerGene_summary.txt"
    params:
        samples = list(config['samples'].keys())
    threads: 4
    script: "src/summarize_counts_endogenous.py"
    