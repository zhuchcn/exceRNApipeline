rule genome_count:
    input:
        gencode = 'genomes/gencode_human_annotation.gtf',
        tRNA = 'genomes/tRNA.gtf',
        piRNA = 'genomes/piRNA.gtf',
        bam="output/04-Genome/{sample}/Aligned.out.bam"
    output: 
        gencode = temp("output/04-Genome/{sample}/ReadsPerGene_gencode.txt"),
        tRNA = temp("output/04-Genome/{sample}/ReadsPerGene_tRNA.txt"),
        piRNA = temp("output/04-Genome/{sample}/ReadsPerGene_piRNA.txt")
    threads: 4
    shell: """
    htseq_count -f bam -s no -i gene_id \
        --additional-attr gene_name --additional-attr gene_type \
        {input.bam} {input.gencode} > {output.gencode}
    htseq_count -f bam -s no -i gene_id \
        --additional-attr gene_name --additional-attr gene_type \
        {input.bam} {input.tRNA} > {output.tRNA};
    htseq_count -f bam -s no -i gene_id \
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
    script: "../src/summarize_counts.py"

rule organize_endogenous:
    input:
        gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
        tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
        piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
        summary = "output/04-Genome/ReadsPerGene_summary.txt",
        bam = expand("output/04-Genome/{sample}/Aligned.out.bam", sample=config["samples"]),
    output:
        directory("output/results/endogenous/")
    threads: 1
    shell: """
    cp {input.gencode} {output}
    cp {input.tRNA} {output}
    cp {input.piRNA} {output}
    cp {input.summary} {output}
    cp {input.bam} {output}/bam/
    """