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
    htseq-count -f bam -s no -i gene_id \\
        --additional-attr gene_name --additional-attr gene_type \\
        {input.bam} {input.gencode} > {output.gencode}
    htseq-count -f bam -s no -i gene_id \\
        --additional-attr gene_name --additional-attr gene_type \\
        {input.bam} {input.tRNA} > {output.tRNA};
    htseq-count -f bam -s no -i gene_id \\
        --additional-attr gene_name --additional-attr gene_type \\
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

rule organize_endogenous_counts:
    input:
        gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
        tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
        piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
        summary = "output/04-Genome/ReadsPerGene_summary.txt"
    output:
        gencode = "output/results/endogenous/ReadsPerGene_gencode.txt",
        tRNA = "output/results/endogenous/ReadsPerGene_tRNA.txt",
        piRNA = "output/results/endogenous/ReadsPerGene_piRNA.txt",
        summary = "output/results/endogenous/ReadsPerGene_summary.txt"
    threads: 1
    shell: """
    cp {input.gencode} {output.gencode}
    cp {input.tRNA} {output.tRNA}
    cp {input.piRNA} {output.piRNA}
    cp {input.summary} {output.summary}
    """

rule endogenous_count_barplot:
    input: "output/results/endogenous/ReadsPerGene_summary.txt"
    output: 
        report("output/results/endogenous/summary.png",
               caption="../reports/endogenous_summary.rst",
               category="Endogenous Genome")
    threads: 1
    script: "../src/endogenous_summary_plot.py"