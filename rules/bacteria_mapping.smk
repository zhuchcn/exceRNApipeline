import os
import re
from src.ensembl import  EnsemblFTP


BACTERIA_COLLECTION_INDICES = EnsemblFTP().get_bacteria_collection_numbers()

rule bacteria_download_genomes:
    output: 
        # FIXME: The bacteria genomes are not marked as temp for now, but
        # should be removed in the end of the day.
        directory(expand(
            "genomes/bacteria/fasta/bacteria_{bacteria_ind}_collection",
            bacteria_ind=BACTERIA_COLLECTION_INDICES
        ))
    params:
        version = config.get("version") or "current",
        output_dir = "genomes/bacteria/fasta",
        scratch = config.get('scratch')
    threads: config["bacteria_download_cpu"]
    script: '../src/bacteria_download.py'

def bacteria_index_prev(wildcards):
    ind = int(wildcards.bacteria_ind)
    jobs = config["bacteria_align_parallel"]
    indices = BACTERIA_COLLECTION_INDICES
    inds_each_job = [indices[start::jobs] for start in range(jobs)]
    cuts = [len(x) for x in inds_each_job]
    cuts = [sum(cuts[:i]) for i in range(len(cuts))]
    if ind in cuts:
        return f"genomes/bacteria/fasta/bacteria_{ind}_collection"
    return expand(
        f"output/07-Bacteria/{{sample}}/Aligned_{ind-1}.txt.gz",
        sample=config["samples"]
    )

rule bacteria_combine_genome:
    input:
        genome_collection = "genomes/bacteria/fasta/bacteria_{bacteria_ind}_collection",
        prev = bacteria_index_prev
    output:
        temp("genomes/bacteria/fasta/bacteria_{bacteria_ind}_collection.fasta")
    threads: 2
    shell: "zcat {input.genome_collection}/*.fa.gz > {output}"

rule bacteria_index:
    input:
        genome = "genomes/bacteria/fasta/bacteria_{bacteria_ind}_collection.fasta" 
    output: 
        temp(directory("genomes/bacteria/star_index/bacteria_{bacteria_ind}_collection"))
    params:
        extra = "--genomeSAindexNbases 13 --genomeChrBinNbits 16",
        scratch = config.get('scratch'),
        mem = config["bacteria_index_ram_gb"]
    threads: config["bacteria_index_cpu"]
    script: "../src/star_index.py"

rule bacteria_mapping:
    input:
        genome = "genomes/bacteria/star_index/bacteria_{bacteria_ind}_collection",
        sample = "output/06-SILVA/{sample}/Unmapped.fastq.gz"
    output: 
        temp("output/07-Bacteria/{sample}/Aligned_{bacteria_ind}.txt.gz")
    params:
        prefix = lambda wildcards: f"06-SILVA/{wildcards.sample}/Aligned_{wildcards.bacteria_ind}/",
        scratch = config.get('scratch')
    threads: config["bacteria_align_cpu"]
    script: "../src/star_align_bacteria.py"

rule bacteria_combine_result:
    input: 
        expand(
            "output/07-Bacteria/{{sample}}/Aligned_{bacteria_ind}.txt.gz",
            bacteria_ind = BACTERIA_COLLECTION_INDICES
        )
    output: temp("output/07-Bacteria/{sample}/Aligned.txt.gz")
    threads: 2
    shell: """
    zcat {input} | \\
    awk -F '\t' '{{split($2, taxa, ":"); print($1"\t"taxa[1])}}' | \\
    uniq | gzip \\
    > {output}
    """

rule bacteria_count_taxa:
    input:
        sample = "output/07-Bacteria/{sample}/Aligned.txt.gz",
        taxdump = "genomes/taxdump"
    output:
        expand(
            "output/07-Bacteria/{{sample}}/exogenousAligned_bacteria_{tax}.txt",
            tax=config["tax_levels"]
        )
    params:
        prefix = "output/07-Bacteria/{sample}/exogenousAligned_bacteria_"
    threads: 4
    shell: """
    hts_taxonomy count-taxa \\
        --input-file {input.sample} \\
        --output-prefix {params.prefix} \\
        --nodes-dump  {input.taxdump}/nodes.dmp \\
        --names-dump {input.taxdump}/names.dmp
    """

rule bacteria_count_combine:
    input: 
        expand(
            "output/07-Bacteria/{sample}/exogenousAligned_bacteria_{{tax}}.txt",
            sample=config["samples"]
        )
    output: "output/07-Bacteria/bacteria_count_{tax}.txt"
    threads: 4
    run: 
        from pandas import read_csv

        counts = None
        for sample in config["samples"].keys():
            path = 'output/07-Bacteria/' + sample + '/exogenousAligned_bacteria_' + \
                   wildcards.tax + '.txt'
            new_counts = read_csv(
                path, sep='\t', index_col=0, names=[sample]
            )
            if counts is None:
                counts = new_counts
            else:
                counts = counts.join(new_counts, how='outer')
        counts.to_csv(output[0], sep='\t')


rule organize_exogenous:
    input:
        silva = expand("output/06-SILVA/SILVA_count_{tax}.txt",tax=config["tax_levels"]),
        bacteria = expand("output/07-Bacteria/bacteria_count_{tax}.txt",tax=config["tax_levels"])
    output:
        silva = directory("output/results/silva/"),
        bacteria = directory("output/results/bacteria/")
    threads: 1
    shell: """
    cp -r {input.silva} {output.silva}
    cp -r {input.bacteria} {output.bacteria} 
    """
    