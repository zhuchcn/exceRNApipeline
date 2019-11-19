import os
import re
from src.ensembl import  EnsemblFTP


BACTERIA_COLLECTION_INDICES = EnsemblFTP().get_bacteria_collection_numbers()

rule bacteria_download_genomes:
    output: 
        temp(expand(
            "genomes/bacteria/fasta/bacteria_{bacteria_ind}_collection",
            bacteria_ind=BACTERIA_COLLECTION_INDICES
        ))
    params:
        version = config.get("version") or "current",
        output_dir = "genomes/bacteria/fasta"
    threads: 24
    script: '../src/bacteria_download.py'

def bacteria_index_prev(wildcards):
    ind = int(wildcards.bacteria_ind)
    jobs = config["parallel"]["bacteria"]
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
        temp("geomes/bacteria/fasta/bacteria_{bacteria_ind}_collection.fasta")
    threads: 2
    shell: "zcat {input.genome_collection}/*.fa.gz > {output}"

rule bacteria_index:
    input:
        genome = "geomes/bacteria/fasta/bacteria_{bacteria_ind}_collection.fasta" 
    output: 
        temp(directory("genomes/bacteria/star_index/bacteria_{bacteria_ind}_collection/"))
    params:
        extra = "--genomeSAindexNbases 13 --genomeChrBinNbits 16",
        use_scratch = config['use_scratch']
    threads: 12
    script: "../src/star_index.py"

rule bacteria_mapping:
    input:
        genome_index = "genomes/bacteria/star_index/bacteria_{bacteria_ind}_collection/",
        sample = "output/06-SILVA/{sample}/Unmapped.fastq.gz"
    output: 
        temp("output/07-Bacteria/{sample}/Aligned_{bacteria_ind}.txt.gz")
    params:
        prefix="06-SILVA/{sample}/Aligned_{silva_ind}/",
        use_scratch = config['use_scratch'],
        ram="34359738368"
    threads: 16
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
                   wildcards.tax + 'txt'
            new_counts = read_csv(
                path, sep='\t', index_col=0, names=[sample]
            )
            if counts is None:
                counts = new_counts
            else:
                counts = counts.join(new_counts, how='outer')
        counts.to_csv(output, sep='\t')

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