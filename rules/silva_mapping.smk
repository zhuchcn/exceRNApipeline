rule silva_mapping:
    input:
        sample = "output/05-RE/{sample}/Unmapped.fastq.gz",
        genome = "genomes/silva_tax/star_index_silva_{silva_ind}"
    output: temp("output/06-SILVA/{sample}/Aligned_{silva_ind}.txt.gz")
    wildcard_constraints:
        silva_ind="[0-9]+"
    params:
        STAR=config['softwares']['STAR'],
        samtools=config['softwares']['samtools'],
        prefix="06-SILVA/{sample}/Aligned_{silva_ind}/",
        useScratch = config['useScratch'],
        ram="64424509440"
    threads: 12
    script: '../src/star_align_silva.py'

rule silva_combine:
    input: 
        expand(
            "output/06-SILVA/{{sample}}/Aligned_{silva_ind}.txt.gz",
            silva_ind=range(1, config['silva_split'] + 1)
        )
    output: temp("output/06-SILVA/{sample}/Aligned_combine.txt.gz")
    threads: 3,
    shell: """
    zcat {input} | \\
    sed -E 's/:[LS]SU$//g' | \\
    gzip -c > {output}
    """

rule download_ncbi_taxdump:
    output: directory("genomes/taxdump")
    params:
        url = config["urls"]["ncbi_taxdump"]
    threads: 1
    shell: """
    wget {params.url} -O genomes/taxdump.tar.gz
    tar -zxf genomes/taxdump.tar.gz -c {output} 
    """

rule silva_count_taxa:
    input: 
        sample = "output/06-SILVA/{sample}/Aligned_combine.txt.gz",
        taxdump = "genomes/taxdump"
    output: 
        expand(
            "output/06-SILVA/{{sample}}/exogenousAligned_SILVA_{tax}.txt",
            tax=config["tax_levels"]
        )
    params: 
        hts_taxonomy = config["softwares"]['hts_taxonomy'],
        prefix = "output/06-SILVA/{sample}/exogenousAligned_SILVA_"
    threads: 4
    shell: """
    {params.hts_taxonomy} count-taxa \\
        --input-file {input.sample} \\
        --output-prefix {params.prefix} \\
        --nodes-dump  {input.taxdump}/nodes.dmp \\
        --names-dump {input.taxdump}/names.dmp
    """

rule silva_count_combine:
    input: 
        expand(
            "output/06-SILVA/{sample}/exogenousAligned_SILVA_{{tax}}.txt",
            sample=config["samples"]
        )
    output: "output/06-SILVA/{sample}/SILVA_count_{tax}.txt"
    threads: 4
    run: 
        from pandas import read_csv

        counts = None
        for sample in config["samples"].keys():
            path = 'output/06-SILVA/' + sample + '/exogenousAligned_SILVA_' + \
                   wildcards.tax + 'txt'
            new_counts = read_csv(
                path, sep='\t', index_col=0, names=[sample]
            )
            if counts is None:
                counts = new_counts
            else:
                counts = counts.join(new_counts, how='outer')
        counts.to_csv(output, sep='\t')
        
rule silva_extract_unmapped:
    input:
        aligned = "output/06-SILVA/{sample}/Aligned_combine.txt.gz",
        fastq = "output/05-RE/{sample}/Unmapped.fastq.gz"
    output: 
        unmapped = "output/06-SILVA/{sample}/Unmapped.fastq.gz",
    params:
        hts_fastx = config["softwares"]["hts_fastx"],
        namelist_path = "output/06-SILVA/{sample}/aligned_reads.txt",
        useScratch = config["useScratch"]
    threads: 3
    script: '../src/silva_extract_unmapped.py'