rule download_silva:
    output:
        silva_lsu = temp("genomes/silva_lsu.fa.gz"),
        silva_ssu = temp("genomes/silva_ssu.fa.gz"),
        taxmap_lsu = temp("genomes/silva_lsu_taxmap.txt"),
        taxmap_ssu = temp("genomes/silva_ssu_taxmap.txt"),
    params:
        lsu = config['urls']['silva_lsu'],
        ssu = config['urls']['silva_ssu'],
        taxmap_lsu = config['urls']['silva_taxmap_lsu'],
        taxmap_ssu = config['urls']['silva_taxmap_ssu']
    shell: """
    wget {params.lsu} -O genomes/silva_lsu.fa.gz
    wget {params.ssu} -O genomes/silva_ssu.fa.gz
    wget {params.taxmap_lsu} -O genomes/silva_taxmap_lsu.txt.gz
    zcat genomes/silva_taxmap_lsu.txt.gz > {output.taxmap_lsu}
    wget {params.taxmap_ssu} -O genomes/silva_taxmap_ssu.txt.gz
    zcat genomes/silva_taxmap_ssu.txt.gz > {output.taxmap_ssu}
    """

rule filter_silva:
    input:
        silva_lsu = "genomes/silva_lsu.fa.gz",
        silva_ssu = "genomes/silva_ssu.fa.gz"
    output:
        silva_lsu = temp("genomes/silva_lsu_filter.fa.gz"),
        silva_ssu = temp("genomes/silva_ssu_filter.fa.gz"),
        silva_lsu_namelist = temp("genomes/silva_lsu_namelist.txt"),
        silva_ssu_namelist = temp("genomes/silva_ssu_namelist.txt")
    threads: 4
    shell: """
    zcat {input.silva_lsu} | grep '^>' |\\
        grep -e ' Bacteria;' -e ' Archaea;' \\
             -e ' Eukaryota;Opisthokonta;Nucletmycea;Fungi;' |\\
        sed 's/^>//g' > {output.silva_lsu_namelist}
    zcat {input.silva_ssu} | grep '^>' |\\
        grep -e ' Bacteria;' -e ' Archaea;' \\
             -e ' Eukaryota;Opisthokonta;Nucletmycea;Fungi;' |\\
        sed 's/^>//g' > {output.silva_ssu_namelist}
    hts_fastx extract-fastx \\
        --input-file {input.silva_lsu} \\
        --output-file {output.silva_lsu} \\
        --namelist-file {output.silva_lsu_namelist} \\
        --format fasta \\
        --verbose
    hts_fastx extract-fastx \\
        --input-file {input.silva_ssu} \\
        --output-file {output.silva_ssu} \\
        --namelist-file {output.silva_ssu_namelist} \\
        --format fasta \\
        --verbose
    """

rule solve_silva_taxa:
    input: 
        fasta = "genomes/silva_{su}_filter.fa.gz",
        taxmap = "genomes/silva_{su}_taxmap.txt"
    output: temp("genomes/silva_{su}_filter_ncbiTaxa.fa.gz")
    threads: 4
    script: "../src/silva_solve_taxa.py"

rule make_silva:
    input: expand("genomes/silva_{su}_filter_ncbiTaxa.fa.gz", su=["ssu", "lsu"])
    output:
        combine = temp("genomes/silva_filter_ncbiTaxa.fasta"),
        silva = temp(expand(
            "genomes/silva_tax/silva_{silva_ind}.fasta",
            silva_ind=range(1, config["parallel"]["silva"] + 1)
        ))
    params:
        parallel = config["parallel"]["silva"],
        output_prefix = "genomes/silva_tax/silva_"
    threads: 4
    shell: """
    zcat {input} | sed -E 's/ .+$/:SSU/' > {output.combine}
    hts_fastx split-fasta \\
        --input-file  {output.combine} \\
        --output-prefix {params.output_prefix} \\
        --n-batch {params.parallel}
    """


rule silva_index:
    input:
        genome = "genomes/silva_tax/silva_{silva_ind}.fasta"
    output:
        temp(directory("genomes/silva_tax/star_index_silva_{silva_ind}"))
    wildcard_constraints:
        silva_ind="[0-9]+"
    params:
        extra = "--genomeSAindexNbases 10 --genomeChrBinNbits 7",
        use_scratch = config['use_scratch']
    threads: 12
    script: '../src/star_index.py'