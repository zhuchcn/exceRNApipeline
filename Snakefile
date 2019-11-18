shell.prefix("set +o pipefail; ")

configfile: "pipeline_config.yml"

# TODO: figure this out
report: "reports/workflow.rst"

if not config["exogenous_mapping"]:
    rule all:
        input:
            gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
            tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
            piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
            summary = "output/04-Genome/ReadsPerGene_summary.txt",
            bam = expand("output/04-Genome/{sample}/Aligned.out.bam", sample=config["samples"]),    
else:
    rule all:
        input: 
            gencode = "output/04-Genome/ReadsPerGene_gencode.txt",
            tRNA = "output/04-Genome/ReadsPerGene_tRNA.txt",
            piRNA = "output/04-Genome/ReadsPerGene_piRNA.txt",
            summary = "output/04-Genome/ReadsPerGene_summary.txt",
            bam = expand("output/04-Genome/{sample}/Aligned.out.bam", sample=config["samples"]),
            silva = expand("output/06-SILVA/SILVA_count_{tax}.txt",tax=config["tax_levels"]),
            bacteria = expand("output/07-Bacteria/bacteria_count_{tax}.txt",tax=config["tax_levels"])

# Step 1: preprocessing
include: "rules/preprocess.smk"

# Step 2: map to the UniVec library
include: "rules/univec.smk"

# Step 3: map to the rRNA sequences.
include: "rules/rrna.smk"

# Step 4: map to the human genome. 
include: "rules/endogenous_genome.smk"
# count endogenous genes
include: "rules/endogenous_count.smk"

# Step 5: map to repetitive elements
include: "rules/re.smk"

# Step 6: map to SILVA
include: "rules/silva_make_db.smk"
include: "rules/silva_mapping.smk"

# Step 7: map to Bacteria
include: "rules/bacteria_mapping.smk"