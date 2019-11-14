shell.prefix("set +o pipefail; ")

configfile: "pipeline_config.yml"

# TODO: figure this out
report: "reports/workflow.rst"

rule all:
    input: 
        silva = expand(
            "output/06-SILVA/{sample}/exogenousAligned_SILVA_{tax}.txt",
            tax=config["tax_levels"], sample=config["samples"]
        ),
        unmapped = expand(
            "output/06-SILVA/{sample}/Unmapped.fastq.gz",
            sample=config["samples"]
        )

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
# include: "rules/bacteria_make_db.smk"