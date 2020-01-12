shell.prefix("set +o pipefail; ")

singularity: "docker://zhuchcn/exce-rna-pipeline:latest"

configfile: "default_config.yml"
configfile: "pipeline_config.yml"

report: "reports/workflow.rst"

def all_input():
    inputs = [
        "output/results/endogenous/summary.png", 
        "output/results/qc/preprocess.tsv"
    ]
    inputs.extend(expand(
        "output/results/qc/{sample}_SE_fastqc/Images/per_base_quality.png",
        sample=config["samples"]
    ))
    if config["exogenous_mapping"]:
        inputs.extend(["output/results/silva/", "output/results/bacteria/"])
    return inputs

rule all:
    input: all_input()

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