# Human Extracellular small RNA-seq pipeline

Data processing pipeline for extracellular small RNA-seq from human specimen. This pipeline is designated for running on HPC with the job management system [SLURM](https://slurm.schedmd.com/sbatch.html). The pipeline has four steps:

1. Preprocess: remove adapters and trimming low quality nucleotides, using [HTStream](https://github.com/ibest/HTStream).
2. UniVec: map to the NCBI's [UniVec](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) database to remove vector origins.
3. RiboRNA: map to the rRNA sequences.
4. Human Genome: map to human genome.

The aligner STAR is used to align reads to database/genome. The pipeline [exceRpt](https://github.com/gersteinlab/exceRpt) was used as reference. The advantage of this pipeline is that it can take the advantage of HPC's full power by submitting jobs in parallel using SLURM.

## Requirement

```
python==3.7
snakemake==5.7.4
pandas==0.25.2
```

## Installation

Install directly from github

```
git clone 
```

## pipeline configuration

### samples

The pipeline configuration is set up in the `pipeline_config.yml` file. The samples need to be given in key value paires in the samples section. If relative path is given, it must refer to the directory of the `Snakefile`. 

### scratch

Depends on your HPC setting, some allows users to use a `/scratch` folder on each node to avoid too much IO. If your HPC does not use `/scratch`, turn the `useScratch` off by setting it to `false`.

### softwares

Softwares must be given as absolute path in the softwares section. Your HPC should have all the softwares installed already. For example, in order to find the path ot the `STAR`, you first find all available STAR by:

```
module avail | grep star
# star/2.4.2a
# star/2.5.2b
# star/2.6.1a
# star/2.7.2a
```

Choose one and load the module by:

```
module load star/2.7.2a
# Module star/2.7.2a loaded
```

Find the path to `STAR`

```
which STAR
# /share/apps/star-2.7.2a/bin/STAR
```

### genomes and annotations

All genome and annotation files can be automatically downloaded by the pipeline, except the complete human ribosomal RNA sequence, U13369. If must be manually download from NCBI's website and put into the `genomes` folder. The path is blow and also can be found in the `pipeline_config.yml`:

https://ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta

## slurm configuration

Set up the slurm options by modifying the `slurm_config.yml` file. The default slurm options are set in the `__default__` seciont, but it can be overwritten for each rule in its own section. The default of the `mail-user` is reading from the environmental variable `$USER_EMAIL`. You can either change it to your own email address, or setting it in bash like:

```
export USER_EMAIL=you@email.com
```

## Run pipeline

When all configurations are set up, run the pipeline using the following command.

```
./snakemakeslurm
```

## Results

The results are in the `output` folder. The per gene read counts for gencode mapped genes, tRNAs, and piRNAs, as well as a gene type summary table, are in a subfolder named `04-Genome`.

## clean scratches

If scratch is used, the files in scratch should be removed by the pipeline automatically. However in case that jobs are stopped or cancelled in the middle, the scratch files will stays no the node in the `/scratch` folder. In those above mentioned cases, use the command below to clean the scratches.

```
./snakemakeslurm --clean-scratch
```

A `slurm_job_status` folder will be created after jobs are submitted. So basically, if the `slurm_job_status` is not empty, run the commdan above.

## TODO

- [ ] Add workflow report.
- [ ] Add conda environment support?
- [ ] Add exogenous genomes.