FROM continuumio/miniconda3
LABEL authors="Chenghao Zhu" \
      description="Docker image containing all requirements for zhuchcn/exceRNAseq pipeline"

RUN apt-get update --fix-missing && \
    apt-get clean && \
    apt-get install unzip && \
    apt-get install -y gawk && \
    rm -rf /var/lib/apt/lists/*

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/exceRNApipeline/bin:$PATH