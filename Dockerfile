FROM continuumio/miniconda3
LABEL authors="Chenghao Zhu" \
      description="Docker image with all dependencies for excellular RNA sequencing data analysis pipeline (https://www.github.com/zhuchcn/exceRNApipeline.git)"

RUN apt-get update --fix-missing && \
    apt-get clean && \
    apt-get install unzip && \
    apt-get install -y gawk && \
    apt-get install -y gcc && \
    apt-get install -y g++ && \
    apt-get install -y swig && \
    apt-get install -y zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

COPY . /home/

RUN conda env create -f /home/environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/exceRNApipeline/bin:$PATH

RUN cd /home && \
    /opt/conda/envs/exceRNApipeline/bin/pip install . --install-option="--fully-install"
