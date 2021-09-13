# Use an official Python runtime as a parent image
FROM nfcore/base:1.13

# Install any needed packages specified in requirements.txt
#RUN apt-get update && apt-get install -y --no-install-recommends build-essential samtools
#RUN apt-get install bedtools
#RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN conda create -n aptardi --channel conda-forge -c bioconda aptardi samtools bedtools && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/aptardi/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name aptardi > aptardi.yml

# Copy the current directory contents into the container at /
#COPY . /
