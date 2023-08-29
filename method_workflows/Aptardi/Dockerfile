# Use an official Python runtime as a parent image
FROM nfcore/base:1.13

RUN conda create -n aptardi --channel conda-forge -c bioconda aptardi samtools bedtools pyranges pandas && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/aptardi/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name aptardi > aptardi.yml
