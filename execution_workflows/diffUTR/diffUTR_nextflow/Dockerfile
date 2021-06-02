# USAGE: Replace all "[execution_workflow]" instances with the name of your workflow.
FROM nfcore/base:1.13.3
LABEL authors="Chi-Lam Poon" \
      description="Docker image containing all software requirements for the diffUTR pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -qf /environment.yml && \
    sed -i 's/^conda activate base$/conda activate diffUTR_nextflow/' ~/.bashrc && \
    conda clean -qafy

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name diffUTR_nextflow > diffUTR_nextflow.yml

# Instruct R processes to use these empty files instead of clashing with a local version
# please unmute the following line if you need R
RUN touch .Rprofile .Renviron
