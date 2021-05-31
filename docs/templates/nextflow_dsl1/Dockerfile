# USAGE: Replace all "[execution_workflow]" instances with the name of your workflow.
FROM nfcore/base:1.13.3
LABEL authors="your name" \
      description="Docker image containing all software requirements for the [execution_workflow] pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -qf /environment.yml && \
    sed -i 's/^conda activate base$/conda activate [execution_workflow]/' ~/.bashrc && \
    conda clean -qafy

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name [execution_workflow] > [execution_workflow].yml

# Instruct R processes to use these empty files instead of clashing with a local version
# please unmute the following line if you need R
#RUN touch .Rprofile .Renviron
