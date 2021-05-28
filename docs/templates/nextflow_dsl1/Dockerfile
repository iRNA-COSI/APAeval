FROM nfcore/base:1.13.3
LABEL authors="your name" \
      description="Docker image containing all software requirements for the testing pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -qf /environment.yml && \
    sed -i 's/^conda activate base$/conda activate testing/' ~/.bashrc && \
    conda clean -qafy

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/testing/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name testing > testing.yml

# Instruct R processes to use these empty files instead of clashing with a local version
# please unmute the following line if you need R
#RUN touch .Rprofile .Renviron
