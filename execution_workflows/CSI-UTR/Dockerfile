FROM nfcore/base:1.13.3
LABEL authors="Melina Klostermann" \
      description="Docker image containing all software requirements for the testing pipeline of CSI-UTR"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f environment.yml && \
    sed -i 's/^conda activate base$/conda activate nf-core-CSI-UTR-1.1.0/' ~/.bashrc && \
    conda clean -qafy

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH "/opt/conda/envs/nf-core-CSI-UTR-1.1.0/bin:$PATH"

# add some installation tools
RUN  apt-get update -y && \
     apt-get -y autoremove && \
     apt-get clean

RUN apt-get install -y build-essential


# Add tools not in conda channals: perl packages
RUN apt-get install -y perl && \
apt-get install -y cpanminus
RUN cpanm Text::NSP::Measures::2D::Fisher::twotailed && \
 cpanm Statistics::Multtest && \
 cpanm Statistics::TTest && \
 cpanm List::MoreUtils && \
 cpanm Getopt::Long && \
 cpanm MIME::Base64 && \
 cpanm File::Which



# Add csi-utr
RUN git clone https://github.com/UofLBioinformatics/CSI-UTR.git
RUN mv CSI-UTR/CSI-UTR_v1.1.0/bin/CSI-UTR usr/bin/
## && mv CSI-UTR/CSI-UTR_v1.1.0/data usr/




# Dump the details of the installed packages to a file for posterity
RUN conda env export --name testing > testing.yml

# Instruct R processes to use these empty files instead of clashing with a local version
# please unmute the following line if you need R
RUN apt-get install -y r-base-core
RUN touch .Rprofile .Renviron
