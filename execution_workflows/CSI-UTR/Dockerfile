# Dockerfile to create docker container used to run DaPars execution workflow
FROM ubuntu:20.04

LABEL authors="Farica Zhuang" \
      description="Docker image containing all software requirements for the CSI-UTR execution workflow pipeline"

# Set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

# Install required R packages
RUN apt-get update
RUN apt-get install -y libxml2-dev \
			r-cran-xml \
			libcurl4-openssl-dev \
			libssl-dev \
			wget \
			lsb-release \
			software-properties-common
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update
RUN apt install -y r-base=4.1.3-1.2004.0
RUN Rscript -e 'install.packages("BiocManager",  repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("Rcpp"),  ask=FALSE, update = FALSE, dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("locfit"), ask=FALSE, update = FALSE, dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("biomaRt"), ask=FALSE, update = FALSE, dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("GenomicFeatures"), ask=FALSE, update = FALSE, dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("DESeq2"), ask=FALSE, update = FALSE, dependencies=TRUE)'
RUN Rscript -e 'BiocManager::install(c("DEXSeq"), ask=FALSE, update = FALSE, dependencies=TRUE)'

# Install required tools and modules for CSI-UTR
RUN apt-get update && apt-get install -y curl unzip perl git bedtools samtools liblist-moreutils-perl python3
RUN curl -L http://cpanmin.us | perl - Statistics::TTest
RUN curl -L http://cpanmin.us | perl - Statistics::Multtest
RUN curl -L http://cpanmin.us | perl - File::Which
RUN curl -L http://cpanmin.us | perl - Text::NSP::Measures::2D::Fisher::twotailed

# Install CSI-UTR
RUN git clone https://github.com/UofLBioinformatics/CSI-UTR.git
RUN chmod +x /CSI-UTR/CSI-UTR_v1.1.0/bin/CSI-UTR
ENV PATH=$PATH:/CSI-UTR/CSI-UTR_v1.1.0/bin/

# Download tools to convert gtf gene model file to bed
RUN apt-get install wget
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod 755 gtfToGenePred genePredToBed
RUN mv gtfToGenePred /usr/local/bin
RUN mv genePredToBed /usr/local/bin


