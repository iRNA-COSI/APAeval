# Dockerfile to create docker container used to run APAlyzer execution workflow

##### BASE IMAGE #####
FROM r-base:latest

##### METADATA #####
LABEL base.image="r-base:latest"
LABEL version="1.0.6"
LABEL software="APAlyzer_R"
LABEL software.version="1.9.1"
LABEL software.description="APAlyzer bioconductor library"
LABEL software.website="https://github.com/RJWANGbioinfo/APAlyzer"
LABEL software.documentation="https://github.com/RJWANGbioinfo/APAlyzer#readme"
LABEL software.license="https://github.com/RJWANGbioinfo/APAlyzer/blob/master/LICENSE"
LABEL authors="Farica Zhuang and Maria Katsantoni" \
      description="Docker image containing all software requirements for APAlyzer execution workflow pipeline"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.9.1

# ##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y --allow-downgrades build-essential \
    curl \
    zlib1g-dev \
    software-properties-common \
    gcc \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    apt-transport-https

# Install r packages needed for APAlyzer
RUN Rscript -e 'install.packages(c("optparse", "BiocManager", "hash", "stringr"), dependencies=TRUE, repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("RJWANGbioinfo/APAlyzer"), dependencies=TRUE)'
