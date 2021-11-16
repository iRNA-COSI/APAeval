# Dockerfile to create docker container used to run APAlyzer execution workflow

##### BASE IMAGE #####
FROM r-base:4.1.0

##### METADATA #####
LABEL base.image="r-base:4.1.0"
LABEL version="1"
LABEL software="APAlyzer_R"
LABEL software.version="1.4.0"
LABEL software.description="APAlyzer bioconductor library"
LABEL software.website="https://github.com/RJWANGbioinfo/APAlyzer"
LABEL software.documentation="https://github.com/RJWANGbioinfo/APAlyzer#readme"
LABEL software.license="https://github.com/RJWANGbioinfo/APAlyzer/blob/master/LICENSE"
LABEL authors="Maria Katsantoni and Farica Zhuang" \
      description="Docker image containing all software requirements for the APAlyzer execution_workflow pipeline"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.4.0

# ##### INSTALL #####
# Install r packages needed for APAlyzer
RUN Rscript -e 'install.packages(c("optparse", "devtools", "RCurl", "BiocManager"), repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("RJWANGbioinfo/APAlyzer"), ask=FALSE, update = FALSE);'
RUN Rscript -e 'BiocManager::install(c("diffloop"), ask=FALSE, update = FALSE);'

# check if all of these are needed
RUN apt-get update -y \
  && apt-get install -y --allow-downgrades \
    build-essential \
    curl  \
    zlib1g-dev  \
    software-properties-common  \
    gcc  \
    libcurl4  \
    libcurl4-openssl-dev  \
    libxml2-dev  \
    libssl-dev  \
    apt-transport-https


