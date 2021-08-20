# Dockerfile to create docker container used to run APAtrap execution workflow
FROM rocker/r-ubuntu:20.04

# set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

LABEL authors="Farica Zhuang" \
      description="Docker image containing all software requirements for the APAtrap execution_workflow pipeline"

# Intstall curl, unzip, and python pip
RUN apt-get update && apt-get install -y curl unzip python3-pip r-base bedtools wget

# Install APAtrap
RUN curl -LJO https://sourceforge.net/projects/apatrap/files/APAtrap_Linux.zip/download
RUN unzip -q "APAtrap_Linux.zip"
ENV PATH=$PATH:/APAtrap/

# Unzip deAPA R package
RUN Rscript -e 'install.packages("stringr", dependencies=TRUE, INSTALL_opts = c("--no-lock"))'
RUN Rscript -e 'install.packages("APAtrap/deAPA_1.0.tar.gz", dependencies=TRUE, repos = NULL, type = "source")'

RUN pip install pandas

# Download tools to convert gtf gene model file to bed
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod 755 gtfToGenePred genePredToBed
RUN mv gtfToGenePred /usr/local/bin
RUN mv genePredToBed /usr/local/bin