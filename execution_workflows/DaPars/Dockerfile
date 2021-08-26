# Dockerfile to create docker container used to run DaPars execution workflow
FROM ubuntu:20.04

LABEL authors="Farica Zhuang" \
      description="Docker image containing all software requirements for the DaPars execution_workflow pipeline"

RUN apt-get update && apt-get install -y curl unzip python python-pip bedtools wget

# Install DaPars
RUN curl -LJO gitreleases.dev/gh/ZhengXia/dapars/latest/dapars_0.9.1.zip
RUN unzip -q "dapars_0.9.1.zip"
RUN chmod +x /dapars/src/*
ENV PATH=$PATH:/dapars/src/

# Download tools to convert gtf gene model file to bed
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod 755 gtfToGenePred genePredToBed
RUN mv gtfToGenePred /usr/local/bin
RUN mv genePredToBed /usr/local/bin

# Install python packages
RUN pip install pandas
RUN pip install configparser