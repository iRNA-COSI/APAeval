# Dockerfile to create docker container used to run DaPars execution workflow
FROM ubuntu:20.04

LABEL authors="Farica Zhuang" \
      description="Docker image containing all software requirements for the DaPars execution_workflow pipeline"

# set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y curl unzip python python-dev python3 python3-pip r-base bedtools wget

# Install DaPars
RUN wget https://github.com/ZhengXia/dapars/archive/546bdd4d809b097fd7a946d42b56e0efa0228472.zip
RUN mv *zip dapars.zip
RUN unzip -q "dapars.zip"
RUN mv dapars-* dapars
RUN chmod +x /dapars/src/*
ENV PATH=$PATH:/dapars/src/

# Download tools to convert gtf gene model file to bed
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod 755 gtfToGenePred genePredToBed
RUN mv gtfToGenePred /usr/local/bin
RUN mv genePredToBed /usr/local/bin

# Install python packages
RUN wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
RUN pip3 install pandas
RUN python2.7 get-pip.py
RUN pip2.7 install numpy
RUN pip2.7 install scipy
RUN pip2.7 install rpy2==2.8.6