# Dockerfile to create docker container used to run IsoSCM execution workflow
FROM ubuntu:20.04

LABEL authors="Farica Zhuang and Chilam Poon" \
      description="Docker image containing all software requirements for IsoSCM execution_workflow pipeline"

# set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    bedtools \
    wget \
    samtools \
    openjdk-8-jdk \
    openjdk-8-jre \
    git \
    build-essential \
    unzip

RUN pip install pandas

# Setup JAVA_HOME & PATH -- in case the above commands not working..
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV PATH /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java:$PATH
RUN export JAVA_HOME
RUN export PATH

# Download isoscm
RUN wget https://github.com/shenkers/isoscm/releases/download/IsoSCM-2.0.12/IsoSCM-2.0.12.jar
RUN mv IsoSCM-2.0.12.jar IsoSCM.jar

# Download STAR
RUN git clone https://github.com/alexdobin/STAR.git
RUN cd STAR/source && make STAR
RUN cp STAR/source/STAR /usr/local/bin
