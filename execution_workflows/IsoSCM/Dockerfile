# Dockerfile to create docker container used to run IsoSCM execution workflow
FROM rocker/r-ubuntu:20.04

LABEL authors="Chilam Poon and Farica Zhuang" \
      description="Docker image containing all software requirements for the DaPars execution_workflow pipeline"

# set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

#RUN apt-get update && apt-get -y install --no-install-recommends --no-install-suggests \
#        ca-certificates software-properties-common gnupg2 gnupg1 \
#      && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

#
RUN apt-get update && apt-get install -y  \
    git \
    wget \
    sudo \
    openjdk-8-jdk \
    openjdk-8-jre \
    python3 \
    python3-dev \
    gcc

#RUN update-alternatives --config java
#RUN update-alternatives --config javac

# Setup JAVA_HOME & PATH -- in case the above commands not working..
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV PATH /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java:$PATH
RUN export JAVA_HOME
RUN export PATH

# Download isoscm
RUN wget https://github.com/shenkers/isoscm/releases/download/IsoSCM-2.0.12/IsoSCM-2.0.12.jar