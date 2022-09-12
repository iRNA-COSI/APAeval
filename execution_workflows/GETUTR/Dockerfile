FROM ubuntu:20.04

MAINTAINER Gregor Rot <gregor.rot@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

RUN apt update -y
RUN apt install software-properties-common -y
RUN apt-get update -y
RUN apt-get install r-base python2.7 curl wget bedtools python-dev samtools -y
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2.7 get-pip.py
RUN pip2 install numpy
RUN pip2 install 'rpy2<2.9.0'
RUN mkdir /software
WORKDIR /software
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN chmod +x gtfToGenePred
ENV PATH="${PATH}:/software"

RUN curl http://big.hanyang.ac.kr/GETUTR/samples/GETUTR.2.0.0.tar.gz --output GETUTR.2.0.0.tar.gz
RUN tar xfz GETUTR.2.0.0.tar.gz
RUN mv GETUTR.2.0.0_for_rpy2 GETUTR