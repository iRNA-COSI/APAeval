FROM ubuntu:20.04

MAINTAINER Gregor Rot <gregor.rot@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

RUN apt update -y
RUN apt install software-properties-common -y
RUN apt-get install r-base python2.7 curl wget bedtools python-dev -y
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python2.7 get-pip.py
RUN pip2 install numpy
RUN pip2 install 'rpy2<2.9.0'
RUN mkdir /software
WORKDIR /software
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN chmod +x gtfToGenePred

CMD ["bash"]

#USER root
#RUN apt-get install python2.7
