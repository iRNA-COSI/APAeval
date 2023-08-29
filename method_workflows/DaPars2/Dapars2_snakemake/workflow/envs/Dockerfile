# Template Dockerfile with python and pandas.
##### BASE IMAGE #####
FROM python:2.7.9

##### METADATA #####
LABEL base.image="python:2.7.9"
LABEL version="1"
LABEL software="python"
LABEL software.version="2.7.9"
LABEL software.description="Python programming language"
LABEL software.website="https://www.python.org/"
LABEL software.documentation="https://docs.python.org/2/"
LABEL software.license="https://docs.python.org/2/license.html"
LABEL software.tags="General"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.7.9
ENV DEBIAN_FRONTEND noninteractive

##### INSTALL #####

# Note - gfortran libopenblas-dev liblapack-dev are needed for numpy 1.11.3, get pip error without pre-installing
# https://github.com/scipy/scipy/issues/9005

# Install Python dependencies

RUN apt-get update -y \
  && apt-get install -y git gfortran libopenblas-dev liblapack-dev \
  #&& pip install --upgrade pip \
  && pip install numpy==1.11.3 && pip install scipy==0.17.1


# Install Dapars2
RUN git clone https://github.com/3UTR/DaPars2 && cd DaPars2 && git checkout 23d89d1

# Install bedtools from Github (specific version)

WORKDIR /tmp/bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz \
  && tar -xvzf bedtools-2.30.0.tar.gz \
  && cd bedtools2 \
  && make \
  && make install

# Install samtools (v 1.12)
WORKDIR /tmp/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 \
  && tar -xf samtools-1.12.tar.bz2 \
  && cd samtools-1.12 \
  && ./configure \
  && make \
  && make install


# install gff3ToGenePred & genePredToBed
WORKDIR /tmp/kent
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred  http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed \
  && chmod a+x * \
  && mv * /usr/local/bin


# Clean up
RUN apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /DaPars2
