FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3 python3-setuptools python3-dev python3-pip

RUN apt-get install bedtools

COPY . /qapa
WORKDIR /qapa

RUN Rscript /qapa/scripts/install_packages.R

RUN python3 /qapa/setup.py install

ENTRYPOINT ["qapa"]
