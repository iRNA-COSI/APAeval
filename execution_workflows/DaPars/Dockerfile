# Dockerfile to create docker container used to run APAtrap execution workflow
FROM ubuntu:20.04

# set this to not get asked for geographic area when downloading R
ENV DEBIAN_FRONTEND noninteractive

LABEL authors="Farica Zhuang" \
      description="Docker image containing all software requirements for the APAtrap execution_workflow pipeline"

RUN apt-get update && apt-get install -y curl unzip python bedtools

# Install DaParsls
RUN curl -LJO gitreleases.dev/gh/ZhengXia/dapars/latest/dapars_0.9.1.zip
RUN unzip -q "dapars_0.9.1.zip"
#RUN mv dapars/src/DaPars_Extract_Anno.py /usr/local/bin/.
#RUN mv dapars/src/DaPars_main.py .
RUN chmod +x /dapars/src/*
ENV PATH=$PATH:/dapars/src/

# Unzip deAPA R package
#RUN Rscript -e 'install.packages("stringr", dependencies=TRUE, INSTALL_opts = c("--no-lock"))'
#RUN Rscript -e 'install.packages("APAtrap/deAPA_1.0.tar.gz", dependencies=TRUE, repos = NULL, type = "source")'
#
#RUN pip install scipy
#RUN pip install numpy