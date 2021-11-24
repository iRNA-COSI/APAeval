# Template Dockerfile with python and pandas.
##### BASE IMAGE #####
FROM python:3.9.4

##### METADATA #####
LABEL base.image="python:3.9.4"
LABEL version="1"
LABEL software="python"
LABEL software.version="3.9.4"
LABEL software.description="Python programming language"
LABEL software.website="https://www.python.org/"
LABEL software.documentation="https://docs.python.org/3/"
LABEL software.license="https://docs.python.org/3/license.html"
LABEL software.tags="General"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.9.4

##### INSTALL #####
RUN apt-get update -y \
  && pip install --upgrade pip \
  && pip install pandas==1.2.4 \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \