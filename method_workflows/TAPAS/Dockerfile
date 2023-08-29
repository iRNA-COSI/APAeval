# Please run this Dockerfile in a cloned TAPAS git repository: https://github.com/arefeen/TAPAS  
FROM rocker/r-base:latest

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Set the working directory to /app
WORKDIR /

## COPY
COPY ./package_install.R /package_install.R
COPY ./Finding_APA_Sites/APA_sites_detection /APA_sites_detection
COPY ./Differential_APA_Site_Analysis/Diff_APA_site_analysis /Diff_APA_site_analysis

RUN apt-get update && apt-get install -y procps
RUN Rscript /package_install.R
