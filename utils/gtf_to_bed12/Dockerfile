# Use an official Python runtime as a parent image
FROM ubuntu:latest

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Set the working directory to /app
WORKDIR /

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3 python3-setuptools python3-dev python3-pip libcurl4
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod +x gtfToGenePred genePredToBed

# Copy the current directory contents into the container at /
COPY . /
