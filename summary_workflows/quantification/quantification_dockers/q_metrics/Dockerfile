# Use an official Python runtime as a parent image
FROM ubuntu:20.04

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Set the working directory to /app
WORKDIR /app

# Copy the requirements.txt into the container at /app
COPY [ "requirements.txt","constraints.txt", "/app/" ]

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3 python3-setuptools python3-dev python3-pip git
RUN apt-get install bedtools
RUN apt-get update && apt-get -y install procps && rm -rf /var/lib/apt/lists/* && pip install --no-cache-dir  --trusted-host pypi.python.org -r requirements.txt -c constraints.txt

# Copy the current directory contents into the container at /app
COPY . /app
