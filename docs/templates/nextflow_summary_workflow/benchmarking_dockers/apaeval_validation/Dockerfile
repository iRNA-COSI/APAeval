# Use an official Python runtime as a parent image
FROM python:2.7-slim

# Set the working directory to /app
WORKDIR /app

# Copy the requirements.txt into the container at /app
COPY [ "requirements.txt","constraints.txt", "/app/" ]

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get -y install procps && rm -rf /var/lib/apt/lists/* && pip install --no-cache-dir --trusted-host pypi.python.org -r requirements.txt -c constraints.txt

# Copy the current directory contents into the container at /app
COPY . /app

ENTRYPOINT ["python","validation.py"]
