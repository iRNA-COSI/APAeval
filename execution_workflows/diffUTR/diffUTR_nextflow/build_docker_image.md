## Build Docker image from this directory with enviroment.yml and Dockerfile
```
docker build . -t [execution_tool]:1.0
```
If you want to update the docker container, please remove your original image first:
```
docker image ls #look for the IMAGE_ID of your docker image
docker rmi [IMAGE_ID]
```
