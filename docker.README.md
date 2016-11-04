Development with Docker
===
This document contains instructions illustrating how to set up a [Docker](https://www.docker.com/) image that contains a compiled version of wnd-charm and other recommended python packages.

Table of Contents
===
* [Docker and WNDCHARM](#Docker-and-WNDCHARM)

Docker and WNDCHARM
===
In the [docker](./docker) folder, there is a `Dockerfile` that compiles and provides and interface to wnd-charm through
[Jupyterhub](https://github.com/jupyterhub/jupyterhub) and
[Jupyterlab](https://github.com/jupyterlab/jupyterlab).

Build the image
-----
From the parent folder of this repo:
```
docker build -t jupyterhub/wndcharm ./docker/.
```
Using the image
-----
Once the image is built, create a container:
```
docker run --name wndcharm -d -p 8000:8000 jupyterhub/wndcharm

docker exec -it wndcharm bash
passwd newuser
```
An interactive prompt will now open, choose a passwd for the newuser. After entering the password and confirming, you can exit the interactive session:
```
exit
```
wndcharm is now available for use at: http://localhost:8000. After logging in using `newuser` and the password entered above, choose `Notebook` and then choose the `Python 2` kernel when prompted (**not `Python 3`**). In the notebook you can run the following code to see the environment setup:
```
from wndcharm import diagnostics
print diagnostics
```
