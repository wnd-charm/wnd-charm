Development with Docker
===
[Docker](https://www.docker.com/) to run, test, and debug our application.
The following documents how to install and use Docker on a Mac. There are
[instructions](https://docs.docker.com/installation) for installing and using
docker on other operating systems.

On the mac, we use [docker for mac](https://docs.docker.com/engine/installation/mac/#/docker-for-mac).

Table of Contents
===
* [Docker and WNDCHARM](#Docker-and-WNDCHARM)

Docker and WNDCHARM
===
In the [docker](./docker) file in this repo, we've included the necessary
files to compile and build a docker container that exposes wnd-charm using
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
