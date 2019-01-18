[![Build Status](https://travis-ci.org/mklarqvist/tomahawk.svg?branch=master)](https://travis-ci.org/mklarqvist/tomahawk)
[![Release](https://img.shields.io/badge/Release-beta_0.7.0-blue.svg?logo=docker)](https://github.com/mklarqvist/tomahawk/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

# Getting started: Docker
## About
[Docker](https://docs.docker.com/install/) uses operating-system-level
virtualization into "containers". With Docker we can create virtual environments
that isolate a Tomahawk installation from the rest of the system. Tomahawk
programs are run within this virtual environment that can share resources with
its host machine. The Tomahawk [Docker
images](https://hub.docker.com/r/mklarqvist/tomahawk/) are tested for each release.

## Prerequisites

To follow this tutorial, you will need the following:

* One Ubuntu 18.04 server set up by following the Ubuntu 18.04 initial server setup guide, including a sudo non-root user and a firewall.
* An account on Docker Hub if you wish to create your own images and push them to Docker Hub, as shown in Steps 7 and 8.

## Getting started
You can pull the latest Tomahawk image using
```
docker pull mklarqvist/tomahawk
```

You can then run Tomahawk directly from that image. 
```bash
docker run -v /var/data/:/data latest mklarqvist/tomahawk:latest tomahawk calc -pi input.twk -o output.two
```
In this example we assume your data is available in the local folder
`/var/data/`. 
!!! note "Tutorial complete"
    
    Note: You now know how to use containerized Tomahawk.

You can run detach the docker instance so that it runs in the
background by passing the `-d` flag.
```bash
docker run -d -v /var/data/:/data latest mklarqvist/tomahawk:latest tomahawk calc -pi input.twk -o output.two
```
Let's demonstrate some more Tomahawk Docker recipes. You can start a `bash`
shell session within a Tomahawk container:
```bash
docker run -it -v /var/data/:/data latest mklarqvist/tomahawk:latest
```
Within the container you can interactively run the `tomahawk` CLI or develop
using the C++ API.

If you want to stop and cancel a running docker instance then run `docker ps` to
retrieve the target container id (first column) and execute `docker stop
CONTAINER_ID`.

!!! success "Tutorial complete"
    
    Success: You now know how to use containerized Tomahawk.