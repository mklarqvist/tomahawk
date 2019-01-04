Dockerized Tomahawk
================

[Tomahawk](https://github.com/mklarqvist/tomahawk) is available as a [Docker
image](https://hub.docker.com/r/mklarqvist/tomahawk/). You can pull the Tomahawk
image using

`docker pull mklarqvist/tomahawk`

You can then run Tomahawk directly from that image. 
```bash
docker run -v /var/data/:/data latest mklarqvist/tomahawk:latest tomahawk calc -pi input.twk -o output.two
```
In this example we assume your data is available in the local folder 
`/var/data/`
You can run detach the docker instance so that it runs in the background by passing the `-d` flag.
```bash
docker run -d -v /var/data/:/data latest mklarqvist/tomahawk:latest tomahawk calc -pi input.twk -o output.two
```
Lastly, you can also run the docker instance interactively. By default, the container will run an instance of `/bin/bash`.
```bash
docker run -it -v /var/data/:/data latest mklarqvist/tomahawk:latest
```

If you want to stop and cancel a running docker instance then run `docker ps` to
retrieve the target container id (first column) and execute `docker stop
CONTAINER_ID`.