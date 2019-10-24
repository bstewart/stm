#### Dockerfile

```
### Get base image
FROM rocker/r-ubuntu:18.04
# As long as the above image doesn't get updated this should fix all of the package versions, which is important for reproducibility.


### Install additional packages (some of them are needed by devtools)
RUN apt-get update && \
    apt-get install -y \
    build-essential=12.4ubuntu1 \
    wget=1.19.4-1ubuntu2.2 \
    vim=2:8.0.1453-1ubuntu1.1 \
    libxml2-dev=2.9.4+dfsg1-6.1ubuntu1.2 \
    libssl-dev=1.1.1-1ubuntu2.1~18.04.4 \
    libcurl4-openssl-dev=7.58.0-2ubuntu3.8
# The above package versions are all fixed so these should not create a problem for stm reproducibility 


### Create directory      
RUN mkdir -p /install
RUN mkdir -p /run
RUN mkdir -p /input
RUN mkdir -p /output
# /install - holds the R install script
# /run - holds the default run script
# /input - can be used to mount a directory from the host system for input files
# /output - can be used to mount a directory from the host system for output files
# /input and /output allow for easy file transfers between the container and the host system

### Copy files
COPY docker_install.R /install/install.R
# Copy the desired install script named docker_install.R to the container's /install directory. The docker_install.R needs to reside in the PATH specified in the build command. E.g. if your build command is "docker build -t name/image ." then it needs to reside in the directory where docker build command is executed from.
COPY docker_run.R /run/run.R
# Copy the desired default run script named docker_run.R to the container's /run directory. The docker_run.R needs to reside in the PATH specified in the build command. E.g. if your build command is "docker build -t name/image ." then it needs to reside in the directory where docker build command is executed from.
RUN Rscript /install/install.R
# Run the install script to install desired R packages
CMD ["Rscript", "/run/run.R"]
# Specify a default command. This command will be executed if the "docker run" command has no arguments
```
#### R install script

```
install.packages(c("devtools", "stringr"))
library(devtools)
install_github("bstewart/stm",dependencies=TRUE)
``` 

#### Default R script

```
library(stringr)
library(stm)

mod <- stm(poliblog5k.docs, poliblog5k.voc, K=10, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100)

fname = str_c("../output/",format(Sys.time(), "%H-%M-%S"), ".mod")
save(mod, file = fname)
```

#### Building a custom image
- Using the above `Dockerfile` you can customize the image you build by modifying the `docker_install.R` script and the `docker_run.R` script. The former controls what R packages will be included in the container, the latter what the default run script should be.
- Note it is possible to install additional R packages in a live image but the changes may not be saved unless you do some extra steps.
- To build the container execute the command:

```
$ docker build -t <your-name>/<your-container-name> <PATH>
```

- Where you replace `<your-name>` with your [dockerhub](https://hub.docker.com/) username and `<your-container-name>` with a desired name for your container, and `<PATH>` with the path to the directory containing the `docker_install.R` and `docker_run.R` scripts.

- Further instructions on the `docker build` command can be found on [the docker website](https://docs.docker.com/engine/reference/commandline/build/).

#### Running an image
- Once you have an image there is many ways you can use it. The most important ways for this project would be the following:
    - **Default** run: `docker run --rm <your-name>/<your-container-name>`, which will run the default run, which is the last `[CMD]` statement in the `Dockerfile` in the above case it would be: `Rscript /run/run.R`
    - **With arguments**:  `docker run --rm <your-name>/<your-container-name> ARG1 ARG2 ...` For example, you could do: `Rscript /path/to/<my-other-run-script.R>` as `ARG1` and `ARG2`. `<my-other-run-script.R>` needs to be accessible from within the container.
    - **Interactive mode**:



    - **Mounting volumes mode**:
