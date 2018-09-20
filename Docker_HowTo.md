
# Deploying and Running QuickReduce via Docker

Steps: 
- configure the container
- build
- run/execute Quickreduce

## Configuring the container before building

There are two steps that determine the execution of QR from within the Docker
container:
- QR-specific parameters in the `podi_sitesetup.py` file
- Docker-specific configuration in the `docker-compose.yaml` file.

The most important aspects of the Docker configuration is the volume section in the docker-compose.yaml file. 
- */some/dir1:/input* mounts the local path */some/dir* into the Docker container;
- */some/dir2:/output* determines where the output files will be written to. 
Assuming the default sitesetup configuration this directory will also receive the debug 
and execution log files.
- */some/dir3:/cals* holds the main calibration files (bias, dark, and flat-fields)
- */some/dir4:/mastercal* are the typically WIYN-supplied long-term calibration products 
including fringe and pupilghost templates.
 
## Compiling/Building the Docker container

To build the Docker container and all its components, run

`docker-compose build qr`

from within the source directory. This will download the base image, install 
all required dependencies (python and python packages) and tools (SExtractor and Swarp).
Note that the building step might take a little while.


## Running QR from Docker

To execute the main QuickReduce task collectcells you can use the following command 
(maybe modified a bit to reflect input and output preferences).

`docker-compose run qr /app/podi_collectcells.py /input/o20171027T231924.2/ /output/docker_standalone.fits -fixwcs -photcalib -cals=/cals`
