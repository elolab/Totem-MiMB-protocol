#!/bin/bash
docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \
	-v $PWD/results:/home/rstudio/results \
	elolab/repro-totem-ti 

