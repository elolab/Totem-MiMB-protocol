#!/bin/bash
docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \
	-v $PWD/data:/home/rstudio/data -v $PWD/scripts:/home/rstudio/scripts \
	-v $PWD/notebooks:/home/rstudio/notebooks -v $PWD/results:/home/rstudio/results \
	repro-totem-ti 

