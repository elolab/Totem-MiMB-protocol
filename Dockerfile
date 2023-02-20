FROM rocker/tidyverse:4.2.1

# Create directory structure:
RUN mkdir -p /home/rstudio/results /home/rstudio/scripts /home/rstudio/notebooks /home/rstudio/data

# Install software packages:
RUN apt update && apt upgrade -y \
	build-essential \
	libglpk40
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('slingshot')"
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('elolab/Totem', ref='fc2d20614a1d82459ca38bb98f22bd212922aec3')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('DelayedMatrixStats')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('Matrix')"
RUN R --no-echo --no-restore --no-save -e "install.packages('markdown')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('BiocStyle')" 
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('scran')" 
RUN R --no-echo --no-restore --no-save -e "install.packages('uwot')"
RUN R --no-echo --no-restore --no-save -e "install.packages('optparse')"

# Descriptions: 
MAINTAINER António Sousa (aggode@utu.fi)
LABEL Description: Container for single-cell RNA-seq trajectory inference analyses with Totem
LABEL Author(s): António Sousa (aggode@utu.fi), Johannes Smolander (johannes.smolander@utu.fi), Sini Junttila (sini.junttila@utu.fi), Laura L Elo (laura.elo@utu.fi)
LABEL Version 1.0.0

