# Base image from Docker Hub:
FROM rocker/tidyverse:4.2.1

# Create directory structure:
RUN mkdir -p /home/rstudio/results /home/rstudio/scripts \
	/home/rstudio/notebooks /home/rstudio/data

# Copy data, Rmd notebooks & scripts: 
COPY data/human_cd34_bm_rep1.rds /home/rstudio/data
COPY notebooks/QuickStart_Totem.Rmd notebooks/GuidedStart_Totem.Rmd \
	notebooks/references.bib /home/rstudio/notebooks/
COPY scripts/QuickStart_Totem.R scripts/GuidedStart_Totem.R \
	/home/rstudio/scripts/

# Install software packages:
RUN apt update && apt upgrade -y \
	build-essential \
	libglpk40 \
	libcairo2-dev libxt-dev # 'scater' dependencies
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
RUN R --no-echo --no-restore --no-save -e "install.packages('hexbin')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('zellkonverter')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('scater')"

# Descriptions: 
MAINTAINER António Sousa (aggode@utu.fi)
LABEL Description: Container to reproduce single-cell RNA-seq trajectory inference analyses with Totem
LABEL Author(s): António Sousa (aggode@utu.fi), Johannes Smolander (johannes.smolander@utu.fi), Sini Junttila (sini.junttila@utu.fi), Laura L Elo (laura.elo@utu.fi)
LABEL Version 1.0.0

