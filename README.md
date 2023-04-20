# Totem protocols for single-cell trajectory inference

Repository with data analyses regarding the cell inference trajectory protocol with _Totem_ to be published at Springer Nature's Methods in Molecular Biology: _User-friendly protocols for inferring tree-shaped single-cell trajectories with Totem_

Author(s): António G.G. Sousa (<aggode@utu.fi>), Johannes Smolander (<johannes.smolander@utu.fi>), Sini Junttila (<simaju@utu.fi>), Laura L. Elo (<laura.elo@utu.fi>) 

<br>

## Contents:

* [Data Set](#data-set)

* [Directory Structure](#directory-structure)

* [Build Docker Image](#build-docker-image)

* [Launch Container Locally](#launch-container-locally)

* [Launch Container Remotely](#launch-container-remotely)

<br>

## Data Set

The data set `human_cd34_bm_rep1.rds` was parsed with the R script `download_h5ad_to_SCE_rds_script.R` (under the `scripts` folder). It is a parsed `SingleCellExperiment` `RDS` object corresponding to the anndata h5ad `human_cd34_bm_rep1.h5ad` available on [HCA Portal]() and published by [Setty et al., 2019](https://www.nature.com/articles/s41587-019-0068-4).

This parsed data set is distributed with the docker image ([elolab/repro-totem-ti](https://hub.docker.com/r/elolab/repro-totem-ti) - see section [Launch Container Locally](#launch-container-locally)). Alternatively, it can be downloaded from Zenodo: [10.5281/zenodo.7845709](https://doi.org/10.5281/zenodo.7845709).

<br>

## Directory Structure

   + `docs`: folder comprising the material used by GitHub Pages to display the website documentation: [https://elolab.github.io/Totem-protocol](https://elolab.github.io/Totem-protocol)

   + `sphinx-website`: folder comprising files used to build the website documentation with the [Sphinx](https://www.sphinx-doc.org/en/master/) application that are used in the `docs` folder to display the website
   
   + `Dockerfile`: file with the instructions to build docker image

   + `scripts` folder: 
   
      + `download_h5ad_to_SCE_rds_script.R`: R script used to obtain `human_cd34_bm_rep1.rds` object
      
      + `QuickStart_Totem.R`: **QuickStart** R script analysis protocol
      
      + `GuidedStart_Totem.R`: **GuidedStart** R script analysis protocol

   + `notebooks` folder: 
   
      + `QuickStart_Totem.Rmd`: **QuickStart** R markdown notebook analysis protocol

      + `QuickStart_Totem.html`: **QuickStart** R markdown analysis protocol compiled into `html` format
      
      + `GuidedStart_Totem.Rmd`: **GuidedStart** R markdown notebook analysis protocol 

      + `GuidedStart_Totem.html`: **GuidedStart** R markdown analysis protocol compiled into `html` format

   + `run_docker.sh`: bash script to launch ([elolab/repro-totem-ti](https://hub.docker.com/r/elolabs41587-019-0068-4) docker container  

   + `run_slurm_singularity.sh`: bash script to submit a job to Slurm workload manager and launch the ([elolab/repro-totem-ti](https://hub.docker.com/r/elolabs41587-019-0068-4) container through singularity in a cluster environment

   + `rmd_to_rscript.R`: bash script to convert Rmd notebooks into R scripts

<br>

## Build Docker Image

### Clone through SSH password-protected key this repository: 
`git clone git@github.com:elolab/Totem-protocol.git`

### Enter the repository folder cloned: 
`cd Totem-protocol`

### Create `data` folder: 
`mkdir data results`

### Download data from Zenodo repository: 
`wget https://zenodo.org/record/7845709/files/human_cd34_bm_rep1.rds?download=1 \`

`    -O data/human_cd34_bm_rep1.rds

### Build `repro-totem-ti` Docker image from `Dockerfile`:
`docker build -t repro-totem-ti .` 

### List the docker image `repro-totem-ti` (under REPOSITORY):
`docker images -a`

### Launch the container (see more instructions at: https://rocker-project.org/images/versioned/rstudio.html ):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/results:/home/rstudio/results \`
	
`	repro-totem-ti`

### Type the following hyperlink in the browser: 
`http://localhost:8787/`

### Use the following credentials: 
`Username: rstudio`

`Password: Totem`

<br>

## Launch Container Locally

### Create folder to save results: 
`mkdir -p repro-totem-ti/results`

`cd repro-totem-ti`

### Launch the container (see more instructions at [rocker](https://rocker-project.org/images/versioned/rstudio.html)):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/results:/home/rstudio/results \`
	
`	elolab/repro-totem-ti`

### Run the bash script `run_docker.sh` as an alternative to the command above (optional):
`./run_docker.sh`

### Type the following hyperlink in the browser: 
`http://localhost:8787/`

### Use the following credentials: 
`Username: rstudio`

`Password: Totem`

<br>

## Launch Container Remotely

A docker image can be converted into a Singularity image locally (1) or just pushed from Docker Hub (2). Below it is presented both solutions. 

<br>

>(1) Convert `repro-totem-ti` docker image created locally into a Singularity image

### Check the IMAGE ID of `repro-totem-ti` created above: 
`docker images`

### Create 'imgs' directory and save the tarball `repro-totem-ti.tar`:
`mkdir imgs`

`sudo docker save <IMAGE ID> -o imgs/repro-totem-ti.tar`

### Create Singularity image from tarball `repro-totem-ti.tar`:
`cd imgs`

`sudo singularity build repro-totem-ti.sif docker-archive://repro-totem-ti.tar`

### Upload this repository to the server environment  

<br>

>(2) Pull image from Docker Hub `elolab/repro-totem-ti`

### Clone through SSH password-protected key this repository to your server/cluster environment: 
`git clone git@github.com:elolab/Totem-protocol.git`

### Enter the repository folder cloned: 
`cd Totem-protocol`

### Create directory structure: 
`mkdir results imgs` 

### Pull the `elolab/repro-totem-ti` image from Docker Hub into the `imgs` folder: 
`singularity pull --dir imgs/ repro-totem-ti.sif docker://elolab/repro-totem-ti:latest`

<br>

>Submit Slurm job to launch RStudio server through the Singularity image and access it locally

### Execute the 'run_slurm_singularity.sh' bash script with Slurm to launch the Singularity container:
`sbatch ./run_slurm_singularity.sh Totem`

### Open the file created 'slurm-<slurm-job.id>.out' and check the ssh command to type in a new local shell

### Type the following in the browser and use your user name in the cluster as user name and as password `Totem`:
`http://localhost:8787/`

