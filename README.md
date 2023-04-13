
<br>

---

<br>

Project folder: B23003_Totem_protocol

Description: repository with data analyses regarding the cell inference trajectory protocol with Totem to be published at Springer Nature's Methods in Molecular Biology

Author(s): António Sousa, Johannes Smolander, Sini Junttila, Laura Elo

Date: 16/02/2023

Archived date: 

<br>

---

<br>

### Content:

<br>

   + scripts: `download_h5ad_to_SCE_rds_script.R` R script used to obtain `human_cd34_bm_rep1.rds` object and `QuickStart_Totem.R` and `GuidedStart_Totem.R` R analysis scripts

   + notebooks: `QuickStart_Totem.Rmd` and `GuidedStart_Totem.Rmd` R markdown notebooks (and respective compliled `.html` files) used in the analyses

   + `run_docker.sh`: bash script to launch `repro-totem-ti` docker image  

   + `run_slurm_singularity.sh`: bash script to submit a job to Slurm workload manager and launch the `repro-totem-ti` singularity image 

   + `rmd_to_rscript.R`: bash script to convert Rmd notebooks into R scripts

<br>

---

<br>

### Data 

<br>

The data set `human_cd34_bm_rep1.rds` was generated with the R script `download_h5ad_to_SCE_rds_script.R`. It is a parsed `SingleCellExperiment` `RDS` object corresponding to the anndata h5ad `human_cd34_bm_rep1.h5ad` available on [HCA Portal]() and published by [Setty et al., 2019](https://www.nature.com/articles/s41587-019-0068-4).

<br>

---

<br>

>Build `repro-totem-ti` locally:

### Clone through SSH password-protected key this repository: 
`git clone git@github.com:elolab/Totem-protocol.git`

#### Build `repro-totem-ti` Docker image from `Dockerfile`:
`docker build -t repro-totem-ti .` 

#### List the Docker image `repro-totem-ti` (under REPOSITORY):
`docker images -a`

#### Create folder to save results: 
`mkdir -p repro-totem-ti/results`

`cd repro-totem-ti`

#### Launch the container (see more instructions at: https://rocker-project.org/images/versioned/rstudio.html ):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/results:/home/rstudio/results \`
	
`	repro-totem-ti`

#### Type the following hyperlink in the browser: 
`http://localhost:8787/`

#### Use the following credentials: 
`Username: rstudio`

`Password: Totem`

<br>

---

<br>

>Launch container (locally) from Docker Hub `elolab/repro-totem-ti` 

#### Create folder to save results: 
`mkdir -p repro-totem-ti/results`

`cd repro-totem-ti`

#### Launch the container (see more instructions at: https://rocker-project.org/images/versioned/rstudio.html ):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/results:/home/rstudio/results \`
	
`	elolab/repro-totem-ti`

#### Run the bash script `run_docker.sh` as an alternative to the command above (optional):
`./run_docker.sh`

#### Type the following hyperlink in the browser: 
`http://localhost:8787/`

#### Use the following credentials: 
`Username: rstudio`

`Password: Totem`

<br>

---

<br>

>Create Singularity image from docker 'repro-totem-ti' Docker image locally to running it remotely in a server 

#### Check the IMAGE ID of 'repro-totem-ti': 
`docker images`

#### Create 'imgs' directory and save the tarball 'repro-totem-ti.tar':
`mkdir imgs`

`sudo docker save <IMAGE ID> -o imgs/repro-totem-ti.tar`

#### Create Singularity image from tarball 'repro-totem-ti.tar':
`cd imgs`

`sudo singularity build repro-totem-ti.sif docker-archive://repro-totem-ti.tar`

#### Upload this repository to the server

#### Execute the 'run_slurm_singularity.sh' bash script with Slurm to launch the Singularity container:
`sbatch ./run_slurm_singularity.sh Totem`

#### Open the file created 'slurm-<slurm-job.id>.out' and check the ssh command to type in a new local shell

#### Type the following in the browser and use your user name in the cluster as username and as password 'Totem':
`http://localhost:8787/`

<br>

---

<br>

>Create cross-platform (amd64,arm64) `repro-totem-ti` image for compatibility with arm64 M1 MacBooks (ongoing)

#### Build and push the image to temporary Docker Hub repository (first install `buildx`, `binfmt` and set up `mybuilder`):
`docker buildx build --no-cache --push --platform linux/amd64,linux/arm64 -t aggode/repro-totem-ti:cross-platform -f Dockerfile-cross-platform .`

