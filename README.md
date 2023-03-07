
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

   + data: data sets used in the analysis protocol

   + scripts: code used in the analysis protocol

   + notebooks: Rmd notebooks used in the analysis protocol

   + results: results generated during analyses

   + `run_docker.sh`: bash script to launch `repro-totem-ti` docker image  

   + `run_slurm_singularity.sh`: bash script to submit a job to Slurm workload manager and launch the `repro-totem-ti` singularity image 
   + `rmd_to_rscript.R`: bash script to convert Rmd notebooks into R scripts

<br>

---

<br>

>Create docker image & launch container (locally) 

#### Build 'repro-totem-ti' Docker image: 
`docker build -t repro-totem-ti .`

#### List the Docker image - 'repro-totem-ti' (under REPOSITORY):
`docker images -a`

#### Create directory structure: 
`mkdir data scripts notebooks results`

#### Launch the container (see more instructions at: https://rocker-project.org/images/versioned/rstudio.html ):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/data:/home/rstudio/data -v $PWD/scripts:/home/rstudio/scripts \`

`	-v $PWD/notebooks:/home/rstudio/notebooks -v $PWD/results:/home/rstudio/results \`
	
`	repro-totem-ti`

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

>Create Singularity image from docker 'repro-totem-ti' Docker image to running it remotely in a server 

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

#### Open the file createdi 'slurm-<slurm-job.id>.out' and check the ssh command to type in a new local shell

#### Type the following in the browser and use your user name in the cluster as username and as password 'Totem':
`http://localhost:8787/`

