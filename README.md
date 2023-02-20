
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

<br>

---

<br>

### Create docker image & launch container (locally): 

<br>

#### Build 'scrna-cell-traject' Docker image: 
`docker build -t scrna-cell-traject .`

#### List the Docker image - 'scrna-cell-traject' (under REPOSITORY):
`docker images -a`

#### Create directory structure: 
`mkdir data scripts notebooks results`

#### Launch the container (see more instructions at: https://rocker-project.org/images/versioned/rstudio.html ):
`docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \`

`	-v $PWD/data:/home/rstudio/data -v $PWD/scripts:/home/rstudio/scripts \`

`	-v $PWD/notebooks:/home/rstudio/notebooks -v $PWD/results:/home/rstudio/results \`
	
`	scrna-cell-traject`

#### Type the following hyperlink in the browser: 
`http://localhost:8787/`

#### Use the following credentials: 
`Username: rstudio`

`Password: Totem`

<br>

---

<br>

## Create Singularity image from docker 'scrna-cell-traject' Docker image to running it remotely in a server: 

<br>

#### Check the IMAGE ID of 'scrna-cell-traject': 

`docker images`

#### Create 'imgs' directory and save the tarball 'scrna-cell-ti.tar'

`mkdir imgs`

`sudo docker save <IMAGE ID> -o imgs/scrna-cell-ti.tar`

#### Create Singularity image from tarball 'scrna-cell-ti.tar'

`cd imgs`

`sudo singularity build scrna-cell-ti.sif docker-archive://scrna-cell-ti.tar`

#### Upload this repository to the server 

#### Execute the 'run_slurm_singularity.sh' bash script with Slurm to launch the Singularity container

`sbatch ./run_slurm_singularity.sh Totem`

#### Open the file createdi 'slurm-<slurm-job.id>.out' and check the ssh command to type in a new local shell

#### Type the following in the browser and use your user name in the cluster as username and as password 'Totem'

`http://localhost:8787/`

