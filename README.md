
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

### Create docker image & launch container: 

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
	`-v $PWD/notebooks:/home/rstudio/notebooks -v $PWD/results:/home/rstudio/results \`
	`scrna-cell-traject`

#### Type the following hyperlink in the browser: 
`http://localhost:8787/`

#### Use the following credentials: 
`User: rstudio`
`Password: Totem`

<br>

---

<br>

