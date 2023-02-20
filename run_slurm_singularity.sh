#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=3-0:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --output=rstudio-server.job.%j

# COMMENT: this script is based on rocker project: 
# https://rocker-project.org/use/singularity

workdir=${PWD}
echo -e "Creating working directory for R session at: ${workdir}\n"
mkdir -p -m 700 ${workdir}/rsessions/run ${workdir}/rsessions/tmp ${workdir}/rsessions/var/lib/rstudio-server
cat > ${workdir}/rsessions/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export R_LIBS_USER=R/rocker-rstudio/4.2.1
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export SINGULARITY_BIND="${workdir}:$HOME,${workdir}/rsessions/run:/run,${workdir}/rsessions/tmp:/tmp,${workdir}/rsessions/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/rsessions/var/lib/rstudio-server:/var/lib/rstudio-server"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export SINGULARITYENV_RSTUDIO_SESSION_TIMEOUT=0

export SINGULARITYENV_USER=$(id -un)
export USER=$SINGULARITYENV_USER
if [[ $1 == "" ]]
then
	export SINGULARITYENV_PASSWORD=$(openssl rand -base64 15)
else
	export SINGULARITYENV_PASSWORD=$1 
fi
export SINGULARITY_CACHEDIR=~/.singularity/cache
# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
#readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
export PORT=8787
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:${HOSTNAME}:${PORT} ${SINGULARITYENV_USER}@<server-ip-address>

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${SINGULARITYENV_USER}
   password: ${SINGULARITYENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END

## Pulling rstudio image if not present
#The following image was built with the def file under the directory: 'cluster_slurm_nb/imgs/single-cell.def' 
#ROCKER=/wrk/local2/B21027_scRNAseq_signatures/261121_tcell_ref/b21027_t_cell_ref/test/cluster_slurm_nb/imgs/single-cell-v1.sif       
ROCKER=imgs/scrna-cell-ti.sif       

echo -e "\nLaunching RStudio...\n"
singularity exec --cleanenv ${ROCKER} \
    /usr/lib/rstudio-server/bin/rserver --server-user ${USER} \
    	    --www-port ${PORT} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh
printf 'rserver exited' 1>&2

