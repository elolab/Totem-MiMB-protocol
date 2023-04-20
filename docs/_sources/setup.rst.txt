Setup
+++++

Setup
=====

Create directory structure: ::

   mkdir -p repro-totem-ti/results

   cd repro-totem-ti

Then, launch the ``repro-totem-ti`` container in the Linux terminal binding the `results` folder created above to the container’s folder with the command: :: 
   
   docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \
        -v $PWD/results:/home/rstudio/results \
        elolab/repro-totem-ti

Finally, go to the browser of your choice (e.g., Chrome, Firefox) and type the following link: ::
   
   http://localhost:8787

Use the following credentials to login in RStudio Server: :: 
   
   Username: rstudio 
   
   Password: Totem


.. figure:: gifs/setup.gif

