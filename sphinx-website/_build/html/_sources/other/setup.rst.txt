Setup
+++++

Setup
=====

Create directory structure: ::

   mkdir data scripts notebooks results

Then, launch the ``repro-totem-ti`` container in the Linux terminal binding the folders created above to the container’s folders with the command: :: 
   
   docker run --rm -ti -e PASSWORD=Totem -p 8787:8787 \
        -v $PWD/data:/home/rstudio/data \
        -v $PWD/scripts:/home/rstudio/scripts \
        -v $PWD/notebooks:/home/rstudio/notebooks \
        -v $PWD/results:/home/rstudio/results \
        compbiomed/repro-totem-ti

Finally, go to the browser of your choice (e.g., Chrome, Firefox) and type the following link: ::
   
   http://localhost:8787

Use the following credentials to login in RStudio Server: :: 
   
   Username: rstudio 
   
   Password: Totem

