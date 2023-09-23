# Exercise 1
## What you will find in this folder
This folder contains all the source files of the hybrid MPI+OpenMP program which implements Conway's game of life and the bash scripts used to test its correctness and scalability.
List of the contents:
* src/: a folder containing the following source files:
    * init.c: used to initialize MPI parallel region and decompose the domain.
    * read_write.c: read and write functions used for the playground and to record performance.
    * static.c: implements static evolution.
    * ordered.c: implements ordered evolution.

* obj/: folder created to host the objects files when the code is compiled

* Makefile: used to compile the source code

* EPYC: contains all the jobs performed to test scaling on EPYC nodes of Orfeo cluster. Jobs are divided in the subfolders:
    * omp: contains the scripts for OpenMP scaling test of static evolution on 1 socket, 1 node and 2 nodes saturated with threads.
    * strong: contains the scripts for strong MPI scaling of static evolution on 3 nodes with processes mapped by core, and 4 nodes with processes mapped by node and bound to socket. Furthermore it contains also a strong MPi scaling test on ordered evolution.
    * weak: contains the scripts for weak scalability on 4 nodes with processes mapped by node and bound to socket.
* THIN: it has all the jobs performed to test scaling on THIN nodes of Orfeo cluster. The structure of the subdirectories and the filenames are the same of those on EPYC nodes.

* create_image.job: bash script that calls the initialization of the playground.
* evo_image.job: bash script that calls the evolution of the playground.

## Compiling the code and running the scripts.
To compile the code just execute `make` from the home folder of this repository.
To initialize the playground execute `sbatch create_image.job`. The playground will be output in the main directory as a pgm file. It is possible to set the playground size and the output file editing the values of the script's variables. Modify the value of `k`if you prefer creating a square matrix or substitute `x`  and `y` for setting number of rows separately. In this last case modify the command `mpirun` changing `-k $k` with `-x $x -y $y`. Initalized playground is saved in the main directory

```
mpirun ./main.x -i -k $k -f $filename
mpirun ./main.x -i -x $x -y $y -f $filename
```

To evolve the playground execute `sbatch evo_image.job`. Again, modifying the script it is possible to set the evolution kind with `evo`, the number of evolution iterations with `n`, the number of evolution steps to perform between each writing on a file of the partial evolution with `steps`, the file to read with `filename`, and the file where to print number of processes, threads, rows of the playground, number of evolutions and runtime with `output`. Set `evo` to 1 for static evolution or 0 to for ordered one. Evolved playgrounds are saved in the subdirectory of the main directory `evos/`.
To run scaling tests, all of them must be run from the main directory related to the architecture. All the paths are of the kind `<node_type>/jobs/{omp|strong|omp}/<script_name>`. Filenames are rather explicative of the kind of scaling testes. So, if for example you wich to run a strong MPI-scaling test on THIN nodes run:
```
cd THIN
sbatch jobs/strong/str_3n.job
```
The output file produced will be saved in the same folder of the script. Runtimes are saved in the `res/` folder located inside `EPYC` or `THIN`, depending on the architecture you run, on a csv file where on each line you can find:
```
#processes,#threads,#rows,#evolutions,#runtime
```
We warmly invite you to read the Report in the home of this repository for further information.