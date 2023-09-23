# EXERCISE 2
Here we uploaded all the files used to perform measurements for the FHPC project exercise 2.

### What does this directory contains

* the markdown file: a brief overview of the main files in the directory
* gemm.c: a code to call the gemm function and to measure its performance in time and Gflops
* Makefile: a Makefile to compile gemm.c with different libraries and precisions
* Makefile_THIN: a Makefile to compile gemm.c with different libraries and precisions (the pathe to the compiled BLIS library on THIN nodes was added) the various path need to be modified.


* size_scalability/: the folder contains data collected to perform the scalability at an increasing number of threadds, it is organized in:
    * `EPYC/`: contains all the size scalability tests performed on EPYC nodes.
    * `THIN/`: contains all the size scalability tests performed on THIN nodes.

Each of the previous folders contains the same sub folders, each one for a different policy for thread allocation: 
* `CLOSE/`
* `SPREAD/`

The EPYC folder also includes two additional folder for the tests performed using `numactl` command:
* `NUMA_CLOSE/`
* `NUMA_SPREAD/`

`core_scalability/` has the same structure, but `NUMA_SPREAD` folder has not been produced for core scalability. 

`NOTEBOOK/`: contains the jupyter notebooks used to draw the graphs used to compare performances of different scenarios for core scalability.
`NOTEBOOK_SIZE/`: contains the jupyter notebooks used to draw the graphs used to compare performances of different scenarios for size scalability.
