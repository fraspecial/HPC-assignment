# EXERCISE 2
here we uploaded all the files used to perform measurements for the FHPC project exercise 2

### What does this directory contains


the markdown file: a brief overview of the main files in the directory
gemm.c: a code to call the gemm function and to measure its performance in time and Gflops
Makefile: a Makefile to compile gemm.c with different libraries and precisions
Makefile: a Makefile to compile gemm.c with different libraries and precisions (the pathe to the compiled BLIS library on THIN nodes was added) the various path need to be modified.


size_scalability/: the folder contains data collected to perform the scalability at an increasing number of threadds, it is organized in:
EPYC/
THIN/

each of the previous folder contains the same sub folders each one for a different policy for thread allocation: 
CLOSE/
SPREAD/

the EPYC folder also includes two additional folder for the NUMACTL function :
NUMA_CLOSE/
NUMA_SPREAD/


core_scalability has the same structure, but NUMA_SPREAD folder has not been produced for core scalability. 


NOTEBOOK/: contains the jupyter notebooks used to draw the graphs used to compare performances of different scenarios for CORE SCALABILITY

NOTEBOOK_SIZE/: contains the jupyter notebooks used to draw the graphs used to compare performances of different scenarios for SIZE SCALABILITY
