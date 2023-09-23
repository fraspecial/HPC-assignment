# FHPC final exam project - course 2022/2023
In this repository we (Sara Cocomello and Francesco Speciale) present our work for the Foundations of High Performance Computing final exam course of  Data Science and Scientific Computing master degree at University of Trieste.

## What you will find in this repository
* Report.pdf: a detailed document about everything we did and why.
* `Exercise1/`: a directory containing a implementation in C language of Conway's game of life (https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life) for two different possible kinds of evolution. The first one is the so called "static", where the status of all the cells is evaluated at first, freezing the system, and updating the status of thecells happens only afterwards. The second one is the "ordered" evolution, where cells are evaluated and updated one by one in row-major order. Since code is realized using MPI and OpenMP APIs and the task is highly parallelizable for its own nature, various scaling test were executed on its performance. We tested OpenMP scalability, strong MPI scalability and weak MPI scalability.
* `Exercise2/`: a directory containing all files related to the second exercise, where we were asked to 
to compare the performance of three HPC math libraries (MKL, openBLAS and BLIS) on a level 3 BLAS function called gemm, which performs matrix-matrix multiplications. Also here we tested scaling, for increasing matrix size (at fixed number of CPUs) and for incresing number of CPUs (at fixed matrix size). We did this both for single and double point precision floating point number operations and with different threads allocation policies. Regarding BLIS library, that had to be downloaded and compiled by the student.

All the scalability studies were done running the codes on ORFEO, the cluster hosted at Area Science Park, Trieste, on two different partitions, EPYC (8 nodes available, with 2 sockets each and 64 cores per socket) and THIN (9 nodes available, with 2 sockets each and 12 cores per socket).

For further details we invite you to read the report and the README files inside the two folders.
