#!/bin/bash
#SBATCH --job-name="TEST"
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:10:00
#SBATCH --exclusive
#SBATCH --partition=EPYC
#SBATCH --error=job.%J.err
###SBATCH --output=job.%J.out


module load architecture/AMD
module load mkl
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/scocom00/myblis/lib:$LD_LIBRARY_PATH
module load hwloc


folder_save=$(pwd)
cd ../../..

make float folder=$folder_save
make double folder=$folder_save

cd $folder_save

partition=EPYC 
cores=64
places=cores
allocation_policy=close 

export OMP_NUM_THREADS=$cores
export BLIS_NUM_THREADS=$cores

export OMP_PLACES=$places
export OMP_PROC_BIND=$allocation_policy

#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > mkl.csv
#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > oblas.csv
echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > stat_blis.csv

#echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > dmkl.csv
#echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > doblas.csv
#echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > stat_dblis.csv

for size in {2000..20000..1000}
do


	for j in {1..5..1}
	do

		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3 ./gemm_mkl.x $size $size $size >> mkl.csv
		#echo -n "${cores}," >> oblas.csv
		#numactl --interleave=0,1,2,3 ./gemm_oblas.x $size $size $size >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		numactl --interleave=0,1,2,3 ./gemm_blis.x $size $size $size >> stat_blis.csv
		numastat -m

		#echo -n "${cores}," >> dmkl.csv
		#numactl --interleave=0,1,2,3 ./dgemm_mkl.x $size $size $size >> dmkl.csv
		#echo -n "${cores}," >> doblas.csv
		#numactl --interleave=0,1,2,3 ./dgemm_oblas.x $size $size $size >> doblas.csv
		#echo -n "${cores}," >> dblis.csv
		#numactl --interleave=0,1,2,3 ./dgemm_blis.x $size $size $size >> stat_dblis.csv

	done
done


echo REMOVING COMPILED PROGRAMS AND UNLOADING MODULES...
echo
cd ../../..
make clean folder=$folder_save
module purge
