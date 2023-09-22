#!/bin/bash
#SBATCH --job-name="TEST"
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:30:00
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
allocation_policy=spread 

export OMP_NUM_THREADS=$cores
export BLIS_NUM_THREADS=$cores

export OMP_PLACES=$places
export OMP_PROC_BIND=$allocation_policy

#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > spread_mkl.csv
#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > spread_oblas.csv
echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > spread_blis.csv

#echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > spread_dmkl.csv
#echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > spread_doblas.csv
echo "DOUBLE, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > pread_dblis.csv

for size in {2000..20000..1000}
do


	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_mkl.x $size $size $size >> spread_mkl.csv
		#echo -n "${cores}," >> oblas.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_oblas.x $size $size $size >> spread_oblas.csv
		#echo -n "${cores}," >> blis.csv
		numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_blis.x $size $size $size >> spread_blis.csv

		#echo -n "${cores}," >> dmkl.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./dgemm_mkl.x $size $size $size >> spread_dmkl.csv
		#echo -n "${cores}," >> doblas.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./dgemm_oblas.x $size $size $size >> spread_doblas.csv
		#echo -n "${cores}," >> dblis.csv
		numactl --interleave=0,1,2,3,4,5,6,7 ./dgemm_blis.x $size $size $size >> spread_dblis.csv

	done
done


echo REMOVING COMPILED PROGRAMS AND UNLOADING MODULES...
echo
cd ../../..
make clean folder=$folder_save
module purge
