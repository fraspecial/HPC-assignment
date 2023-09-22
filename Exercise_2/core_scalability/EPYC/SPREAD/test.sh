#!/bin/bash
#SBATCH --job-name="TEST"
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --partition=EPYC
#SBATCH --error=job.%J.err
###SBATCH --output=job.%J.out


module load architecture/AMD
module load mkl
module load openBLAS/0.3.23-omp
export LD_LIBRARY_PATH=/u/dssc/scocom00/myblis/lib:$LD_LIBRARY_PATH

folder_save=$(pwd)
cd ../../..

make float folder=$folder_save
make double folder=$folder_save

cd $folder_save

partition=EPYC 
places=cores
allocation_policy=spread

export OMP_PLACES=$places
export OMP_PROC_BIND=$allocation_policy


#echo "FLOAT, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > mkl.csv
#echo "FLOAT, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > oblas.csv
echo "FLOAT, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > blis.csv

#echo "DOUBLE, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > dmkl.csv
#echo "DOUBLE, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > doblas.csv
echo "DOUBLE, ${partition}, ${allocation_policy}, \n cores,size,time in seconds,GFLOPS" > dblis.csv

for cores in {1..128..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		#echo -n "${cores}," >> oblas.csv
		#./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		echo -n "${cores}," >> blis.csv
		./gemm_blis.x 11000 11000 11000 >> blis.csv

		#echo -n "${cores}," >> dmkl.csv
		#./dgemm_mkl.x 11000 11000 11000 >> dmkl.csv
		#echo -n "${cores}," >> doblas.csv
		#./dgemm_oblas.x 11000 11000 11000 >> doblas.csv
		echo -n "${cores}," >> dblis.csv
		./dgemm_blis.x 11000 11000 11000 >> dblis.csv

	done
done


echo REMOVING COMPILED PROGRAMS AND UNLOADING MODULES...
echo
cd ../../..
make clean folder=$folder_save
module purge
