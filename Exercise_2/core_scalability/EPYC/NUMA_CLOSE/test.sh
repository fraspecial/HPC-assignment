#!/bin/bash
#SBATCH --job-name="TEST"
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=01:30:00
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

cd $folder_save

partition=EPYC 
cores=128
places=cores
allocation_policy=close 


export OMP_PLACES=$places
export OMP_PROC_BIND=$allocation_policy

#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > mkl.csv
echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > oblas.csv
#echo "FLOAT, ${partition}, ${allocation_policy}, size,time in seconds,GFLOPS" > blis.csv

for cores in {1..16..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores
	numastat -m

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done


for cores in {17..32..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done


for cores in {33..48..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done

for cores in {49..64..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2,3 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2,3 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done

for cores in {65..80..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3,4 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2,3,4 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2,3,4 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done

for cores in {81..96..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x
	numastat -m
	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3,4,5 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2,3,4,5 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2,3,4,5 ./gemm_blis.x 11000 11000 11000 >> blis.csv
		

	done
done

for cores in {97..112..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x

	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3,4,5,6 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2,3,4,5,6 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2,3,4,5,6 ./gemm_blis.x 11000 11000 11000 >> blis.csv

	done
done

for cores in {113..128..1}
do
	export OMP_NUM_THREADS=$cores
	export BLIS_NUM_THREADS=$cores

	gcc /u/dssc/scocom00/git_FHPC/FHP_myrepo/Basic/OpenMP/threads_affinity/00_where_I_am.c -o topo.x -fopenmp
	./topo.x
	numastat -m
	for j in {1..5..1}
	do
		#echo -n "${cores}," >> mkl.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_mkl.x 11000 11000 11000 >> mkl.csv
		echo -n "${cores}," >> oblas.csv
		numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_oblas.x 11000 11000 11000 >> oblas.csv
		#echo -n "${cores}," >> blis.csv
		#numactl --interleave=0,1,2,3,4,5,6,7 ./gemm_blis.x 11000 11000 11000 >> blis.csv
		

	done
done



echo REMOVING COMPILED PROGRAMS AND UNLOADING MODULES...
echo
cd ../../..
make clean folder=$folder_save
module purge
