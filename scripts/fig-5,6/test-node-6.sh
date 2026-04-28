#!/bin/sh
#PBS -o sys_mesg.log
#PBS -N LACOS

#PBS -l select=6:ncpus=52:mpiprocs=52
#PBS -q full

cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate mpilacos

ulimit -S unlimited

for i in 312
do
	cp -r quaternary-scratch quaternary-core-$i
	cd quaternary-core-$i
	for j in 1000 10000 100000
	do
		cd $j
		CPyMonC<inp
		mpiexec -np $i /SC26/Lacos/CPyMonC/lib/MPI_C_Driver.out
		cd ../
	done
	cd ../
done
       

