#!/bin/sh
#PBS -o sys_mesg.log
#PBS -N LACOS

#PBS -l select=1:ncpus=64:mpiprocs=64
#PBS -q full64

cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate mpilacos

ulimit -S unlimited


cp ../../datasets/quaternary domain_size_scaling
cd domain_size_scaling

for i in 1000 10000 100000 1000000
do
    cd $i
    CPyMonC<inp
    mpiexec -np 64 /SC26/Lacos/CPyMonC/lib/MPI_C_Driver.out
    cd ../
done


