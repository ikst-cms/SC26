#!/bin/sh
#PBS -o sys_mesg.log
#PBS -N LACOS

#PBS -l select=1:ncpus=64:mpiprocs=64
#PBS -q full64

cd $PBS_O_WORKDIR

eval "$(conda shell.bash hook)"
conda activate mpilacos

ulimit -S unlimited


mkdir cores_scaling
cd cores_scaling
cp -r ../../../datasets/quaternary/10000 scratch

for i in 1 2 4 8 16 32 64
do
    cp -r scratch core-$i
    cd core-$i
    CPyMonC<inp
    mpiexec -np $i /SC26/Lacos/CPyMonC/lib/MPI_C_Driver.out
    cd ../
done



