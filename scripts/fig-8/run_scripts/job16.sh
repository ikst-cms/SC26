#!/bin/sh
#PBS -N KM
#PBS -V
#PBS -q normal
#PBS -A inhouse
#PBS -l select=1:ncpus=64:mpiprocs=4
#PBS -l walltime=48:00:00

module purge

module load craype-network-opa craype-mic-knl intel/19.1.2 impi/19.1.2
EXE="/home01/r735krs/LACOS/CPyMonC/lib/MPI_C_Driver.out"


cd $PBS_O_WORKDIR


mpirun -np 16 $EXE

# Rename MC.log to core number
if [ -f "MC.log" ]; then
    mv MC.log 16.log
fi

exit 0
