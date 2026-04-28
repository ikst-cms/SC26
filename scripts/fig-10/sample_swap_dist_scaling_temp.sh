#!/bin/bash

# Check if env variable is set
if [ -z "$CPYMONC_MPI_OUT" ]; then
    echo "Error: CPYMONC_MPI_OUT is not set."
    echo "Please compile Lacos/CPyMonC/C_Functions_PyMonC/MPI_C_Driver.c, by Lacos/CPyMonC/C_Functions_PyMonC/mpi_compile and export it, eg:"
    echo "export CPYMONC_MPI_OUT=/path/to/Lacos/CPyMonC/lib/MPI_C_Driver.out"
    exit 1
fi

for i in 0.9; do
    cd "swap-dist-$i-temp" || { echo "Directory swap-dist-$i-temp not found"; exit 1; }

    CPyMonC < inp

    mpiexec -np 32 "$CPYMONC_MPI_OUT"

    cd .. || exit 1
done