for i in 0.9 3
do
cd swap-dist-$i
CPyMonC<inp
mpiexec -np 32 /SC26/Lacos/CPyMonC/lib/MPI_C_Driver.out
cd ../
done
