cp -r ../../datasets/quaternary domain_size_scaling
cd domain_size_scaling
for i in 1000 10000 100000 1000000
do
    cd $i
    CPyMonC<inp
    mpiexec -np 64 /SC26/Lacos/CPyMonC/lib/MPI_C_Driver.out
    cd ../
done