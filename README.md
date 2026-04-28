# Lacos Software/Simulation Package

install_requires= ["python=3.8"]

Install by 

`python -m pip install .`

To Execute

From Command prompt write `DBMaker` to lunch the module similarly for `PyCLEX`, `CPyMonC` and `DBRunner`. 

# CPyMonC

Compilation:

- from the main directory `cd Lacos/CPyMonC/C_Frunctions_PyMonC`
- `chmod +x mpi_compile`
- `./mpi_compile`
- `cd ../../..`

Flags:

- Logs:
    - `TIME`: add this to get time in time in hrs for each MC step
    - `CTIME`: add this to get total time for swaps and number of swaps in file:`custom_log` for `mpi`, `custom_no_op` otherwise
    - `CLIST`: add this to get list sizes in hybrid in file:`custom_list`

- `Hybrid`: to activate Hybrid optimisation
- `OPT_PSITE`: optimise the number of psites from 1 site to two sites.
- `OPT_PRE_CLUSTER`: optimises iteration over clusters

Notes:
- For `OpenMP` add flags in `Setup.py` and run setup again
- For `MPI` add flags in `Lacos/CPyMonC/C_Frunctions_PyMonC` and run the `mpi_compile` script again
- `MPI Driver binary`: `Lacos/CPyMonC/lib/MPI_C_Driver.out`
