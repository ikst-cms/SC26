""" This module is a lite Monte Carlo code with main driver written in C and wrapped by Python. This code is to be used
 with output coming from the PyClex (Cluster Expansion) code to perform finite-T simulation.
 CPyMonC - (C)-based (Py)thon-wrapped (Mo)nte (Ca)rlo - to be pronounced as 'See'-'Pi'-'Monk'
 The operation are a) Calculate energy of any give POSCARs/dir of POSCARs using given ECIs.
 b) Performs a random search through the configurationa space for lowest energy state  and
 c) Finite-T Monte Carlo which has Python-based and C-based (partially parallelized using OpenMP)
 This module is tested with python 3.8.3 at IIT Bhubaneswar."""

import time

from Lacos.CPyMonC.Drivers import ConfigurationSearch, DirectoryCalc
from Lacos.CPyMonC.Energy import Initalize_Correlation_Calculator
from Lacos.CPyMonC.MonteCarlo import MC_Driver
from Lacos.CPyMonC.MPI_handler import MPI_WRITE
from Lacos.CPyMonC.Utility import Load_Template_Structure
from Lacos.CPyMonC.miscfunctions import *


def main():
    """
    This is the main method .
    :return: None.
    """
    os.system( 'clear' )
    write_header(
        "Hello User.... Welcome to Cluster Expansion - application (CPyMonC) Module"
    )
    Start_Time = time.time()
    Project_Name = read_string_input(
        "Give a name for the simulation? (\"Enter\" if not required): ", "Allow Blank"
    )
    Base_Lattice = Load_Template_Structure()
    ClusterList, ECIList, ElementRep = Initalize_Correlation_Calculator( Base_Lattice )
    print(
        "\n What do you want to do? Options are:\n   (a) Calculate energies of structures in dir (POSCAR format)"
    )
    print(
        "(b) Search configuration space for low-energy structures by random-walk (fixed composition and atom "
        "numbers)" )
    print( "   (c) Perfrom atom-exchange Monte Carlo simulation" )
    print( "   (d) Prepare files for MPI based Monte Carlo simulation" )
    ToDo = read_string_input( "\nInput a|b|c|d: " )
    if ToDo.lower() == 'a':
        DirectoryCalc( Base_Lattice, ClusterList, ECIList, ElementRep )
    elif ToDo.lower() == 'b':
        print_warning(
            "This is still Pythonic and very slow. Monte-Carlo (at high-T) is equivalent to pure Random-walk (all "
            "states equally probable)."
        )
        print_warning( "This feature will be discontinued in later version." )
        ConfigurationSearch( Base_Lattice, ClusterList, ECIList, ElementRep )
    elif ToDo.lower() == 'c':
        MC_Driver( Base_Lattice, ClusterList, ECIList, ElementRep )
    elif ToDo.lower() == 'd':
        MPI_WRITE( Base_Lattice, ClusterList, ECIList, ElementRep )
    else:
        exit_script( "Don't understand what you want to do. Exit." )
    print( "\nTotal time taken by the code: %f (sec)\n" % (time.time() - Start_Time) )
    write_header( "Successfully Completed !!!" )
    return


if __name__ == '__main__':
    main()
