#include <stdio.h>
#define MPI_ELEMENT_REP_FILE "./Output/Element_Rep.txt"
#define MPI_BASE_LATTICE_FILE "./Output/Base_Lattice.txt"
#define MPI_SUPER_LATTICE_FILE "./Output/Super_Lattice.txt"
#define MPI_MC_PARAMETERS_FILE "./Output/MC_Parameters.txt"
#define MPI_CLUSTER_INFO_FILE "./Output/Cluster_Info.txt"
#define MPI_EXCHANGE_DIST_OUTPUT_FILE "./Output/ExchangeDist_OutputFile.txt"

FILE *MPI_Element_Rep_File = NULL; 
FILE *MPI_Base_Lattice_File = NULL;
FILE *MPI_Super_Lattice_File = NULL;
FILE *MPI_MC_Parameters_File = NULL;
FILE *MPI_Cluster_Info_File = NULL;
FILE *MPI_Exchange_Dist_Output_File = NULL;