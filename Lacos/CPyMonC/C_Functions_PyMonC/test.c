#include "mpi.h"
#include "time.h"
#include "stdio.h"
#include "stdlib.h"
#include "MPI_C_MonteCarlo.h"
#include "C_GenFunctions.c"
#include "C_MCmisc.c"
#include "C_MChash.c"
#include "C_MChash_new.c"
#include "MPI_C_MonteCarlo.c"
#include "MPI_C_File.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  FILE *F=NULL;
  F=fopen("./custom", "a");
  pcg32_random_t rng;
  pcg32_srandom_r(&rng,42u,54u);
  for(int i=0;i<5;++i){
    fprintf(F,"\n%d %d %d\n",rank,size,pcg32_boundedrand_r(&rng,5));
  }
  MPI_Finalize();
	return 0;
}