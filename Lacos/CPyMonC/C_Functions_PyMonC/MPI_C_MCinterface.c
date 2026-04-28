//-----------------------------------------------------------------------------------------------------+
// Metropolis-based atom-exchange Monte Carlo code written in C and accessed through Python.The code is|
// parallelized using MPI and implements Loop parallelisation.                            |
// 												       |
// This is the main driver and creates interface with Python.					       |
//-----------------------------------------------------------------------------------------------------+

#include "MPI_C_MonteCarlo.h"
#include "C_GenFunctions.c"
#include "C_MCmisc.c"
#include "MPI_C_File.h"
#include <stdio.h>

void First_Entry_and_Receive_Elements (struct Element *, int, char *); 
void Check_for_OpenMP (int *);
void Receive_BaseLattice (struct Labeled_Coord *, int, char *, int, struct Coord *);
void Receive_SupperLattice (struct AtomLabeled_Coord *, int, struct Coord *, int *);
void Receive_MC_Parameters (float *, int, int, float, int, int, int, int);
void Receive_ClusterInfo (struct ClusterFamily *, int, double);
int Create_Neighbour_List_inC (float);
void Create_Perturbed_SiteList_inC ();


void C_MonteCarlo_Driver (char *);
void C_MonteCarlo_Driver__omp (char *);
void Free_Global_Memory_inC (void);


struct Coord Cart2Frac_coordinates (struct Coord *, struct Coord *);


struct tmpAtomLabeled_Coord {
  struct Coord R;
  char *Atom;                        
  int iAtom;                         
  int iBaseLat;  
};


void First_Entry_and_Receive_Elements (struct Element *elementrep, int nelms, char *log_filename)        
{

  if (!(MPI_Element_Rep_File = fopen(MPI_ELEMENT_REP_FILE, "w"))) {
    Exit_Error ("Error in opening element rep file for Monte Carlo (from C).");
  }
  
  int i,j;
  nElements = nelms;
  clock_t tt;
  
  tt = clock();

  // printing logfile name
  fprintf(MPI_Element_Rep_File, "%s\n", log_filename);
  fprintf(MPI_Element_Rep_File, "%d\n", nElements);

  // printing element names
  for (i=0; i<nelms; i++) {
  	fprintf(MPI_Element_Rep_File, "%s ", elementrep[i].Name);
  }
  fprintf(MPI_Element_Rep_File, "\n");

  // printing Theta Rep
  for (i=0; i<nelms; i++) {
  	for (j=0; j<nelms; j++) {
      fprintf(MPI_Element_Rep_File, "%lf ", elementrep[i].Theta[j]);
    }
    fprintf(MPI_Element_Rep_File, "\n");
  }

  tt = clock() - tt;
  if (!(MC_LogFile=fopen(log_filename, "w"))) {
    Exit_Error ("Error in opening log file for Monte Carlo (from C).");
  }
   
  fprintf (MC_LogFile, "\nSuccessfully received elements from Python ...: %.10lf Secs\n",((double)tt)/CLOCKS_PER_SEC);

  fclose(MPI_Element_Rep_File);
  return;
}


void Receive_BaseLattice (struct Labeled_Coord *coordarray, int ncoords, char *labels, int nlabels, struct Coord *latvec) 
{
  if (!(MPI_Base_Lattice_File = fopen(MPI_BASE_LATTICE_FILE, "w"))) {
    Exit_Error ("Error in opening base lattice file for Monte Carlo (from C).");
  }
  int i;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");
  tt = clock();
  
  // printing labels
  fprintf(MPI_Base_Lattice_File, "%d %s\n", nlabels, labels);

  // printing axis
  for (i = 0; i < 3; i++) {
    fprintf(MPI_Base_Lattice_File, "%f %f %f\n", latvec[i].x, latvec[i].y, latvec[i].z);
  }
  
  Base_Lattice.nSites = ncoords;
  // priting site labels and coardinates
  fprintf(MPI_Base_Lattice_File, "%d\n", ncoords);
  for (i=0; i<ncoords; i++){
    fprintf(MPI_Base_Lattice_File, "%s %f %f %f\n", coordarray[i].Label, 
            coordarray[i].R.x, coordarray[i].R.y, coordarray[i].R.z);
  }

  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Base Lattice info from Python... : %.10lf Secs\n", ((double)tt)/CLOCKS_PER_SEC);
  fclose(MPI_Base_Lattice_File);
  return;
}


void Receive_SuperLattice (struct tmpAtomLabeled_Coord *coordarray, int ncoords, struct Coord *latvec, int *ncell, struct CharList *frozen_tmp, int nfrozen)
{
  if (!(MPI_Super_Lattice_File = fopen(MPI_SUPER_LATTICE_FILE, "w"))) {
    Exit_Error ("Error in opening super lattice file for Monte Carlo (from C).");
  }

  int i, j;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");
  tt = clock();

  // printing frozen elements
  fprintf(MPI_Super_Lattice_File, "%d\n", nfrozen);
  for(i=0; i<nfrozen; i++){
    fprintf(MPI_Super_Lattice_File, "%s ", frozen_tmp[i].p);
  }
  fprintf(MPI_Super_Lattice_File, "\n");

  // printing super cell size
  fprintf(MPI_Super_Lattice_File, "%d %d %d\n", ncell[0], ncell[1], ncell[2]);
  
  // printing axis
  for (i=0; i<3; i++) {
    fprintf(MPI_Super_Lattice_File, "%f %f %f\n", latvec[i].x, latvec[i].y, latvec[i].z);
  }

  // printing number if sites
  fprintf(MPI_Super_Lattice_File, "%d\n", ncoords);
    
  for (i=0; i<ncoords; i++) {
    fprintf(MPI_Super_Lattice_File, "%f %f %f %d %s\n", 
      coordarray[i].R.x, coordarray[i].R.y, coordarray[i].R.z,
      coordarray[i].iBaseLat, coordarray[i].Atom
    );
  }

  tt = clock() - tt;
  
  fprintf (MC_LogFile, "Successfully received SuperLattice info from Python; generated hash-table for coordinates...: %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
  fclose(MPI_Super_Lattice_File);
  return;
}


void Receive_MC_Parameters (float *tkseq, int ntkseq, int eqsteps, float eqtol, int nhist, int calcsteps, int blocksize, int min_mc_steps)
{
  if (!(MPI_MC_Parameters_File = fopen(MPI_MC_PARAMETERS_FILE, "w"))) {
    Exit_Error ("Error in opening mc param file for Monte Carlo (from C).");
  }
  int i;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");

  tt = clock();
  fprintf(MPI_MC_Parameters_File, "%d\n", ntkseq);
  for (i=0; i<ntkseq; i++){
    fprintf(MPI_MC_Parameters_File,"%f ", tkseq[i]);
  }
  fprintf(MPI_MC_Parameters_File, "\n");

  fprintf(MPI_MC_Parameters_File, "%d\n", eqsteps);
  fprintf(MPI_MC_Parameters_File, "%f\n", eqtol);
  fprintf(MPI_MC_Parameters_File, "%d\n", calcsteps);
  fprintf(MPI_MC_Parameters_File, "%d\n", nhist);
  fprintf(MPI_MC_Parameters_File, "%d\n", blocksize);
  fprintf(MPI_MC_Parameters_File, "%d\n", min_mc_steps);
  
  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Monte-Carlo parameters from Python...: %.10f Secs \n",((double)tt)/CLOCKS_PER_SEC);
  fclose(MPI_MC_Parameters_File);
  return;
}


void Print_Cluster(FILE* Output_File, struct Cluster curCluster, int Nk)
{
  int i;
  for (i=0; i<Nk; i++) {
    fprintf(Output_File, "%f %f %f\n", curCluster.R[i].x, 
            curCluster.R[i].y, curCluster.R[i].z);
  }
}


void Receive_ClusterInfo (struct ClusterFamily *cflist, int ncfl, double eci0)          
{
  if (!(MPI_Cluster_Info_File = fopen(MPI_CLUSTER_INFO_FILE, "w"))) {
    Exit_Error ("Error in opening cluster info file for Monte Carlo (from C).");
  }

  int i, j, k, cls;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");
                                        
  fprintf(MPI_Cluster_Info_File, "%d\n", ncfl);
  tt = clock();
  for (cls=0; cls<ncfl; cls++){

    fprintf(MPI_Cluster_Info_File, "%d\n", cflist[cls].Ksites);
    fprintf(MPI_Cluster_Info_File, "%s\n", cflist[cls].ClusterID);
    
    fprintf(MPI_Cluster_Info_File, "%d\n", cflist[cls].nECIs);              
    for (j=0; j<cflist[cls].nECIs; j++) {
      fprintf(MPI_Cluster_Info_File, "%lf\n", cflist[cls].ECI_List[j].ECI);

      for (k=0; k<cflist[cls].Ksites; k++) {
        fprintf(MPI_Cluster_Info_File, "%d ", cflist[cls].ECI_List[j].Index[k]);
      }

      fprintf(MPI_Cluster_Info_File, "\n");
    }         

    fflush(MPI_Cluster_Info_File);
    Print_Cluster(MPI_Cluster_Info_File, cflist[cls].Parent, cflist[cls].Ksites);
    fflush(MPI_Cluster_Info_File);

    for (j=0; j<Base_Lattice.nSites; j++) {

      fprintf(MPI_Cluster_Info_File, "%d\n", cflist[cls].ChildGroup_List[j].nChild);
      fprintf(MPI_Cluster_Info_File, "%d\n", cflist[cls].ChildGroup_List[j].BL_Index);

      for (k=0; k<cflist[cls].ChildGroup_List[j].nChild; k++) {
        Print_Cluster(MPI_Cluster_Info_File, cflist[cls].ChildGroup_List[j].Child[k], cflist[cls].Ksites);
      }     
    }
    fflush(MPI_Cluster_Info_File);  
  }
  fprintf(MPI_Cluster_Info_File, "%lf\n", eci0);
  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Cluster and ECI information from Python... : %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
  fclose(MPI_Cluster_Info_File);
  return;
}


void Get_Exchange_Dist_and_Output_File(float ext_minDist, char *output_file, int numList, int listCap){
  if (!(MPI_Exchange_Dist_Output_File = fopen(MPI_EXCHANGE_DIST_OUTPUT_FILE, "w"))) {
    Exit_Error ("Error in opening exchange dist output file for Monte Carlo (from C).");
  }
  fprintf(MPI_Exchange_Dist_Output_File, "%f\n%s\n%d\n%d", ext_minDist, output_file, numList, listCap);
  fclose(MPI_Exchange_Dist_Output_File);
}