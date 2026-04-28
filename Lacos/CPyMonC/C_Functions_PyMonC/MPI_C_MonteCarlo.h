//-----------------------------------------------------------------------------------------------------+
// Header for the Monte Carlo code written in C and accessed through Python. The code is parallelized  |
// using OpenMP and implements Domain Decomposition technique.                                         |
//-----------------------------------------------------------------------------------------------------+

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pcg_basic.c"



#define TINY 0.0
#define LARGE 1.0e23

#define COORD_TOL 1.0e-3
#define CLUSTER_TOL 1.0e-5 


#define COULOMB 1.60217646e-19 
#define ANGS 1.0e-10
#define eV2kJ  96.48538
#define mJFac 1000

#define PI 3.1415926535897932384626433832795
#define kBOLTZMANN 8.617343e-5              
#define kBOLTZMANN_kJmol 0.00831447         


float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

float signarga, signargb;
#define SIGN(signarga,signargb) ((signargb) >= 0.0 ? fabs(signarga): -fabs(signarga))


struct CharList {
  char *p;
};

struct Stat {
  double Mean, Median, StdDev;
  double SpHeat, dSpHeat;
};


struct Coord {
  float x, y, z;
};

struct IntCoord {
  int x, y, z;
};

struct Labeled_Coord {
  struct Coord R;
  char *Label;
};


struct AtomLabeled_Coord {
  struct Coord R;
  int iAtom;                         
  int iBaseLat;  
  int Active;                        
  
  int *Site_Neighbors;               
  int nNeighbors;
  
  int *Perturbed_Sites;
  int nPerturbed_Sites;

  int *correspondingCluster;
};

struct Coord origin={0.0,0.0,0.0};   
struct Coord DiagVec={1.0,1.0,1.0};  

struct Element {  
  char *Name;                        
  double *Theta;                     
};                                   

int nElements;                       
struct Element *ElementList=NULL;          

int nFrozen_Elements;
int *Frozen_Elements=NULL;

struct UnitCell {
  struct Labeled_Coord *Sites;
  struct Coord Lattice_Vectors[3];
  int nSites;

  char *Labels;
  int nLabels;
};

struct UnitCell Base_Lattice;        
struct SuperCell {
  struct AtomLabeled_Coord *Sites;   
  int nSites;                        
  int *nAtoms;                       

  struct Coord Lattice_Vectors[3];    
  struct Coord Ncell;                

  int nActiveSites;
};

struct SuperCell MCSystem;           

float *TKelvin_Sequence=NULL;
int nTKelvin_Sequence;

int Equilibrium_Steps, Calculation_Steps, BLOCK_SIZE;
float Equilibrium_Tolerance;
int NHISTORY, minEQUIL_STEPS;
int DEFAULT_EQUILIBRIUM_STEPS=1000000;

FILE *MC_LogFile=NULL;               


struct Cluster {                     
  struct Coord *R;                   
};

struct ChildGroup {
  struct Cluster *Child;
  int nChild, BL_Index;
};

struct Indexed_ECI {
  double ECI;
  int *Index;
};

struct ClusterFamily {
  int Ksites;
  char *ClusterID;
  
  struct Cluster Parent;                 
  struct ChildGroup *ChildGroup_List;   
  int nChildGroups;                     
  
  struct Indexed_ECI *ECI_List;
  int nECIs;
};

int nClusters;
struct ClusterFamily *ClusterList=NULL;

int nCE_TERMS;
double *ECI_Vec=NULL;

int LargestKsites, maxClusterComponents; 

int Use_Local_EnergyDiff_Algo;
float Perturbation_Distance;                
float maxAtom_Exchng_Dist;
int Perturbed_Sites_List_Allocated=0;       


int OpenMP_ON=0;                            
int Domain_Decomposition=0;                 
int Loop_Parallelization=0;                

struct IntCoord nSuperBlocks, nSubBlocks;   
int Total_nSuperBlocks, Total_nSubBlocks;   
struct IntCoord dSuperBlock, dSubBlock;     

int grank, gsize;
int *mpi_total_psites, max_psites;

// psp -> partial sequence preserving algorithm

int *psp_labels;
int psp_curList,psp_maxList,psp_maxCap,psp_totalList,psp_actCap;
int **psp_list_jump,**psp_list_neigh,*psp_act_jump,*psp_act_neigh;
int *psp_curCap;
double *psp_local_energy,*psp_final_energy;
int psp_max_psites;
int *psp_act_psites,*psp_swap_ind;

pcg32_random_t rng,rng_select;

long double custom_swaps,custom_time;