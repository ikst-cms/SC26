//-----------------------------------------------------------------------------------------------------+
// Metropolis-based atom-exchange Monte Carlo code written in C and accessed through Python.The code is|
// parallelized using OpenMP and implements Domain Decomposition technique.                            |
// 												       |
// This is the main driver and creates interface with Python.					       |
//-----------------------------------------------------------------------------------------------------+

#include "C_MonteCarlo.h"
#include "C_GenFunctions.c"
#include "C_MCmisc.c"
#include "C_MChash.c"
#include "C_MChash_new.c"
#include "time.h"

#ifdef _OPENMP
#include "C_DD_MonteCarlo.c"
#endif

#include "C_MonteCarlo.c"


void First_Entry_and_Receive_Elements (struct Element *, int, char *); 
void Check_for_OpenMP (int *);
void Receive_BaseLattice (struct Labeled_Coord *, int, char *, int, struct Coord *);
void Receive_SupperLattice (struct AtomLabeled_Coord *, int, struct Coord *, int *);
void Receive_MC_Parameters (float *, int, int, float, int, int, int, int);
void Receive_ClusterInfo (struct ClusterFamily *, int, double);
int Create_Neighbour_List_inC (float);
int Create_Perturbed_SiteList_inC (int);



void C_MonteCarlo_Driver (char *);
void C_MonteCarlo_Driver__omp (char *);
void Free_Global_Memory_inC (void);


void Output_AtomCoordinates_poscar (float);
void Output_CheckPoint_File (float);
struct Coord Cart2Frac_coordinates (struct Coord *, struct Coord *);


struct tmpAtomLabeled_Coord {
  struct Coord R;
  char *Atom;                        
  int iAtom;                         
  int iBaseLat;  
};



void First_Entry_and_Receive_Elements (struct Element *elementrep, int nelms, char *log_filename)        
{
  int i, j;
  clock_t tt;
  nElements = nelms;
 
  if (!(ElementList=(struct Element *) malloc (nElements*sizeof(struct Element)))) 
    Exit_Error ("Unable to allocate ElementList.");

  tt = clock();
  for (i=0; i<nElements; i++) ElementList[i].Name  = alloc_1Dchar_array(2); 
  for (i=0; i<nElements; i++) ElementList[i].Theta = alloc_1Ddouble_array(nElements);
    
  for (i=0; i<nElements; i++) {
  	for (j=0; j<strlen(elementrep[i].Name)+1; j++) ElementList[i].Name[j] = elementrep[i].Name[j];	
  	for (j=0; j<nElements; j++) ElementList[i].Theta[j] = elementrep[i].Theta[j];
  }

  tt = clock() - tt;
  if (!(MC_LogFile=fopen(log_filename, "w"))) Exit_Error ("Error in opening log file for Monte Carlo (from C).");
  
  print_header ("Starting Monte Carlo in C", MC_LogFile);  
  fprintf (MC_LogFile, "\nSuccessfully received elements from Python ...: %.10lf Secs\n",((double)tt)/CLOCKS_PER_SEC);
  
  return;
} 


void Check_for_OpenMP (int *pBlocks)         
{
  int ncores, nthreads;

  if (pBlocks[0] == 0 && pBlocks[1] == 0 && pBlocks[2] == 0) { 
#ifdef _OPENMP                            
    ncores = omp_get_num_procs ();        
    nthreads = omp_get_max_threads ();    
#else
    ncores = 1;
    nthreads = 1;
#endif
    Max_OpenMP_Threads = ncores;          
    pBlocks[0] = ncores;
    pBlocks[1] = nthreads;
    return;}
  else {
    OpenMP_ON = 1;

    if (pBlocks[0] != -1 && pBlocks[1] != -1 && pBlocks[2] != -1)  {
      Domain_Decomposition = 1;
      nSuperBlocks.x = pBlocks[0];
      nSuperBlocks.y = pBlocks[1];
      nSuperBlocks.z = pBlocks[2];}
    else Loop_Parallelization=1;

    return;
  }

}



void Receive_BaseLattice (struct Labeled_Coord *coordarray, int ncoords, char *labels, int nlabels, struct Coord *latvec) 
{
  int i;
  clock_t tt;

  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");
  tt = clock();
  Base_Lattice.nSites = ncoords;
  Base_Lattice.nLabels =nlabels;
  Base_Lattice.Labels = alloc_1Dchar_array (Base_Lattice.nLabels+1);
  strcpy (Base_Lattice.Labels, labels);
  for (i=0; i<3; i++) Base_Lattice.Lattice_Vectors[i] = (struct Coord) {.x=latvec[i].x, .y=latvec[i].y, .z=latvec[i].z};
  
  Base_Lattice.Sites = (struct Labeled_Coord *) malloc (Base_Lattice.nSites*sizeof(struct Labeled_Coord));
  for (i=0; i<Base_Lattice.nSites; i++) Base_Lattice.Sites[i].Label = alloc_1Dchar_array (2);
  
  for (i=0; i<Base_Lattice.nSites; i++) {
    Base_Lattice.Sites[i].R = (struct Coord) {.x=coordarray[i].R.x, .y=coordarray[i].R.y, .z=coordarray[i].R.z};
    strcpy (Base_Lattice.Sites[i].Label, coordarray[i].Label);
  }
  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Base Lattice info from Python... : %.10lf Secs\n", ((double)tt)/CLOCKS_PER_SEC);

  return;
}


void Receive_SuperLattice (struct tmpAtomLabeled_Coord *coordarray, int ncoords, struct Coord *latvec, int *ncell, struct CharList *frozen_tmp, int nfrozen)
{
  int i, j;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");
  tt = clock();
  nFrozen_Elements = nfrozen;
  if (nFrozen_Elements > 0) {
    Frozen_Elements = alloc_1Dint_array (nFrozen_Elements); 
    for (i=0; i<nFrozen_Elements; i++) 
      for (j=0; j<nElements; j++) if (strcmp (frozen_tmp[i].p, ElementList[j].Name) == 0) Frozen_Elements[i] = j;
  }
 
  MCSystem.Ncell = (struct Coord) {.x=(float)(ncell[0]), .y=(float) (ncell[1]), .z=(float) (ncell[2])};
  
  for (i=0; i<3; i++) MCSystem.Lattice_Vectors[i] = (struct Coord) {.x=latvec[i].x, .y=latvec[i].y, .z=latvec[i].z};
    
  MCSystem.nAtoms = alloc_1Dint_array (nElements); 
  for (i=0; i<nElements; i++) MCSystem.nAtoms[i]=0;
  MCSystem.nSites = ncoords;
  MCSystem.Sites = (struct AtomLabeled_Coord *) malloc (MCSystem.nSites*sizeof(struct AtomLabeled_Coord));
  MCSystem.nActiveSites=0;
    
  for (i=0; i<MCSystem.nSites; i++) {
    MCSystem.Sites[i].R = (struct Coord) {.x=coordarray[i].R.x, .y=coordarray[i].R.y, .z=coordarray[i].R.z};
    MCSystem.Sites[i].iBaseLat = coordarray[i].iBaseLat;
    
    for (j=0; j<nElements; j++) 
      if (strcmp (coordarray[i].Atom, ElementList[j].Name) == 0) {
        MCSystem.Sites[i].iAtom = j; 
        MCSystem.nAtoms[j]++;
        
        if (nFrozen_Elements > 0)
          if (exist_in_list (Frozen_Elements, nFrozen_Elements, j) == -1) MCSystem.Sites[i].Active = 1;
          else MCSystem.Sites[i].Active = 0;
        else MCSystem.Sites[i].Active = 1;    
        break;
      }
    
    if (MCSystem.Sites[i].Active) {                             
      MCSystem.nActiveSites++;                                 
      MCSystem.Sites[i].nPerturbed_Sites = 0;}                  
    else {
      MCSystem.Sites[i].nPerturbed_Sites = -1;                  
    }
  }
  tt = clock() - tt;

#ifdef OLD_HASH_KEY  
  Generate_Old_HashTable (&MCSystem, &Base_Lattice);      
#else
  Generate_New_HashTable (&MCSystem);                                   
#endif
  
  fprintf (MC_LogFile, "\nFound %d active sites out of %d sites in the system.\n", MCSystem.nActiveSites, MCSystem.nSites);
  fprintf (MC_LogFile, "Successfully received SuperLattice info from Python; generated hash-table for coordinates...: %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
  
  return;
}



void Receive_MC_Parameters (float *tkseq, int ntkseq, int eqsteps, float eqtol, int nhist, int calcsteps, int blocksize, int min_mc_steps)
{
  int i;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");

  tt = clock();
  nTKelvin_Sequence = ntkseq;
  TKelvin_Sequence = alloc_1Dfloat_array (nTKelvin_Sequence);
  for (i=0; i<nTKelvin_Sequence; i++) TKelvin_Sequence[i] = tkseq[i];

  Equilibrium_Steps = eqsteps;
  Equilibrium_Tolerance = eqtol;
  Calculation_Steps =calcsteps;
  NHISTORY = nhist;
  BLOCK_SIZE = blocksize;
  minEQUIL_STEPS = min_mc_steps;
  
  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Monte-Carlo parameters from Python...: %.10f Secs \n",((double)tt)/CLOCKS_PER_SEC);
  return;
}



void Receive_ClusterInfo (struct ClusterFamily *cflist, int ncfl, double eci0)          
{
  int i, j, k, cls;
  clock_t tt;
  if (MC_LogFile == NULL) Exit_Error ("Incorrect entry to C module of Monte Carlo. Check");

  nClusters = ncfl;                                                                    
  if (!(ClusterList = (struct ClusterFamily *) malloc (sizeof(struct ClusterFamily)*(nClusters)))) Exit_Error ("Unable to allocate memory for Cluster List.");

  nCE_TERMS=1;
  tt = clock();
  for (cls=0; cls<ncfl; cls++){
    ClusterList[cls].Ksites = cflist[cls].Ksites;
    ClusterList[cls].ClusterID = alloc_1Dchar_array (15);
    strcpy (ClusterList[cls].ClusterID, cflist[cls].ClusterID);
    
    ClusterList[cls].nECIs = cflist[cls].nECIs;                          
    ClusterList[cls].ECI_List = (struct Indexed_ECI *) malloc (sizeof(struct Indexed_ECI)*cflist[cls].nECIs);
    for (j=0; j<ClusterList[cls].nECIs; j++) {
      ClusterList[cls].ECI_List[j].ECI = cflist[cls].ECI_List[j].ECI;   

      ClusterList[cls].ECI_List[j].Index = alloc_1Dint_array (ClusterList[cls].Ksites);    
      for (k=0; k<ClusterList[cls].Ksites; k++)                          
        ClusterList[cls].ECI_List[j].Index[k] = cflist[cls].ECI_List[j].Index[k]; 
    }
    nCE_TERMS += ClusterList[cls].nECIs;
    ClusterList[cls].Parent.R = alloc_1DCoord_array (ClusterList[cls].Ksites);
    Copy_Cluster (&ClusterList[cls].Parent, &cflist[cls].Parent, ClusterList[cls].Ksites);

    ClusterList[cls].nChildGroups = Base_Lattice.nSites;
    ClusterList[cls].ChildGroup_List = (struct ChildGroup *) malloc(sizeof(struct ChildGroup)*ClusterList[cls].nChildGroups);
    for (j=0; j<ClusterList[cls].nChildGroups; j++) {
      ClusterList[cls].ChildGroup_List[j].nChild = cflist[cls].ChildGroup_List[j].nChild;
      ClusterList[cls].ChildGroup_List[j].BL_Index = cflist[cls].ChildGroup_List[j].BL_Index;
      
      ClusterList[cls].ChildGroup_List[j].Child = (struct Cluster *) malloc(sizeof(struct Cluster)*ClusterList[cls].ChildGroup_List[j].nChild);
      for (k=0; k<ClusterList[cls].ChildGroup_List[j].nChild; k++) {
        ClusterList[cls].ChildGroup_List[j].Child[k].R = alloc_1DCoord_array (ClusterList[cls].Ksites);
        Copy_Cluster (&ClusterList[cls].ChildGroup_List[j].Child[k], &cflist[cls].ChildGroup_List[j].Child[k], ClusterList[cls].Ksites);
      }     
    }  
  }
  tt = clock() - tt;
  fprintf (MC_LogFile, "\nSuccessfully received Cluster and ECI information from Python... : %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
 
  ECI_Vec = alloc_1Ddouble_array (nCE_TERMS);
  ECI_Vec[0] = eci0;
  for (i=1, cls=0; cls<nClusters; cls++) 
    for (j=0; j<ClusterList[cls].nECIs; j++, i++) ECI_Vec[i] = ClusterList[cls].ECI_List[j].ECI;
 
  for (cls=0, LargestKsites=1; cls<nClusters; cls++) 
    if (ClusterList[cls].Ksites > LargestKsites) LargestKsites = ClusterList[cls].Ksites;

  for (cls=0, maxClusterComponents=1; cls<nClusters; cls++) 
    if (ClusterList[cls].nECIs > maxClusterComponents) maxClusterComponents = ClusterList[cls].nECIs;

  return;
}


int Create_Neighbour_List_inC (float ext_minDist)
{
  int i, j;
  float *minDist=NULL, pdist, max_minDist;
  int **tmp_NNb, *tmp_nNNb=NULL;
  int MAX_NNB, avgNNB, n;
  time_t tt;
  tt = clock();
 
  minDist = alloc_1Dfloat_array (Base_Lattice.nSites);
  if ((ext_minDist) == -1.0) {
    if (Base_Lattice.nSites == 1) minDist[0]=1.0;
    else {
      for (i=0; i<Base_Lattice.nSites; i++) minDist[i]=LARGE;

      for (i=0; i<Base_Lattice.nSites; i++) {
        for (j=i+1; j<Base_Lattice.nSites; j++) {
          pdist = Periodic_Distance (&Base_Lattice.Sites[i].R, &Base_Lattice.Sites[j].R, &DiagVec); 
          if (pdist < minDist[i]) minDist[i] = pdist;
          if (pdist < minDist[j]) minDist[j] = pdist;
        }
      }
    }}
  else 
    for (i=0; i<Base_Lattice.nSites; i++) minDist[i]=(ext_minDist);
  
  max_minDist=0.0;
  for (i=0; i<Base_Lattice.nSites; i++) if (minDist[i] > max_minDist) max_minDist=minDist[i];
  MAX_NNB = 10 + (int) (pow(2*max_minDist,3.0) * Base_Lattice.nSites);             
  tmp_NNb = alloc_2Dint_array (MCSystem.nSites, MAX_NNB);
  tmp_nNNb = alloc_1Dint_array (MCSystem.nSites);         
  for (i=0; i<MCSystem.nSites; i++) tmp_nNNb[i]=0;

  for (i=0; i<MCSystem.nSites; i++) {
    if (MCSystem.Sites[i].Active) {                                             

      for (j=i+1; j<MCSystem.nSites; j++) {
        if (MCSystem.Sites[j].Active) {                                         
          if (Periodic_Distance (&MCSystem.Sites[i].R, &MCSystem.Sites[j].R, &MCSystem.Ncell) < (minDist[MCSystem.Sites[i].iBaseLat] + COORD_TOL)) {
            tmp_NNb[i][tmp_nNNb[i]]=j; tmp_nNNb[i]+=1;                          
            tmp_NNb[j][tmp_nNNb[j]]=i; tmp_nNNb[j]+=1;                          
          }
        }
      }}
    else 
      tmp_nNNb[i]=-1;
  }
  tt = clock() - tt;
 
  for (i=0; i<MCSystem.nSites; i++) 
    if (tmp_nNNb[i] == 0) {
      free (minDist); free (tmp_nNNb); 
      free_2Dint_array (tmp_NNb, MCSystem.nSites);
      return 0;
    }
 
  avgNNB=0;
  for (i=n=0; i<MCSystem.nSites; i++) {
    if (tmp_nNNb[i] == -1) 
      MCSystem.Sites[i].nNeighbors=0;                                                    
    else {
      MCSystem.Sites[i].nNeighbors=tmp_nNNb[i];                                           
      avgNNB+=tmp_nNNb[i]; n+=1;
      
      MCSystem.Sites[i].Site_Neighbors = alloc_1Dint_array (tmp_nNNb[i]);         
      for (j=0; j<tmp_nNNb[i]; j++) MCSystem.Sites[i].Site_Neighbors[j] = tmp_NNb[i][j];  
    }
  }
  
  free (minDist); 
  free (tmp_nNNb); 
  free_2Dint_array (tmp_NNb, MCSystem.nSites);
  avgNNB/= (float) (n);
  printf ("\n  Average number of neighbors (from C-code): %d\n", (int)(avgNNB));
  
  fprintf (MC_LogFile, "\nSuccessfully created nearest-neighbor list for all active coordinates... : %.10lf Secs \n", ((double)tt)/CLOCKS_PER_SEC);
  fprintf (MC_LogFile, "Average number of neighbors (from C-code): %d\n", (int)(avgNNB));

  maxAtom_Exchng_Dist = max_minDist;           

  return 1; 
}



int Create_Perturbed_SiteList_inC (int visit)
{
  int i, j, p;
  struct Coord Rdiff;
  float pnorm=0.0, d;
  int nMax_PerturbedSites, *tmparray=NULL;
  clock_t tt;
  tt = clock();

  if (visit == 1) {                                       
    for (i=0; i<nClusters; i++)                          
      if (ClusterList[i].Ksites > 1) 
        for (j=0; j<(ClusterList[i].Ksites-1); j++) {
          Rdiff = Subtract_Coords (&ClusterList[i].Parent.R[j], &ClusterList[i].Parent.R[j+1]);
          if ((d=Norm3d(&Rdiff)) > pnorm) pnorm=d; 
        }

#ifdef OPT_PSITE
    Perturbation_Distance = pnorm*1.01;
#else     
    Perturbation_Distance = (pnorm + maxAtom_Exchng_Dist)*1.2;
#endif

    fprintf (MC_LogFile, "\n  Perturbation length uesd in calculating energy difference: %4.2f\n", Perturbation_Distance);  fflush(MC_LogFile);
    printf ("\n  Perturbation length uesd in calculating energy difference: %4.2f\n", Perturbation_Distance);

    fprintf (MC_LogFile, "  Attempting to find number of sites in perturbed region. If more than 60%% of total sites, global energy calculation to be used for Metropolis step.\n"); fflush(MC_LogFile);
    printf ("  Attempting to find number of sites in perturbed region. If more than 60%% of total sites, global energy calculation to be used for Metropolis step.\n");
    
    nMax_PerturbedSites=0;  
    for (i=0; i<Base_Lattice.nSites; i++) {             
      p=0;
      if (MCSystem.Sites[i].Active) {
        for (j=0; j<MCSystem.nSites; j++) 
          if (Periodic_Distance (&MCSystem.Sites[i].R, &MCSystem.Sites[j].R, &MCSystem.Ncell) < Perturbation_Distance) p++;
      }
      if (p > nMax_PerturbedSites) nMax_PerturbedSites=p;
    }

#ifdef OPT_PSITE
    total_psites = alloc_1Dint_array(nMax_PerturbedSites);
#endif

    if (nMax_PerturbedSites > (int) (0.6*MCSystem.nSites)) Use_Local_EnergyDiff_Algo=0;      
    else Use_Local_EnergyDiff_Algo=1;
    fprintf (MC_LogFile, "  Approximate Number of sites within perturbed region based on largest cluster in CE (compare it with 60%% of total sites): %d (%d)\n", nMax_PerturbedSites, (int) (0.6*MCSystem.nSites)); fflush(MC_LogFile);
    printf ("  Approximate Number of sites within perturbed region based on largest cluster in CE (compare it with 60%% of total sites): %d (%d)\n", nMax_PerturbedSites, (int) (0.6*MCSystem.nSites));

    if (Use_Local_EnergyDiff_Algo) {
     tmparray = alloc_1Dint_array (nMax_PerturbedSites);             
   
      for (i=0; i<MCSystem.nSites; i++) {
        if (MCSystem.Sites[i].Active) {
          for (j=0; j<MCSystem.nSites; j++) 
            if (Periodic_Distance (&MCSystem.Sites[i].R, &MCSystem.Sites[j].R, &MCSystem.Ncell) < Perturbation_Distance) 
              tmparray[MCSystem.Sites[i].nPerturbed_Sites++]=j;
           
          MCSystem.Sites[i].Perturbed_Sites = alloc_1Dint_array(MCSystem.Sites[i].nPerturbed_Sites);               
          for (j=0; j<MCSystem.Sites[i].nPerturbed_Sites; j++) MCSystem.Sites[i].Perturbed_Sites[j] = tmparray[j];}            
        else
          MCSystem.Sites[i].nPerturbed_Sites=0;

      }
      tt = clock() - tt;
      free (tmparray);
    
      fprintf (MC_LogFile, "\n   Using Local calculation to find energy difference on atom-exchange (with perturbed sites stored in memory)... : %f Secs \n", ((double)tt)/CLOCKS_PER_SEC);  
      fflush(MC_LogFile);
      printf ("\n   Using Local calculations to find energy difference on atom-exchange (with perturbed sites stored in memory)\n");
      
      Perturbed_Sites_List_Allocated =1;
      return 1;}
    else {
      fprintf (MC_LogFile, "  Perturbed region is comparable to system size. Using full system for calculating energy difference.\n");
      fflush(MC_LogFile);
      printf ("  Perturbed region is comparable to system size. Using full system for calculating energy difference.\n");
      return 0;
    }
  }


  if (visit == 2) {                                      
    for (i=0; i<MCSystem.nSites; i++) 
      if (MCSystem.Sites[i].Active) {
        free (MCSystem.Sites[i].Perturbed_Sites);        
        MCSystem.Sites[i].nPerturbed_Sites = -1;         
      }
  
    Perturbed_Sites_List_Allocated = 0;    
    fprintf (MC_LogFile, "\n   Deallocated memory for perturbed-sites.\n");  fflush(MC_LogFile);
    printf ("\n   Deallocated memory for perturbed-sites.\n");
    return (1);
  }

  return -1;
}



void C_MonteCarlo_Driver (char *out_filename)
{
  int i, EQUIL_STEPS_ACHIEVED, EquilSteps;
  struct Stat EnergyStat;
  float TkB, EQUIL_TOL_ACHIEVED;
  FILE *MC_OutputFile=NULL;
  
#if _OPENMP
  if (OpenMP_ON == 1 && Domain_Decomposition == 1) {
    if (! Use_Local_EnergyDiff_Algo)                                         
       Exit_Error (" Cannot use global energy calculation-based Metropolis algorithm with OpenMP - that is very inefficient. Switch-on Local Energy algo and restart.");
    C_MonteCarlo_Driver__omp (out_filename);
    return;
  }
#endif

  Initialize_MonteCarlo ();                                                   

  custom_swaps = custom_time = 0;

  if (!(MC_OutputFile=fopen(out_filename, "w"))) 
    Exit_Error ("Error in opening output file for Monte Carlo (from C).");    
  
  print_header ("Output from Monte Carlo (C-code)", MC_OutputFile);
  fprintf (MC_OutputFile, "\nTemperature (K), Average energy (eV/atom), Median energy (eV/atom) and STD (eV/atom), Specific Heat (eV/K/atom),  Std Cp (eV/K/atom)\n");

  if (Equilibrium_Steps < 0.0) EquilSteps = DEFAULT_EQUILIBRIUM_STEPS;
  else EquilSteps = Equilibrium_Steps;
  
  if (Use_Local_EnergyDiff_Algo)
    if (Perturbed_Sites_List_Allocated)  
      if (Loop_Parallelization)  fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code), local energy calculation and stored perturbed-sites ----------\n", Max_OpenMP_Threads);
      else                       fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop on single processor (C Code) with local energy and stored perturbed-sites ----------\n");  
    else                                 
      if (Loop_Parallelization)  fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code), local energy calculation and perturbed-sites on-the-fly ----------\n", Max_OpenMP_Threads);
      else                       fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop on single processor (C Code) with local energy and perturbed-sites on-the-fly ----------\n");
  else                                   
    if (Loop_Parallelization)    fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code) with global energy ----------\n", Max_OpenMP_Threads);
    else                         fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop on single processor (C Code) with global energy ----------\n");

  if (Use_Local_EnergyDiff_Algo)
    if (Perturbed_Sites_List_Allocated)  
      if (Loop_Parallelization)  printf ("\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code), local energy calculation and stored perturbed-sites ----------\n", Max_OpenMP_Threads);
      else                       printf ("\n----------Starting Monte-Carlo Loop on single processor (C Code) with local energy and stored perturbed-sites ----------\n");  
    else                                 
      if (Loop_Parallelization)  printf ("\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code), local energy calculation and perturbed-sites on-the-fly ----------\n", Max_OpenMP_Threads);
      else                       printf ("\n----------Starting Monte-Carlo Loop on single processor (C Code) with local energy and perturbed-sites on-the-fly ----------\n");
  else                                   
    if (Loop_Parallelization)    printf ("\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d cores (C Code) with global energy ----------\n", Max_OpenMP_Threads);
    else                         printf ("\n----------Starting Monte-Carlo Loop on single processor (C Code) with global energy ----------\n");
  fflush (MC_LogFile);
   
  for (i=0; i<nTKelvin_Sequence; i++) {
    fprintf (MC_LogFile, "\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);
    printf ("\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);  fflush (MC_LogFile);
     
    TkB = TKelvin_Sequence[i] * kBOLTZMANN;

    MonteCarlo (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, 0, 0, &EnergyStat);
    
    if (Equilibrium_Tolerance > 0.0 && EQUIL_STEPS_ACHIEVED == EquilSteps) {
      fprintf (MC_LogFile, "<<<<  WARNING :: Unable to reach equilibrium even after a million steps. >>>>\n");
      fprintf (MC_LogFile, "Re-check parameters, specially Energy Tolerance. Continuing simulation.\n"); fflush (MC_LogFile);
    }
         
    fprintf (MC_LogFile, "\n After Equilibriation at %f K, total energy (eV/atom): %lf\n", TKelvin_Sequence[i], EnergyStat.Mean);
    fprintf (MC_LogFile, "Energy tolerance achieved in %d MC steps: %f\n", EQUIL_STEPS_ACHIEVED, EQUIL_TOL_ACHIEVED);  fflush (MC_LogFile);
   
    MonteCarlo (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, Calculation_Steps, BLOCK_SIZE, &EnergyStat);

    fprintf (MC_LogFile, "Temperature, Average energy, Median energy and STD (eV/atom): %9.4f  %13.8lf  %13.8lf  %13.8lf  %e  %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat); 
    fprintf (MC_OutputFile, "%9.4f   %13.8lf   %13.8lf   %13.8lf   %e   %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat);  
    fflush (MC_LogFile); fflush (MC_OutputFile);

    Output_AtomCoordinates_poscar (TKelvin_Sequence[i]);
    Output_CheckPoint_File (TKelvin_Sequence[i]);
  }

#ifdef CTIME
  FILE *CUSTOM = fopen("custom_no_op","a");
  long double avg = custom_time/custom_swaps;
  fprintf(CUSTOM,"NCORES:%d NSWAPS:%Lf TOTAL_TIME:%Lf AVG_TIME:%Lf\n",1,custom_swaps, custom_time, avg);
  fclose(CUSTOM);
#endif

  print_header ("Successfully completed Monte Carlo simulation (from C.) !", MC_LogFile); 
  
  fclose (MC_LogFile);    
  fclose (MC_OutputFile);

  return;
}





#ifdef _OPENMP
void C_MonteCarlo_Driver__omp (char *out_filename)
{
  int i, EQUIL_STEPS_ACHIEVED, EquilSteps;
  struct Stat EnergyStat;
  float TkB, EQUIL_TOL_ACHIEVED;
  FILE *MC_OutputFile=NULL;
   
  if (!(MC_OutputFile=fopen(out_filename, "w"))) 
    Exit_Error ("Error in opening output file for Monte Carlo (from C).");     
  
  print_header ("Output from Monte Carlo (C-code)", MC_OutputFile);
  fprintf (MC_OutputFile, "\nTemperature (K), Average energy (eV/atom), Median energy (eV/atom) and STD (eV/atom), Specific Heat (eV/K/atom),  Std Cp (eV/K/atom)\n");

  if (Equilibrium_Steps < 0.0) EquilSteps = DEFAULT_EQUILIBRIUM_STEPS;
  else EquilSteps = Equilibrium_Steps;

  if (EquilSteps%NHISTORY < NHISTORY/2) EquilSteps=(EquilSteps/NHISTORY)*NHISTORY;
  else EquilSteps=(1+EquilSteps/NHISTORY)*NHISTORY;
 
  Initialize_MonteCarlo__omp ( );

  if (Perturbed_Sites_List_Allocated) fprintf (MC_LogFile, "  Sites perturbed are stored in list.\n");
  else                                fprintf (MC_LogFile, "  Perturbed sites to be determined on-the-fly in Metropolis algorithm.\n");  
  
  fprintf (MC_LogFile, "\n----------Starting parallel Monte-Carlo Loop on %d processors (C Code)----------\n", Total_nSuperBlocks);  fflush (MC_LogFile);
  if (Use_Local_EnergyDiff_Algo)
    if (Perturbed_Sites_List_Allocated)  printf ("\n----------Starting parallel Monte-Carlo Loop on %d processors (C Code) with local energy and stored perturbed-sites----------\n",Total_nSuperBlocks);  
    else                                 printf ("\n----------Starting parallel Monte-Carlo Loop on %d processors (C Code) with local energy and perturbed-sites determined on-fly----------\n",Total_nSuperBlocks);
  else                                   printf ("\n----------Starting parallel Monte-Carlo Loop on %d processors (C Code) with global energy----------\n",Total_nSuperBlocks);
 
  for (i=0; i<nTKelvin_Sequence; i++) {
    fprintf (MC_LogFile, "\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);
    printf ("\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);  fflush (MC_LogFile);
     
    TkB = TKelvin_Sequence[i] * kBOLTZMANN;

    if (Perturbed_Sites_List_Allocated) MonteCarlo_LocalEnergy_wMemory__omp (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, 0, 0, &EnergyStat);
    else                                MonteCarlo_LocalEnergy__omp (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, 0, 0, &EnergyStat);

    if (Equilibrium_Tolerance > 0.0 && EQUIL_STEPS_ACHIEVED == EquilSteps) {
      fprintf (MC_LogFile, "<<<<  WARNING :: Unable to reach equilibrium even after a million steps. >>>>\n");
      fprintf (MC_LogFile, "Re-check parameters, specially Energy Tolerance. Continuing simulation.\n"); fflush (MC_LogFile);
    }
         
    fprintf (MC_LogFile, "\n After Equilibriation at %f K, total energy (eV/atom): %lf\n", TKelvin_Sequence[i], EnergyStat.Mean);
    fprintf (MC_LogFile, " Energy tolerance achieved in %d MC steps: %f\n", EQUIL_STEPS_ACHIEVED, EQUIL_TOL_ACHIEVED);  fflush (MC_LogFile);
   
    if (Perturbed_Sites_List_Allocated) MonteCarlo_LocalEnergy_wMemory__omp (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, Calculation_Steps, BLOCK_SIZE, &EnergyStat);
    else                                MonteCarlo_LocalEnergy__omp (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, Calculation_Steps, BLOCK_SIZE, &EnergyStat);

    fprintf (MC_LogFile, "Temperature, Average energy, Median energy and STD (eV/atom): %9.4f  %13.8lf  %13.8lf  %13.8lf  %e  %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat); 
    fprintf (MC_OutputFile, "%9.4f   %13.8lf   %13.8lf   %13.8lf   %e   %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat);  
    fflush (MC_LogFile); fflush (MC_OutputFile);

    Output_AtomCoordinates_poscar (TKelvin_Sequence[i]);
    Output_CheckPoint_File (TKelvin_Sequence[i]);
  }  
 
  print_header ("Successfully completed Monte Carlo simulation (from C) !", MC_LogFile);    
  fclose (MC_OutputFile); 
  fclose (MC_LogFile);  
  
  return;
}
#endif


void Free_Global_Memory_inC (void)
{
  int i, j, k;

  
#ifdef _OPENMP  
  if (OpenMP_ON) Free_Global_Memory__omp ();
#else
  Free_Global_Memory ();
#endif

  for (i=0; i<Base_Lattice.nSites; i++) free (Base_Lattice.Sites[i].Label); 
  free (Base_Lattice.Sites);
  free (Base_Lattice.Labels);
    
  for (i=0; i<nElements; i++) {
    free (ElementList[i].Name); 
    free (ElementList[i].Theta);
  }
  free (ElementList);
  if (nFrozen_Elements > 0) free (Frozen_Elements);  
  free (TKelvin_Sequence);
  
  for (i=0; i<nClusters; i++) {
   
    for (j=0; j<ClusterList[i].nChildGroups; j++) {
      for (k=0; k<ClusterList[i].ChildGroup_List[j].nChild; k++) free(ClusterList[i].ChildGroup_List[j].Child[k].R);
      free (ClusterList[i].ChildGroup_List[j].Child);
    }
    free (ClusterList[i].ChildGroup_List);
   
    for (j=0; j<ClusterList[i].nECIs; j++) free(ClusterList[i].ECI_List[j].Index);
    free (ClusterList[i].ECI_List);
   
    free(ClusterList[i].ClusterID);
  }
  free (ClusterList);
  free (ECI_Vec);
  
  for (i=0; i<MCSystem.nSites; i++) 
    if (MCSystem.Sites[i].Active) {
      free (MCSystem.Sites[i].Site_Neighbors);
      if (Perturbed_Sites_List_Allocated) free(MCSystem.Sites[i].Perturbed_Sites);
    }
  free (MCSystem.Sites);
  free (MCSystem.nAtoms);

  return;
}


void Output_AtomCoordinates_poscar (float TempK)
{
  FILE *fptmp=NULL;
  char *filename=NULL, *tmpname=NULL;
  int i, j;
  int **atom_index=NULL, *natom=NULL;
  struct Coord b;
  struct AtomLabeled_Coord *a=NULL;


  filename=alloc_1Dchar_array(40); tmpname=alloc_1Dchar_array(20);
  strcpy (filename, "Config_"); 
  snprintf (tmpname, 20, "%5.1f", TempK); strcat (tmpname, ".poscar");
  strcat (filename, tmpname);

  if (!(fptmp=fopen(filename,"w"))) {              
    fprintf (MC_LogFile, "Error in Output_AtomCoordinates_poscar: Unable to open %s.\nExit \n", filename); 
    exit (EXIT_FAILURE);
  }

  fprintf (fptmp, "Configuration_at_T_%fK\n", TempK);
  fprintf (fptmp, "1.0\n");
  for (i=0; i<3; i++) 
    fprintf (fptmp, "  %14.10f %14.10f %14.10f\n", MCSystem.Lattice_Vectors[i].x, MCSystem.Lattice_Vectors[i].y, MCSystem.Lattice_Vectors[i].z);

  atom_index=alloc_2Dint_array(nElements, MCSystem.nSites);
  natom=alloc_1Dint_array (nElements); for (i=0; i<nElements; i++) natom[i]=0;

  a=MCSystem.Sites;
  for (i=0; i<MCSystem.nSites; i++) {
    atom_index[a[i].iAtom][natom[a[i].iAtom]]=i;    
    natom[a[i].iAtom]++;
  }  

  for (i=0; i<nElements; i++) fprintf (fptmp, "%s  ", ElementList[i].Name); fprintf (fptmp, "\n");
  for (i=0; i<nElements; i++) fprintf (fptmp, "%d  ", natom[i]); fprintf (fptmp, "\n");
  fprintf (fptmp, "Direct\n");

  for (i=0; i<nElements; i++) {
    for (j=0; j<natom[i]; j++) {
      b=Cart2Frac_coordinates (&MCSystem.Sites[atom_index[i][j]].R, &MCSystem.Ncell);
      fprintf (fptmp, "%14.10f %14.10f %14.10f\n", b.x, b.y, b.z);
    }
  }  
  
  fprintf (MC_LogFile, "\n  --->Coordinates at %f written to file %s.\n", TempK, filename);
  
  free (filename); free (tmpname);
  fclose (fptmp);

  return;
} 



void Output_CheckPoint_File (float TempK)
{
  FILE *fptmp=NULL;
  char *filename=NULL;
  int i;
  struct AtomLabeled_Coord *a=NULL;


  filename=alloc_1Dchar_array(40);
  strcpy (filename, "CheckPoint.mc"); 

  if (!(fptmp=fopen(filename,"w"))) {              
    fprintf (MC_LogFile, "Error in Output_AtomCoordinates_poscar: unable to open file %s.\n", filename); 
    exit (EXIT_FAILURE);
  }

  fprintf (fptmp, "Temperature = %f\n", TempK);
  fprintf (fptmp, "Lattice_Vectors =\n");
  for (i=0; i<3; i++) 
    fprintf (fptmp, "  %14.10f %14.10f %14.10f\n", MCSystem.Lattice_Vectors[i].x, MCSystem.Lattice_Vectors[i].y, MCSystem.Lattice_Vectors[i].z);
  
  fprintf (fptmp, "nCell = %d %d %d\n", (int)(MCSystem.Ncell.x), (int)(MCSystem.Ncell.y), (int)(MCSystem.Ncell.z));
  fprintf (fptmp, "Element_List = "); for (i=0; i<nElements; i++) fprintf (fptmp, "%s ", ElementList[i].Name); fprintf (fptmp, "\n");
  fprintf (fptmp, "Frozen_Elements = "); for (i=0; i<nFrozen_Elements; i++) fprintf (fptmp, "%s ", ElementList[Frozen_Elements[i]].Name); fprintf (fptmp, "\n");
  fprintf (fptmp, "Atom_Numbers = %d\n", MCSystem.nSites); 
  fprintf (fptmp, "Coordinates = \n");
  a=MCSystem.Sites;
  for (i=0; i<MCSystem.nSites; i++) {
    fprintf (fptmp, "%14.10f %14.10f %14.10f %s %d\n", a[i].R.x, a[i].R.y, a[i].R.z, ElementList[a[i].iAtom].Name, a[i].iBaseLat);
  }
 
  fprintf (MC_LogFile, "\n  --->CheckPoint file created for T = %fK.\n", TempK);
  
  free (filename);
  fclose (fptmp);

  return;
} 



struct Coord Cart2Frac_coordinates (struct Coord *in, struct Coord *ncell) 
{
  struct Coord p;
  p.x=in->x/ncell->x; p.y=in->y/ncell->y; p.z=in->z/ncell->z;   
  return p;
}



