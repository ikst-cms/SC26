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

void Output_AtomCoordinates_poscar (float);
void Output_CheckPoint_File (float);
struct Coord Cart2Frac_coordinates (struct Coord *, struct Coord *);


void Read_Elements_From_File()
{
  if (!(MPI_Element_Rep_File = fopen(MPI_ELEMENT_REP_FILE, "r"))) {
    Exit_Error ("Error in opening element rep file for Monte Carlo (from C).");
  }

  char log_filename[10];
  fscanf(MPI_Element_Rep_File, "%s", log_filename);

  int i,j;
  clock_t tt;
  
  tt = clock();

  fscanf(MPI_Element_Rep_File, "%d", &nElements);
    
  if (!(ElementList=(struct Element *) malloc (nElements*sizeof(struct Element)))) 
    Exit_Error ("Unable to allocate ElementList.");

  for (i=0; i<nElements; i++) ElementList[i].Name  = alloc_1Dchar_array(2); 
  for (i=0; i<nElements; i++) ElementList[i].Theta = alloc_1Ddouble_array(nElements);

  
  // printing element names
  for (i=0; i<nElements; i++) {
  	fscanf(MPI_Element_Rep_File, "%s ", ElementList[i].Name);
  }

  // printing Theta Rep
  for (i=0; i<nElements; i++) {
  	for (j=0; j<nElements; j++) {
      fscanf(MPI_Element_Rep_File, "%lf ", &ElementList[i].Theta[j]);
    }
  }

  tt = clock() - tt;

  if(grank==0){
    if (!(MC_LogFile=fopen(log_filename, "w"))) {
      Exit_Error ("Error in opening log file for Monte Carlo (from C).");
    }
  }

  tt = clock() - tt;
     
  if(grank==0){
    print_header ("Starting Monte Carlo in C", MC_LogFile); 
    fprintf (MC_LogFile, "\nSuccessfully received elements from File ...: %.10lf Secs\n",((double)tt)/CLOCKS_PER_SEC);
  }
  fclose(MPI_Element_Rep_File);
}


void Read_BaseLattice_From_File() 
{
  if (!(MPI_Base_Lattice_File = fopen(MPI_BASE_LATTICE_FILE, "r"))) {
    Exit_Error ("Error in opening base lattice file for Monte Carlo (from C).");
  }
  int i;
  clock_t tt;
  tt = clock();

  fscanf(MPI_Base_Lattice_File, "%d", &Base_Lattice.nLabels);
  Base_Lattice.Labels = alloc_1Dchar_array (Base_Lattice.nLabels+1);

  fscanf(MPI_Base_Lattice_File, "%s", Base_Lattice.Labels);

  for (i = 0; i < 3; i++) {
    fscanf(MPI_Base_Lattice_File, "%f %f %f", &Base_Lattice.Lattice_Vectors[i].x, 
            &Base_Lattice.Lattice_Vectors[i].y, &Base_Lattice.Lattice_Vectors[i].z);
  }

  fscanf(MPI_Base_Lattice_File, "%d\n", &Base_Lattice.nSites);

  Base_Lattice.Sites = (struct Labeled_Coord *) malloc (Base_Lattice.nSites*sizeof(struct Labeled_Coord));
  for (i=0; i<Base_Lattice.nSites; i++) Base_Lattice.Sites[i].Label = alloc_1Dchar_array (2);

  for (i=0; i<Base_Lattice.nSites; i++) {
    fscanf(MPI_Base_Lattice_File, "%s %f %f %f", Base_Lattice.Sites[i].Label, 
            &Base_Lattice.Sites[i].R.x, &Base_Lattice.Sites[i].R.y, &Base_Lattice.Sites[i].R.z);
  }

  if(grank==0){
    tt = clock() - tt;
    fprintf (MC_LogFile, "\nSuccessfully received Base Lattice info from File... : %.10lf Secs\n", ((double)tt)/CLOCKS_PER_SEC);
  }
  fclose(MPI_Base_Lattice_File);
}


void Read_SuperLattice_From_File() 
{
  if (!(MPI_Super_Lattice_File = fopen(MPI_SUPER_LATTICE_FILE, "r"))) {
    Exit_Error ("Error in opening base lattice file for Monte Carlo (from C).");
  }
  int i,j;
  clock_t tt;
  tt = clock();

  fscanf(MPI_Super_Lattice_File, "%d", &nFrozen_Elements);
  if (nFrozen_Elements > 0) {
    Frozen_Elements = alloc_1Dint_array (nFrozen_Elements); 
    for (i=0; i<nFrozen_Elements; i++) {
      char frozen_tmp[100];
      fscanf(MPI_Super_Lattice_File,"%s",frozen_tmp);
      for (j=0; j<nElements; j++) if (strcmp (frozen_tmp, ElementList[j].Name) == 0) Frozen_Elements[i] = j;
    }
  }

  fscanf(MPI_Super_Lattice_File, "%f %f %f", &MCSystem.Ncell.x, &MCSystem.Ncell.y, &MCSystem.Ncell.z);

  for (i=0; i<3; i++) {
    fscanf(MPI_Super_Lattice_File, "%f %f %f", &MCSystem.Lattice_Vectors[i].x,
            &MCSystem.Lattice_Vectors[i].y, &MCSystem.Lattice_Vectors[i].z);
  }

  MCSystem.nAtoms = alloc_1Dint_array (nElements); 
  for (i=0; i<nElements; i++) MCSystem.nAtoms[i]=0;

  fscanf(MPI_Super_Lattice_File, "%d", &MCSystem.nSites);
  
  MCSystem.Sites = (struct AtomLabeled_Coord *) malloc (MCSystem.nSites*sizeof(struct AtomLabeled_Coord));
  MCSystem.nActiveSites=0;

  for (i=0; i<MCSystem.nSites; i++) {
    char atom[100];
    fscanf(MPI_Super_Lattice_File,"%f %f %f %d %s", &MCSystem.Sites[i].R.x,
            &MCSystem.Sites[i].R.y, &MCSystem.Sites[i].R.z, &MCSystem.Sites[i].iBaseLat,
            atom);
    
    for (j=0; j<nElements; j++) 
      if (strcmp (atom, ElementList[j].Name) == 0) {
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

#ifdef OLD_HASH_KEY  
  Generate_Old_HashTable (&MCSystem, &Base_Lattice);      
#else
  Generate_New_HashTable (&MCSystem);                                   
#endif

  if(grank==0){
    tt = clock() - tt;

    fprintf (MC_LogFile, "\nFound %d active sites out of %d sites in the system.\n", MCSystem.nActiveSites, MCSystem.nSites);
    fprintf (MC_LogFile, "Successfully received SuperLattice info from File; generated hash-table for coordinates...: %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
  }
  
  fclose(MPI_Super_Lattice_File);
}


void Read_MC_Parameters ()
{
  if (!(MPI_MC_Parameters_File = fopen(MPI_MC_PARAMETERS_FILE, "r"))) {
    Exit_Error ("Error in opening mc param file for Monte Carlo (from C).");
  }
  int i;
  clock_t tt;
  tt = clock();
  fscanf(MPI_MC_Parameters_File, "%d\n", &nTKelvin_Sequence);
  TKelvin_Sequence = alloc_1Dfloat_array (nTKelvin_Sequence);
  for (i=0; i<nTKelvin_Sequence; i++){
    fscanf(MPI_MC_Parameters_File,"%f ", &TKelvin_Sequence[i]);
  }

  fscanf(MPI_MC_Parameters_File, "%d\n", &Equilibrium_Steps);
  fscanf(MPI_MC_Parameters_File, "%f\n", &Equilibrium_Tolerance);
  fscanf(MPI_MC_Parameters_File, "%d\n", &Calculation_Steps);
  fscanf(MPI_MC_Parameters_File, "%d\n", &NHISTORY);
  fscanf(MPI_MC_Parameters_File, "%d\n", &BLOCK_SIZE);
  fscanf(MPI_MC_Parameters_File, "%d\n", &minEQUIL_STEPS);
  
  if(grank==0){
    tt = clock() - tt;
    fprintf (MC_LogFile, "\nSuccessfully received Monte-Carlo parameters from File...: %.10f Secs \n",((double)tt)/CLOCKS_PER_SEC);
  }

  fclose(MPI_MC_Parameters_File);
}


void Read_Cluster(FILE* Input_File, struct Cluster *curCluster, int Nk) {
  int i;
  for (i=0; i<Nk; i++) {
    fscanf(Input_File, "%f %f %f", &curCluster->R[i].x, 
            &curCluster->R[i].y, &curCluster->R[i].z);
  }
}


void Read_ClusterInfo_From_File () 
{
  int i, j, k, cls;
  clock_t tt;

  if (!(MPI_Cluster_Info_File = fopen(MPI_CLUSTER_INFO_FILE, "r"))) {
    Exit_Error ("Error in opening mc param file for Monte Carlo (from C).");
  }

  tt = clock();

  fscanf(MPI_Cluster_Info_File, "%d", &nClusters);
                                                                   
  if (!(ClusterList = (struct ClusterFamily *) malloc (sizeof(struct ClusterFamily)*(nClusters)))) {
    Exit_Error ("Unable to allocate memory for Cluster List.");
  }

  nCE_TERMS=1;
  for(cls=0; cls<nClusters; cls++) {
    ClusterList[cls].ClusterID = alloc_1Dchar_array (15);
    fscanf(MPI_Cluster_Info_File, "%d", &ClusterList[cls].Ksites);
    fscanf(MPI_Cluster_Info_File, "%s", ClusterList[cls].ClusterID);
    fscanf(MPI_Cluster_Info_File, "%d", &ClusterList[cls].nECIs);                          
    ClusterList[cls].ECI_List = (struct Indexed_ECI *) malloc (sizeof(struct Indexed_ECI)*ClusterList[cls].nECIs);
    for(j=0; j<ClusterList[cls].nECIs; j++) {
      fscanf(MPI_Cluster_Info_File, "%lf", &ClusterList[cls].ECI_List[j].ECI);
      ClusterList[cls].ECI_List[j].Index = alloc_1Dint_array (ClusterList[cls].Ksites); 
      for (k=0; k<ClusterList[cls].Ksites; k++) {
        fscanf(MPI_Cluster_Info_File, "%d", &ClusterList[cls].ECI_List[j].Index[k]);
      }
    }
    nCE_TERMS += ClusterList[cls].nECIs;
    ClusterList[cls].Parent.R = alloc_1DCoord_array (ClusterList[cls].Ksites);
    Read_Cluster(MPI_Cluster_Info_File, &ClusterList[cls].Parent, ClusterList[cls].Ksites);

    ClusterList[cls].nChildGroups = Base_Lattice.nSites;
    ClusterList[cls].ChildGroup_List = (struct ChildGroup *) malloc(sizeof(struct ChildGroup)*ClusterList[cls].nChildGroups);
    for (j=0; j<ClusterList[cls].nChildGroups; j++) {
      fscanf(MPI_Cluster_Info_File, "%d", &ClusterList[cls].ChildGroup_List[j].nChild);
      fscanf(MPI_Cluster_Info_File, "%d", &ClusterList[cls].ChildGroup_List[j].BL_Index);
      
      ClusterList[cls].ChildGroup_List[j].Child = (struct Cluster *) malloc(sizeof(struct Cluster)*ClusterList[cls].ChildGroup_List[j].nChild);
      for (k=0; k<ClusterList[cls].ChildGroup_List[j].nChild; k++) {
        ClusterList[cls].ChildGroup_List[j].Child[k].R = alloc_1DCoord_array (ClusterList[cls].Ksites);
        Read_Cluster (MPI_Cluster_Info_File, &ClusterList[cls].ChildGroup_List[j].Child[k], ClusterList[cls].Ksites);
      }     
    } 
  }
  
  if(grank==0){
    tt = clock() - tt;
    fprintf (MC_LogFile, "\nSuccessfully received Cluster and ECI information from File... : %.10lf Secs \n",((double)tt)/CLOCKS_PER_SEC);
  }

  ECI_Vec = alloc_1Ddouble_array (nCE_TERMS);
  fscanf(MPI_Cluster_Info_File, "%lf", &ECI_Vec[0]);
  for (i=1, cls=0; cls<nClusters; cls++) 
    for (j=0; j<ClusterList[cls].nECIs; j++, i++) ECI_Vec[i] = ClusterList[cls].ECI_List[j].ECI;
 
  for (cls=0, LargestKsites=1; cls<nClusters; cls++) 
    if (ClusterList[cls].Ksites > LargestKsites) LargestKsites = ClusterList[cls].Ksites;

  for (cls=0, maxClusterComponents=1; cls<nClusters; cls++) 
    if (ClusterList[cls].nECIs > maxClusterComponents) maxClusterComponents = ClusterList[cls].nECIs;
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
  
  if(grank==0){
    printf ("\n  Average number of neighbors (from C-code): %d\n", (int)(avgNNB));
  
    fprintf (MC_LogFile, "\nSuccessfully created nearest-neighbor list for all active coordinates... : %.10lf Secs \n", ((double)tt)/CLOCKS_PER_SEC);
    fprintf (MC_LogFile, "Average number of neighbors (from C-code): %d\n", (int)(avgNNB));
  }
  maxAtom_Exchng_Dist = max_minDist;           

  return 1; 
}


void Create_Perturbed_SiteList_inC ()
{
  int i, j, p;
  struct Coord Rdiff;
  float pnorm=0.0, d;
  int nMax_PerturbedSites, *tmparray=NULL;
  clock_t tt;
  tt = clock();

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
    
  
  if(grank==0){  
    fprintf (MC_LogFile, "\n  Perturbation length uesd in calculating energy difference: %4.2f\n", Perturbation_Distance);  fflush(MC_LogFile);
    printf ("\n  Perturbation length uesd in calculating energy difference: %4.2f\n", Perturbation_Distance);

    fprintf (MC_LogFile, "  Attempting to find number of sites in perturbed region. If more than 60%% of total sites, global energy calculation to be used for Metropolis step.\n"); fflush(MC_LogFile);
    printf ("  Attempting to find number of sites in perturbed region. If more than 60%% of total sites, global energy calculation to be used for Metropolis step.\n");
  }

  nMax_PerturbedSites=0;  
  for (i=0; i<Base_Lattice.nSites; i++) {             
    p=0;
    if (MCSystem.Sites[i].Active) {
      for (j=0; j<MCSystem.nSites; j++) 
        if (Periodic_Distance (&MCSystem.Sites[i].R, &MCSystem.Sites[j].R, &MCSystem.Ncell) < Perturbation_Distance) p++;
    }
    if (p > nMax_PerturbedSites) nMax_PerturbedSites=p;
  }

  Use_Local_EnergyDiff_Algo=1;

  if(grank==0){
    fprintf (MC_LogFile, "  Approximate Number of sites within perturbed region based on largest cluster in CE (compare it with 60%% of total sites): %d (%d)\n", nMax_PerturbedSites, (int) (0.6*MCSystem.nSites)); fflush(MC_LogFile);
    printf ("  Approximate Number of sites within perturbed region based on largest cluster in CE (compare it with 60%% of total sites): %d (%d)\n", nMax_PerturbedSites, (int) (0.6*MCSystem.nSites));
  }

  tmparray = alloc_1Dint_array (nMax_PerturbedSites);             
  max_psites = nMax_PerturbedSites;
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

  if(grank==0){
    fprintf (MC_LogFile, "\n   Using Local calculation to find energy difference on atom-exchange (with perturbed sites stored in memory)... : %f Secs \n", ((double)tt)/CLOCKS_PER_SEC);  
    fflush(MC_LogFile);
    printf ("\n   Using Local calculations to find energy difference on atom-exchange (with perturbed sites stored in memory)\n");
  }

  Perturbed_Sites_List_Allocated =1;
}


void C_MonteCarlo_Driver (char *out_filename)
{
  int i, EQUIL_STEPS_ACHIEVED, EquilSteps;
  struct Stat EnergyStat;
  float TkB, EQUIL_TOL_ACHIEVED;
  FILE *MC_OutputFile=NULL;

  Initialize_MonteCarlo (); 

  if(grank==0){
    custom_swaps = custom_time = 0;
    if (!(MC_OutputFile=fopen(out_filename, "w"))) 
      Exit_Error ("Error in opening output file for Monte Carlo (from C).");   
  } 
  
  if(grank==0){
    print_header ("Output from Monte Carlo (C-code)", MC_OutputFile);
    fprintf (MC_OutputFile, "\nTemperature (K), Average energy (eV/atom), Median energy (eV/atom) and STD (eV/atom), Specific Heat (eV/K/atom),  Std Cp (eV/K/atom)\n");
  }

  if (Equilibrium_Steps < 0.0) EquilSteps = DEFAULT_EQUILIBRIUM_STEPS;
  else EquilSteps = Equilibrium_Steps;

  if(grank==0){
    fprintf (MC_LogFile, "\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d processes (C Code), local energy calculation and stored perturbed-sites  ----------\n", gsize);
    printf ("\n----------Starting Monte-Carlo Loop with Loop Parallelization on %d processes (C Code), local energy calculation and stored perturbed-sites  ----------\n", gsize);
    fflush (MC_LogFile);
  }
  
  for (i=0; i<nTKelvin_Sequence; i++) {
    if(grank==0){
      fprintf (MC_LogFile, "\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);
      printf ("\nPerforming simulation for T = %f K\n", TKelvin_Sequence[i]);  fflush (MC_LogFile);
    }
     
    TkB = TKelvin_Sequence[i] * kBOLTZMANN;

    MonteCarlo (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, 0, 0, &EnergyStat);
    
    if (Equilibrium_Tolerance > 0.0 && EQUIL_STEPS_ACHIEVED == EquilSteps) {
      if(grank==0){
        fprintf (MC_LogFile, "<<<<  WARNING :: Unable to reach equilibrium even after a million steps. >>>>\n");
        fprintf (MC_LogFile, "Re-check parameters, specially Energy Tolerance. Continuing simulation.\n"); fflush (MC_LogFile);
      }
    }

    if(grank==0){
      fprintf (MC_LogFile, "\n After Equilibriation at %f K, total energy (eV/atom): %lf\n", TKelvin_Sequence[i], EnergyStat.Mean);
      fprintf (MC_LogFile, "Energy tolerance achieved in %d MC steps: %f\n", EQUIL_STEPS_ACHIEVED, EQUIL_TOL_ACHIEVED);  fflush (MC_LogFile);
    }

    MonteCarlo (TkB, NHISTORY, minEQUIL_STEPS, Equilibrium_Tolerance, EquilSteps, &EQUIL_TOL_ACHIEVED, &EQUIL_STEPS_ACHIEVED, Calculation_Steps, BLOCK_SIZE, &EnergyStat);

    if(grank==0){
      fprintf (MC_LogFile, "Temperature, Average energy, Median energy and STD (eV/atom): %9.4f  %13.8lf  %13.8lf  %13.8lf  %e  %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat); 
      fprintf (MC_OutputFile, "%9.4f   %13.8lf   %13.8lf   %13.8lf   %e   %e\n", TKelvin_Sequence[i], EnergyStat.Mean, EnergyStat.Median, EnergyStat.StdDev, EnergyStat.SpHeat, EnergyStat.dSpHeat);  
      fflush (MC_LogFile); fflush (MC_OutputFile);

      Output_AtomCoordinates_poscar (TKelvin_Sequence[i]);
      Output_CheckPoint_File (TKelvin_Sequence[i]);
    }
  }  

  if(grank==0){

#ifdef CTIME

#ifdef Hybrid
    FILE *CUSTOM = fopen("custom_mpi_hybrid","a");
    long double avg = custom_time/custom_swaps;
    fprintf(CUSTOM,"NCORES:%d NSWAPS:%Lf TOTAL_TIME:%Lf AVG_TIME:%Lf\n",gsize,custom_swaps, custom_time, avg);
    fclose(CUSTOM);
#else
    FILE *CUSTOM = fopen("custom_no_op","a");
    long double avg = custom_time/custom_swaps;
    fprintf(CUSTOM,"NCORES:%d NSWAPS:%Lf TOTAL_TIME:%Lf AVG_TIME:%Lf\n",gsize,custom_swaps, custom_time, avg);
    fclose(CUSTOM);
#endif

#endif

    print_header ("Successfully completed Monte Carlo simulation (from C.) !", MC_LogFile); 

    fclose (MC_LogFile);    
    fclose (MC_OutputFile);
  }

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

void Read_Exchange_Dist_and_Output_File ()
{
  if (!(MPI_Exchange_Dist_Output_File = fopen(MPI_EXCHANGE_DIST_OUTPUT_FILE, "r"))) {
    Exit_Error ("Error in opening exchange dist output file for Monte Carlo (from C).");
  }
  float ext_minDist;
  char output_file[100];

#ifdef Hybrid
  fscanf(MPI_Exchange_Dist_Output_File, "%f %s %d %d", &ext_minDist, output_file, &psp_maxList, &psp_maxCap);
#else
  fscanf(MPI_Exchange_Dist_Output_File, "%f %s", &ext_minDist, output_file);
#endif

  int check = Create_Neighbour_List_inC(ext_minDist);
  if(check==0) {
    Exit_Error("Sites found with empty neighbour-list and occupied by atom not in frozen-element list.");
  }
  fclose(MPI_Exchange_Dist_Output_File);
  Create_Perturbed_SiteList_inC();

#ifdef Hybrid
  psp_labels = alloc_1Dint_array (MCSystem.nSites);
  psp_curCap = alloc_1Dint_array (psp_maxList);
  psp_list_jump = alloc_2Dint_array(psp_maxList, psp_maxCap);
  psp_list_neigh = alloc_2Dint_array(psp_maxList, psp_maxCap);
  psp_act_jump = alloc_1Dint_array(psp_maxCap);
  psp_act_neigh = alloc_1Dint_array(psp_maxCap);
  psp_local_energy = alloc_1Ddouble_array(psp_maxCap);
  psp_final_energy = alloc_1Ddouble_array(psp_maxCap);

#ifdef OPT_PSITE
  psp_max_psites = (psp_maxCap*2*max_psites);
  psp_act_psites = alloc_1Dint_array(psp_max_psites);
  psp_swap_ind = alloc_1Dint_array(psp_max_psites);  
#endif

#else
  mpi_total_psites = alloc_1Dint_array(max_psites*2);
#endif

  pcg32_srandom_r(&rng,42u,54u);
  pcg32_srandom_r(&rng_select,42u,54u);

  C_MonteCarlo_Driver(output_file);
}

void preprocessClusters(){
  int i,j,ind;
  for(ind=0;ind<MCSystem.nSites;ind++){
    MCSystem.Sites[ind].correspondingCluster = alloc_1Dint_array(nClusters);
    for (i=0; i<nClusters; i++) {
      MCSystem.Sites[ind].correspondingCluster[i] = -1; 
      for (j=0; j<ClusterList[i].nChildGroups; j++) {                             
          if (ClusterList[i].ChildGroup_List[j].BL_Index 
                == MCSystem.Sites[ind].iBaseLat 
                  && ClusterList[i].ChildGroup_List[j].nChild > 0) { 
            MCSystem.Sites[ind].correspondingCluster[i] = j;
            break;           
          }
      }
    }
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &gsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &grank);

#ifdef CLIST
  if(grank==0){
    FILE *CUSTOM = fopen("custom_list_size","w");
    fclose(CUSTOM);
  }
#endif

  if(grank==0){
	  printf("\n\t...Reading data from file to perform Monte Carlo Simulation...\n\n");
  }
	Read_Elements_From_File();
  Read_BaseLattice_From_File();
  Read_SuperLattice_From_File();
  Read_MC_Parameters();
  Read_ClusterInfo_From_File();
#ifdef OPT_PRE_CLUSTER
  preprocessClusters();
#endif
  Read_Exchange_Dist_and_Output_File();
  MPI_Finalize();
	return 0;
}