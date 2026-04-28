//-----------------------------------------------------------------------------------------------------+
// Metropolis-based atom-exchange parallel Monte Carlo code (for shared memory only). The code uses    |
// Domain-decomposition technique wherein the spatial system is divided into blocks and each block is  |
// assigned to a core/thread for computation. The code is not "truly" parallel since detailed balance  |
// condition must be satisfied (and Markov chain is inherently sequential). But all core computation   |
// performed in parallel to give better runtime. This code is interfaced to Python code PyCMonC.       |
//                                                                                                     |
// NOTE:: Comparison of energy calculation between C-code and Python-code shows a difference in 5-th   |
//        decimal position which could be due to (double) precision used in C-code for energy          |
//                                                                                                     |
//-----------------------------------------------------------------------------------------------------+

 

#define ENE_CUTOFF__OMP 1.e-10



void Initialize_MonteCarlo__omp (void);
void Get_Block_parameters (void);
void Do_Domain_Decomposition (void);
void Create_ActiveSites_Domains (void);
void Get_Block_Coordinate_for_Site__omp (struct Coord *, struct IntCoord *, struct IntCoord *, struct IntCoord *, struct IntCoord *, int *, int *);
int Number_of_Metropolis_Steps_perBlock (void);
void Free_Global_Memory__omp (void);
void Seed_RNG__omp (void);


double Configuration_Energy__omp (struct SuperCell *, struct Element *);
double Average_ClusterEnergy__omp (struct SuperCell *, struct ClusterFamily *, struct Element *);
void Cluster_Occupations__omp (struct Coord *, struct Cluster *, int, struct AtomLabeled_Coord *, struct Coord *, int *);

double Average_perSite_dE__omp (struct AtomLabeled_Coord *, struct SuperCell *, struct Element *, int, int, int, int); 
unsigned int Site_Perturbed__omp (struct Coord *, struct Coord *, struct Coord *, float *);
void Xchngd_Cluster_Occupations__omp (struct Coord *, struct Cluster *, int, struct AtomLabeled_Coord *, struct Coord *, int *, int *, int, int, int, int);

double Local_dE__omp (struct SuperCell *, struct Element *, int, int, float *);
void MonteCarlo_LocalEnergy__omp (float, int, int, float, int, float *, int *, int, int, struct Stat *);

double Local_dE_wMemory__omp (struct SuperCell *, struct Element *, int, int, int *, int);
void MonteCarlo_LocalEnergy_wMemory__omp (float, int, int, float, int, float *, int *, int, int, struct Stat *);



pcg32_random_t *rng_parallel=NULL;  


struct SubBlock {
  int *Sites;                         
  int nSites;     
};

struct SuperBlock {
  struct SubBlock *subDomain;
};

struct SuperBlock *Domain=NULL;             
struct SuperBlock *aSiteDomain=NULL;        

int nMetropolis_per_SuperBlock;




void Initialize_MonteCarlo__omp (void)                     
{
  
  Get_Block_parameters ();
  Do_Domain_Decomposition ();
  Create_ActiveSites_Domains ();
  nMetropolis_per_SuperBlock = Number_of_Metropolis_Steps_perBlock ();

  Seed_RNG__omp ();   
  return;
}


void Get_Block_parameters (void)
{
  Total_nSuperBlocks = nSuperBlocks.x * nSuperBlocks.y * nSuperBlocks.z;     
  
  nSubBlocks.x = (nSuperBlocks.x > 1) ? 2 : 1;
  nSubBlocks.y = (nSuperBlocks.y > 1) ? 2 : 1;
  nSubBlocks.z = (nSuperBlocks.z > 1) ? 2 : 1;
  Total_nSubBlocks = nSubBlocks.x * nSubBlocks.y * nSubBlocks.z;             
   
  dSuperBlock.x=((int)(MCSystem.Ncell.x)/nSuperBlocks.x); if (nSuperBlocks.x>1) dSubBlock.x=dSuperBlock.x/2; else dSubBlock.x=dSuperBlock.x;
  dSuperBlock.y=((int)(MCSystem.Ncell.y)/nSuperBlocks.y); if (nSuperBlocks.y>1) dSubBlock.y=dSuperBlock.y/2; else dSubBlock.y=dSuperBlock.y;
  dSuperBlock.z=((int)(MCSystem.Ncell.z)/nSuperBlocks.z); if (nSuperBlocks.z>1) dSubBlock.z=dSuperBlock.z/2; else dSubBlock.z=dSuperBlock.z;

  if (dSuperBlock.x * dSuperBlock.y * dSuperBlock.z * Total_nSuperBlocks != ((int)(MCSystem.Ncell.x * MCSystem.Ncell.y * MCSystem.Ncell.z)))
    Exit_Error ("Unable to perform uniform domain-decomposition. Check the Supercell-size and SuperBlock size inputs. Exit.");

  if (dSubBlock.x * dSubBlock.y * dSubBlock.z * Total_nSubBlocks != (dSuperBlock.x * dSuperBlock.y *dSuperBlock.z))
    Exit_Error ("Unable to perform uniform domain-decomposition. Check the Supercell-size and SuperBlock size inputs (super-block extension should be a multiple of 2). Exit.");

  return;
}


void Do_Domain_Decomposition (void)
{
  int i, j, k, iSup, iSub, Tot;
  struct SuperBlock *tmpDomain=NULL;
 
  tmpDomain = (struct SuperBlock *) malloc (sizeof(struct SuperBlock)*Total_nSuperBlocks); 
  k=(int)(Base_Lattice.nSites * (dSubBlock.x * dSubBlock.y * dSubBlock.z)*1.25);   
 
  for (i=0; i<Total_nSuperBlocks; i++) {
    tmpDomain[i].subDomain = (struct SubBlock *) malloc (sizeof(struct SubBlock)*Total_nSubBlocks);    
    for (j=0; j<Total_nSubBlocks; j++) {
      tmpDomain[i].subDomain[j].nSites=0;
      tmpDomain[i].subDomain[j].Sites = alloc_1Dint_array (k);             
    }
  }
   
#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks)
{
  #pragma omp for schedule(dynamic) private(i, iSup, iSub)  
    for (i=0; i<MCSystem.nSites; i++) { 
      Get_Block_Coordinate_for_Site__omp (&MCSystem.Sites[i].R, &nSuperBlocks, &nSubBlocks, &dSuperBlock, &dSubBlock, &iSup, &iSub);
  #pragma omp critical
      {tmpDomain[iSup].subDomain[iSub].Sites[tmpDomain[iSup].subDomain[iSub].nSites++] = i;}
    }
}
  
 
  Domain = (struct SuperBlock *) malloc (sizeof(struct SuperBlock)*Total_nSuperBlocks);           
  for (i=0; i<Total_nSuperBlocks; i++) {
    Domain[i].subDomain = (struct SubBlock *) malloc (sizeof(struct SubBlock)*Total_nSubBlocks);     
    for (j=0; j<Total_nSubBlocks; j++) {           
      Domain[i].subDomain[j].nSites = tmpDomain[i].subDomain[j].nSites;                           
      Domain[i].subDomain[j].Sites = alloc_1Dint_array (tmpDomain[i].subDomain[j].nSites);               
      for (k=0; k<tmpDomain[i].subDomain[j].nSites; k++) Domain[i].subDomain[j].Sites[k] = tmpDomain[i].subDomain[j].Sites[k];     
      
      free (tmpDomain[i].subDomain[j].Sites);                                                     
    }
    free (tmpDomain[i].subDomain);                                                                
  }
  free (tmpDomain);                                                                               
  
  fprintf (MC_LogFile, "Domain decomposition performed. Total number of SuperBlocks and SubBlocks(per SuperBlock) are %d and %d, respectively.\n", Total_nSuperBlocks, Total_nSubBlocks);
  fprintf (MC_LogFile, "\nNumber of active sites assigned to each processor and active sites per sub-block:\n");
  Tot=0;
  for (i=0; i<Total_nSuperBlocks; i++) {
    for (j=k=0; j<Total_nSubBlocks; j++) k+=Domain[i].subDomain[j].nSites;
    Tot+=k;
    fprintf (MC_LogFile, "Processor %d: %-7d ", i, k);
    fprintf (MC_LogFile, "Sites in each sub-block: "); 
    for (j=0; j<Total_nSubBlocks; j++) fprintf (MC_LogFile, "%d ", Domain[i].subDomain[j].nSites); fprintf (MC_LogFile, "\n");
  }
  
  if (Tot != MCSystem.nSites) 
    Exit_Error (" Error in Domain-decomposition: Total number of active sites allocated to all domains don't match the number of active sites in MCSystem. Check.");

  return;
}
 

void Get_Block_Coordinate_for_Site__omp (struct Coord *R, struct IntCoord *nsup, struct IntCoord *nsub, struct IntCoord *dsup, struct IntCoord *dsub, int *iSup, int *iSub)
{
  struct IntCoord tmp, smp;

  tmp.x=R->x/dsup->x;  
  tmp.y=R->y/dsup->y;  
  tmp.z=R->z/dsup->z; 
  (*iSup)=tmp.x + tmp.y * nsup->x + tmp.z * nsup->x * nsup->y;
  
  smp.x=(R->x - tmp.x*dsup->x)/dsub->x;  
  smp.y=(R->y - tmp.y*dsup->y)/dsub->y;
  smp.z=(R->z - tmp.z*dsup->z)/dsub->z;
  (*iSub)=smp.x + smp.y * nsub->x + smp.z * nsub->x * nsub->y;

  return;
}



void Create_ActiveSites_Domains ()
{
  int i, j, k, iSup, iSub;

  aSiteDomain = (struct SuperBlock *) malloc (sizeof(struct SuperBlock)*Total_nSuperBlocks); 
 
  for (i=0; i<Total_nSuperBlocks; i++) {
    aSiteDomain[i].subDomain = (struct SubBlock *) malloc (sizeof(struct SubBlock)*Total_nSubBlocks);    
    for (j=0; j<Total_nSubBlocks; j++) aSiteDomain[i].subDomain[j].nSites=0;
  }
 
#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks)
{
  #pragma omp for schedule(static, 1) private(i, iSup, iSub)  
    for (iSup=0; iSup<Total_nSuperBlocks; iSup++) { 
      for (iSub=0; iSub<Total_nSubBlocks; iSub++) 
        for (i=0; i<Domain[iSup].subDomain[iSub].nSites; i++) 
          if (MCSystem.Sites[Domain[iSup].subDomain[iSub].Sites[i]].Active) aSiteDomain[iSup].subDomain[iSub].nSites++;
    }
}

  for (i=0; i<Total_nSuperBlocks; i++)  for (j=0; j<Total_nSubBlocks; j++) aSiteDomain[i].subDomain[j].Sites = alloc_1Dint_array (aSiteDomain[i].subDomain[j].nSites);        
  
#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks)
{
  #pragma omp for schedule(static, 1) private(i, iSup, iSub, k)  
    for (iSup=0; iSup<Total_nSuperBlocks; iSup++) { 
      for (iSub=0; iSub<Total_nSubBlocks; iSub++) 
        for (i=0, k=0; i<Domain[iSup].subDomain[iSub].nSites; i++) 
          if (MCSystem.Sites[Domain[iSup].subDomain[iSub].Sites[i]].Active) aSiteDomain[iSup].subDomain[iSub].Sites[k++]=Domain[iSup].subDomain[iSub].Sites[i];
    }
}

  return;
}


 
int Number_of_Metropolis_Steps_perBlock (void)
{
  int i, j, k, kk=0;

  for (i=0; i<Total_nSuperBlocks; i++) {
    for (j=0, k=0; j<Total_nSubBlocks; j++) k+=aSiteDomain[i].subDomain[j].nSites;
    if (k > kk) kk=k;
  }
  return kk;
}



void Free_Global_Memory__omp (void)
{
  int i, j;

  free (rng_parallel);

  for (i=0; i<Total_nSuperBlocks; i++) {
    for (j=0; j<Total_nSubBlocks; j++) free (Domain[i].subDomain[j].Sites);
    for (j=0; j<Total_nSubBlocks; j++) free (aSiteDomain[i].subDomain[j].Sites);  
    free (Domain[i].subDomain);
    free (aSiteDomain[i].subDomain);
  }
  free (Domain);
  free (aSiteDomain);
}



void Seed_RNG__omp (void)                   
{
  int i;

  rng_parallel = (pcg32_random_t *) malloc(sizeof(pcg32_random_t)*Total_nSuperBlocks);
  for (i=0; i<Total_nSuperBlocks; i++)  pcg32_random_float_seed (&rng_parallel[i]);                     
}




double Configuration_Energy__omp (struct SuperCell *FullSys, struct Element *ElementList) 
{
  int i;
  double ene;
  
  ene = ECI_Vec[0];
  for (i=0; i<nClusters; i++) ene += Average_ClusterEnergy__omp (FullSys, &ClusterList[i], ElementList);
    
  return ene;
}



double Average_ClusterEnergy__omp (struct SuperCell *FullSys, struct ClusterFamily *inCluster, struct Element *ElementList) 
{
  int i, j, k, m, isite, iRef, ncls=0;
  int iSup, iSub;
  double **pcorrel=NULL, phi, ClsEne;
  int *Cluster_Atoms=NULL;

  pcorrel=alloc_2Ddouble_array (Total_nSuperBlocks, inCluster->nECIs);    
  for (i=0; i<Total_nSuperBlocks; i++) for (m=0; m<inCluster->nECIs; m++) pcorrel[i][m]=0.0;
   


#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks) private (Cluster_Atoms, i, j, k, m, phi, iSup, iSub, iRef, isite) 
{  
  Cluster_Atoms = alloc_1Dint_array (inCluster->Ksites);

#pragma omp for schedule (static, 1) reduction(+:ncls) 
  for (iSup=0; iSup<Total_nSuperBlocks; iSup++) {
    for (iSub=0; iSub<Total_nSubBlocks; iSub++) {
      
      for (i=0; i<Domain[iSup].subDomain[iSub].nSites; i++) {                                  
        iRef=Domain[iSup].subDomain[iSub].Sites[i];

        
        for (j=0; j<inCluster->nChildGroups; j++)                           
          if (inCluster->ChildGroup_List[j].BL_Index == FullSys->Sites[iRef].iBaseLat && inCluster->ChildGroup_List[j].nChild > 0) {         
            for (k=0; k<inCluster->ChildGroup_List[j].nChild; k++) {                                     
              ncls+=1;   
                
              Cluster_Occupations__omp (&(FullSys->Sites[iRef].R), &inCluster->ChildGroup_List[j].Child[k], inCluster->Ksites, FullSys->Sites, &FullSys->Ncell, Cluster_Atoms);
              
              for (m=0; m<inCluster->nECIs; m++) {                                                      
                for (isite=0, phi=1.0; isite<inCluster->Ksites; isite++) 
                  phi *= ElementList[Cluster_Atoms[isite]].Theta[inCluster->ECI_List[m].Index[isite]];  
                pcorrel[iSup][m] += phi;        
              }
            }            
            break;
          }
      }  
    }    
  }      

  free (Cluster_Atoms);
}     
  
  for (m=0; m<inCluster->nECIs; m++)                         
    for (iSup=1; iSup<Total_nSuperBlocks; iSup++) pcorrel[0][m]+=pcorrel[iSup][m];  
  
  ncls *= inCluster->Ksites;  
  for (m=0; m<inCluster->nECIs; m++) pcorrel[0][m]/=ncls;    

  for (m=0, ClsEne=0.0; m<inCluster->nECIs; m++)             
    ClsEne += inCluster->ECI_List[m].ECI * pcorrel[0][m];
  
  free_2Ddouble_array(pcorrel, Total_nSuperBlocks);
  return ClsEne;
}



void Cluster_Occupations__omp (struct Coord *R, struct Cluster *cls, int ksites, struct AtomLabeled_Coord *FullSys_Sites, struct Coord *ncell, int *Cluster_Occupation)
{
  int i;
  struct Cluster Displaced_Cluster;

  Displaced_Cluster.R = alloc_1DCoord_array (ksites);
  Displace_Cluster_PeriodicBoundary (R, cls, &Displaced_Cluster, ksites, ncell);

#ifdef OLD_HASH_KEY
  for (i=0; i<ksites; i++) Cluster_Occupation[i] = FullSys_Sites[Old_Fetch_Coordinate_Index(&Displaced_Cluster.R[i])].iAtom;
#else  
  for (i=0; i<ksites; i++) Cluster_Occupation[i] = FullSys_Sites[New_Fetch_Coordinate_Index(&Displaced_Cluster.R[i])].iAtom;
#endif

  free (Displaced_Cluster.R);
  return;
}
  




double Average_perSite_dE__omp (struct AtomLabeled_Coord *iSite, struct SuperCell *FullSys, struct Element *ElementList, int iRef, int iExch, int Ref_Atom, int Exch_Atom) 
{
  int i, j, k, m, kk;
  int *Ref_Occup=NULL, *Exch_Occup=NULL;                    
  double *Ref_correl=NULL, *Exch_correl=NULL;
  double phi, Ref_Ene, Exch_Ene;

  Ref_Occup  = alloc_1Dint_array (LargestKsites);   Ref_correl  = alloc_1Ddouble_array (maxClusterComponents);
  Exch_Occup = alloc_1Dint_array (LargestKsites);   Exch_correl = alloc_1Ddouble_array (maxClusterComponents);  

  Ref_Ene=Exch_Ene=0.0;
  
  for (i=0; i<nClusters; i++) {            

    for (j=0; j<maxClusterComponents; j++)  Ref_correl[j]=0.0;
    for (j=0; j<maxClusterComponents; j++)  Exch_correl[j]=0.0;
    
    for (j=0; j<ClusterList[i].nChildGroups; j++)                             
      if (ClusterList[i].ChildGroup_List[j].BL_Index == iSite->iBaseLat && ClusterList[i].ChildGroup_List[j].nChild > 0) {   
      
        for (k=0; k<ClusterList[i].ChildGroup_List[j].nChild; k++) {                                    
          
          Xchngd_Cluster_Occupations__omp (&(iSite->R), &ClusterList[i].ChildGroup_List[j].Child[k], ClusterList[i].Ksites, FullSys->Sites, &FullSys->Ncell, Ref_Occup, Exch_Occup, iRef, iExch, Ref_Atom, Exch_Atom);
          
          for (m=0; m<ClusterList[i].nECIs; m++) {                                                      
            for (kk=0, phi=1.0; kk<ClusterList[i].Ksites; kk++) 
              phi *= ElementList[Ref_Occup[kk]].Theta[ClusterList[i].ECI_List[m].Index[kk]];           
            Ref_correl[m] += phi;                                                                       
          }
          
          
          for (m=0; m<ClusterList[i].nECIs; m++) {                                                      
            for (kk=0, phi=1.0; kk<ClusterList[i].Ksites; kk++)
              phi *= ElementList[Exch_Occup[kk]].Theta[ClusterList[i].ECI_List[m].Index[kk]];           
            Exch_correl[m] += phi;                                                                      
          }          
        }
        for (m=0; m<ClusterList[i].nECIs; m++) 
          Ref_Ene += ClusterList[i].ECI_List[m].ECI * Ref_correl[m] / (ClusterList[i].Ksites * ClusterList[i].ChildGroup_List[j].nChild);
        
        for (m=0; m<ClusterList[i].nECIs; m++) 
          Exch_Ene += ClusterList[i].ECI_List[m].ECI * Exch_correl[m] / (ClusterList[i].Ksites * ClusterList[i].ChildGroup_List[j].nChild);  
        
        break;
      }                    
  }
  
  free (Ref_Occup);   free (Exch_Occup);
  free (Ref_correl);  free (Exch_correl);
  return (Exch_Ene-Ref_Ene);
}



unsigned int Site_Perturbed__omp (struct Coord *R0, struct Coord *R1, struct Coord *ncell, float *p)
{
  if (Periodic_Distance (R0, R1, ncell) < (*p)) return 1;
  else return 0;
}



void Xchngd_Cluster_Occupations__omp (struct Coord *R, struct Cluster *cls, int ksites, struct AtomLabeled_Coord *FullSys_Sites, struct Coord *ncell, int *Ref_Occup, int *Exch_Occup, \
                                      int iRef, int iExch, int Ref_Atom, int Exch_Atom)
{
  int i;
  struct Cluster Displaced_Cluster;

  Displaced_Cluster.R = alloc_1DCoord_array (ksites);
  Displace_Cluster_PeriodicBoundary (R, cls, &Displaced_Cluster, ksites, ncell);

#ifdef OLD_HASH_KEY  
  for (i=0; i<ksites; i++) Exch_Occup[i] = Old_Fetch_Coordinate_Index(&Displaced_Cluster.R[i]);     
#else
  for (i=0; i<ksites; i++) Exch_Occup[i] = New_Fetch_Coordinate_Index(&Displaced_Cluster.R[i]);     
#endif

  for (i=0; i<ksites; i++) Ref_Occup[i] = FullSys_Sites[Exch_Occup[i]].iAtom;                       

  for (i=0; i<ksites; i++)                                                                          
    if (Exch_Occup[i] == iRef) Exch_Occup[i] = Exch_Atom;      
    else if (Exch_Occup[i] == iExch) Exch_Occup[i] = Ref_Atom;
    else Exch_Occup[i] = FullSys_Sites[Exch_Occup[i]].iAtom;  

  free (Displaced_Cluster.R);
  return;
}
  


double Local_dE__omp (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, float *pdist)
{
  int i, Ref_Atom, Exch_Atom;
  double dE;

  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;

  for (i=0, dE=0.0; i<FullSys->nSites; i++) 
    if (Site_Perturbed__omp (&(FullSys->Sites[Ref_Site].R), &(FullSys->Sites[i].R), &(FullSys->Ncell), pdist)) 
      dE += Average_perSite_dE__omp (&(FullSys->Sites[i]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);

  return (dE/FullSys->nSites);
}



void MonteCarlo_LocalEnergy__omp (float TkB, int NHISTORY, int minEQUIL_STEPS, float EQUIL_TOL, int EQUIL_STEPS, float *EQUIL_TOL_ACHIEVED, int *EQUIL_STEPS_ACHIEVED, int CALC_STEPS, int BLOCK_SIZE, struct Stat *EnergyStat)
{
  int  MCSteps_Reached, nBlocks, Block_Counter=0, k, nHist_Max;                    
  double Ecurrent, *Ehist=NULL, *Eblock=NULL, *Spblock=NULL;
  double MeanE, StdE, MedianE, SpHeat, dSpHeat, tmpx, MeanStdE;
  int i, site_key, jump_site_key;
  double dE, Total_dE, X;
  int iSup, iSub;

  Ecurrent = Configuration_Energy__omp (&MCSystem, ElementList);
  if (CALC_STEPS) {                                     
    nHist_Max = BLOCK_SIZE;
    if (BLOCK_SIZE > 0) Ehist = alloc_1Ddouble_array (BLOCK_SIZE);
    else Exit_Error ("BLOCK_SIZE is 0 in MonteCarlo_LocalEnergy__omp. Exit.");

    nBlocks = ((int)(CALC_STEPS/BLOCK_SIZE));
    if (nBlocks == 0) nBlocks=1;
    Eblock = alloc_1Ddouble_array (nBlocks); 
    Spblock = alloc_1Ddouble_array (nBlocks);}
  else {                                               
    nHist_Max = NHISTORY;
    if (NHISTORY > 0) Ehist = alloc_1Ddouble_array (NHISTORY);
    else Exit_Error ("BLOCK_SIZE is 0 in MonteCarlo_LocalEnergy__omp. Exit.");
    MeanStdE=0.0;
    nBlocks=0;
  }
 
  MCSteps_Reached = 0; 

  
  while (1) {                     
        
    for (k=0; k<nHist_Max; k++) {     
      
      iSub=0;
      
      for (i=0; i<nMetropolis_per_SuperBlock; i++) {  
        Total_dE=0.0;
#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks) private (site_key, jump_site_key, X, i, iSup, dE) firstprivate (TkB)
{         
        #pragma omp for schedule (static, 1) reduction (+:Total_dE)                       
          for (iSup=0; iSup<Total_nSuperBlocks; iSup++) {
            dE=0.0;
            site_key = aSiteDomain[iSup].subDomain[iSub].Sites[pcg32_boundedrand(aSiteDomain[iSup].subDomain[iSub].nSites)];
       
            jump_site_key = MCSystem.Sites[site_key].Site_Neighbors[pcg32_boundedrand(MCSystem.Sites[site_key].nNeighbors)];  
            if (MCSystem.Sites[jump_site_key].iAtom == MCSystem.Sites[site_key].iAtom) continue;
            dE = Local_dE__omp (&MCSystem, ElementList, site_key, jump_site_key, &Perturbation_Distance);
    
            if (dE < -ENE_CUTOFF__OMP) {                                                            
              Total_dE+=dE;
              int_swap (&MCSystem.Sites[site_key].iAtom, &MCSystem.Sites[jump_site_key].iAtom);}    
            else if (dE > ENE_CUTOFF__OMP) {                                                        
              X = pcg32_random_float (&rng_parallel[iSup]);                                               
              if (X < exp(-(MCSystem.nSites * dE)/TkB)) {                                           
                Total_dE+=dE;
                int_swap (&MCSystem.Sites[site_key].iAtom, &MCSystem.Sites[jump_site_key].iAtom);}    
            }
            
          }
}         

        iSub++; iSub%=Total_nSubBlocks; 
        Ecurrent += Total_dE;
      }
      

      Ehist[k]=Ecurrent;
      MCSteps_Reached++;
    }  
    
    if (CALC_STEPS) {                        
      Get_Simple_Statistics_double (Ehist, nHist_Max, &MeanE, &MedianE, &StdE);
      Get_Specific_Heat (Ehist, nHist_Max, 1.0/TkB, &SpHeat);
      Eblock[Block_Counter]=MeanE;
      Spblock[Block_Counter]=SpHeat;
      Block_Counter++;
      if (Block_Counter == nBlocks) break;}
    else {                                  
      Get_Simple_Statistics_double (Ehist, nHist_Max, &MeanE, &MedianE, &StdE);
      (*EnergyStat) = (struct Stat) {.Mean=MeanE, .Median=MedianE, .StdDev=StdE};
      (*EQUIL_STEPS_ACHIEVED) = MCSteps_Reached;    
      
      tmpx = (MeanStdE*nBlocks)/(nBlocks+1.0) + StdE/(nBlocks+1.0);
      (*EQUIL_TOL_ACHIEVED) = fabs(tmpx-MeanStdE);
      
      if (nBlocks > minEQUIL_STEPS && EQUIL_TOL > 0.0 && (*EQUIL_TOL_ACHIEVED) < EQUIL_TOL) break;   
      MeanStdE=tmpx;
      nBlocks++;  

      fprintf (MC_LogFile, "Current Energy: %f   Mean: %f   Median: %f    StdE: %f    MeanStdE: %f\n", Ecurrent, MeanE, MedianE, StdE, MeanStdE);
      if (MCSteps_Reached == EQUIL_STEPS) break; 
    }      
  }       

  if (CALC_STEPS) {                                         
    Get_Simple_Statistics_double (Eblock, nBlocks, &MeanE, &MedianE, &StdE);
    Get_Simple_Statistics_double (Spblock, nBlocks, &SpHeat, &tmpx, &dSpHeat);
    (*EnergyStat) = (struct Stat) {.Mean=MeanE, .Median=MedianE, .StdDev=StdE, .SpHeat=SpHeat, .dSpHeat=dSpHeat};
    free (Eblock);
    free (Spblock);
  }
  
  free (Ehist); 
  return;
}


double Local_dE_wMemory__omp (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, int *psites, int np_sites)
{
  int i, Ref_Atom, Exch_Atom;
  double dE;

  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
  
  for (i=0, dE=0.0; i<np_sites; i++) dE += Average_perSite_dE__omp (&(FullSys->Sites[psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
  
  return (dE/FullSys->nSites);
}



void MonteCarlo_LocalEnergy_wMemory__omp (float TkB, int NHISTORY, int minEQUIL_STEPS, float EQUIL_TOL, int EQUIL_STEPS, float *EQUIL_TOL_ACHIEVED, int *EQUIL_STEPS_ACHIEVED, int CALC_STEPS, int BLOCK_SIZE, struct Stat *EnergyStat)
{
  int  MCSteps_Reached, nBlocks, Block_Counter=0, k, nHist_Max;                    
  double Ecurrent, *Ehist=NULL, *Eblock=NULL, *Spblock=NULL;
  double MeanE, StdE, MedianE, SpHeat, dSpHeat, tmpx, MeanStdE;
  int i, site_key, jump_site_key;
  double dE, Total_dE, X;
  int iSup, iSub;

  Ecurrent = Configuration_Energy__omp (&MCSystem, ElementList);
  if (CALC_STEPS) {                                     
    nHist_Max = BLOCK_SIZE;
    if (BLOCK_SIZE > 0) Ehist = alloc_1Ddouble_array (BLOCK_SIZE);
    else Exit_Error ("BLOCK_SIZE is 0 in MonteCarlo_LocalEnergy__omp. Exit.");

    nBlocks = ((int)(CALC_STEPS/BLOCK_SIZE));
    if (nBlocks == 0) nBlocks=1;
    Eblock = alloc_1Ddouble_array (nBlocks); 
    Spblock = alloc_1Ddouble_array (nBlocks);}
  else {                                               
    nHist_Max = NHISTORY;
    if (NHISTORY > 0) Ehist = alloc_1Ddouble_array (NHISTORY);
    else Exit_Error ("BLOCK_SIZE is 0 in MonteCarlo_LocalEnergy__omp. Exit.");
    nBlocks=0;
    MeanStdE=0.0;
  }
  
  MCSteps_Reached=0;


  while (1) {     
    
    for (k=0; k<nHist_Max; k++) {     
      
      iSub=0;
      for (i=0; i<nMetropolis_per_SuperBlock; i++) {  
        Total_dE=0.0;
#pragma omp parallel default(shared) num_threads (Total_nSuperBlocks) private (site_key, jump_site_key, X, iSup, dE) 
{      
        #pragma omp for schedule (static, 1) reduction (+:Total_dE)                                 
          for (iSup=0; iSup<Total_nSuperBlocks; iSup++) {
            site_key = aSiteDomain[iSup].subDomain[iSub].Sites[pcg32_boundedrand(aSiteDomain[iSup].subDomain[iSub].nSites)];  
            
            jump_site_key = MCSystem.Sites[site_key].Site_Neighbors[pcg32_boundedrand(MCSystem.Sites[site_key].nNeighbors)];  
            if (MCSystem.Sites[jump_site_key].iAtom == MCSystem.Sites[site_key].iAtom) continue;
            dE = Local_dE_wMemory__omp (&MCSystem, ElementList, site_key, jump_site_key, MCSystem.Sites[site_key].Perturbed_Sites, MCSystem.Sites[site_key].nPerturbed_Sites);
            
            if (dE < -ENE_CUTOFF__OMP) {                                                            
              Total_dE+=dE;
              int_swap (&MCSystem.Sites[site_key].iAtom, &MCSystem.Sites[jump_site_key].iAtom);}    
            else if (dE > ENE_CUTOFF__OMP) {                                                        
              X = pcg32_random_float (&rng_parallel[iSup]);                                               
              if (X < exp(-(MCSystem.nSites * dE)/TkB)) {                                           
                Total_dE+=dE;
                int_swap (&MCSystem.Sites[site_key].iAtom, &MCSystem.Sites[jump_site_key].iAtom);}                                                                       
            }
            
          }
}           
              
        iSub++; iSub%=Total_nSubBlocks; 
        Ecurrent += Total_dE;     
      }      
      
      Ehist[k]=Ecurrent;
      MCSteps_Reached++;
    }

    if (CALC_STEPS) {                       
      Get_Simple_Statistics_double (Ehist, nHist_Max, &MeanE, &MedianE, &StdE);
      Get_Specific_Heat (Ehist, nHist_Max, 1.0/TkB, &SpHeat);
      Eblock[Block_Counter]=MeanE;
      Spblock[Block_Counter]=SpHeat;
      Block_Counter++;
      if (Block_Counter == nBlocks) break;}
    else {                                  
      Get_Simple_Statistics_double (Ehist, nHist_Max, &MeanE, &MedianE, &StdE);
      (*EnergyStat) = (struct Stat) {.Mean=MeanE, .Median=MedianE, .StdDev=StdE};
      (*EQUIL_STEPS_ACHIEVED) = MCSteps_Reached;   

      tmpx = (MeanStdE*nBlocks)/(nBlocks+1.0) + StdE/(nBlocks+1.0);
      (*EQUIL_TOL_ACHIEVED) = fabs(tmpx-MeanStdE);
               
      if (nBlocks > minEQUIL_STEPS && EQUIL_TOL > 0.0 && (*EQUIL_TOL_ACHIEVED) < EQUIL_TOL) break;   
      MeanStdE=tmpx;
      nBlocks++;  
      fprintf (MC_LogFile, "Current Energy: %f   Mean: %f   Median: %f    StdE: %f    MeanStdE: %f\n", Ecurrent, MeanE, MedianE, StdE, MeanStdE); 
      if (MCSteps_Reached == EQUIL_STEPS) break;      
    }    
  }       
 
  if (CALC_STEPS) {                                         
    Get_Simple_Statistics_double (Eblock, nBlocks, &MeanE, &MedianE, &StdE);
    Get_Simple_Statistics_double (Spblock, nBlocks, &SpHeat, &tmpx, &dSpHeat);
    (*EnergyStat) = (struct Stat) {.Mean=MeanE, .Median=MedianE, .StdDev=StdE, .SpHeat=SpHeat, .dSpHeat=dSpHeat};
    free (Eblock);
    free (Spblock);  
  }
 
  free (Ehist);  
  return;
}




