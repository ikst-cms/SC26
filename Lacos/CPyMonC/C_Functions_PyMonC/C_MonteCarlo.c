//-----------------------------------------------------------------------------------------------------+
// Metropolis-based atom-exchange Monte Carlo code. This code is interfaced to Python code PyCMonC.    |
//                                                                                                     |
// NOTE:: Comparison of energy calculation between C-code and Python-code shows a difference in 5-th   |
//        decimal position which could be due to (double) precision used in C-code for energy          |
//                                                                                                     |
//-----------------------------------------------------------------------------------------------------+
 
#define ENE_CUTOFF 1.e-15

void Initialize_MonteCarlo (void);
void Free_Global_Memory (void);

void MC_Step_wGlobal_dE (double *, struct SuperCell *, struct Element *, float);
double Configuration_Energy (struct SuperCell *, struct Element *);
double Average_ClusterEnergy (struct SuperCell *, struct ClusterFamily *, struct Element *);
void Cluster_Occupations (struct Coord *, struct Cluster *, int, struct AtomLabeled_Coord *, struct Coord *, int *);

double Average_ClusterEnergy_LP (struct SuperCell *, struct ClusterFamily *, struct Element *);
void Site_Correlation_LP (struct Coord *, struct ClusterFamily *, struct Element *, int, struct AtomLabeled_Coord *, struct Coord *, int *, int *, double *);

double Average_perSite_dE (struct AtomLabeled_Coord *, struct SuperCell *, struct Element *, int, int, int, int); 
unsigned int Site_Perturbed (struct Coord *, struct Coord *, struct Coord *, float *);
void Xchngd_Cluster_Occupations (struct Coord *, struct Cluster *, int, struct AtomLabeled_Coord *, struct Coord *, int *, int *, int, int, int, int);

void MC_Step_wLocal_dE (double *, struct SuperCell *, struct Element *, float, float *);
double Local_dE (struct SuperCell *, struct Element *, int, int, float *);
double Local_dE_LP (struct SuperCell *, struct Element *, int, int, float *);

void MC_Step_wLocal_dE_wMemory (double *, struct SuperCell *, struct Element *, float);
double Local_dE_wMemory (struct SuperCell *, struct Element *, int, int, int *, int);
double Local_dE_wMemory_LP (struct SuperCell *, struct Element *, int, int, int *, int);
double Local_dE_wMemory_LP_Psite (struct SuperCell *, struct Element *, int, int);


void MonteCarlo (float, int, int, float, int, float *, int *, int, int, struct Stat *);

pcg32_random_t rng1;  

int *ActiveSites=NULL;           

void Initialize_MonteCarlo (void)
{
  int i, j;

  pcg32_random_float_seed (&rng1);  

  ActiveSites = alloc_1Dint_array (MCSystem.nActiveSites);
  for (i=0, j=0; i<MCSystem.nSites; i++) if (MCSystem.Sites[i].Active) ActiveSites[j++]=i;
  return;
}



void Free_Global_Memory (void)
{
  free (ActiveSites);
}


void MC_Step_wGlobal_dE (double *Ecurrent, struct SuperCell *FullSys, struct Element *ElementList, float TkB) 
{
  int i, site_key, jump_site_key;
  pcg32_random_t *rnglocal;                                                                  
  double Enew, dE, X;
  rnglocal = &rng1;                                                                          
  for (i=0; i<FullSys->nActiveSites; i++) {
    site_key = ActiveSites[pcg32_boundedrand(FullSys->nActiveSites)];
    if (! FullSys->Sites[site_key].Active) continue;
 
    jump_site_key = FullSys->Sites[site_key].Site_Neighbors[pcg32_boundedrand(FullSys->Sites[site_key].nNeighbors)]; 
    if (FullSys->Sites[jump_site_key].iAtom == FullSys->Sites[site_key].iAtom) continue;    
    int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);        
    Enew = Configuration_Energy (FullSys, ElementList);                                      
    dE = (Enew-(*Ecurrent))*FullSys->nSites;                                                 

    if (dE < -ENE_CUTOFF) (*Ecurrent)=Enew;                                                  
    else if (dE > ENE_CUTOFF) {                                                              
      X = pcg32_random_float (rnglocal);                                               
      if (X < exp (-dE/TkB)) (*Ecurrent)=Enew;
      else int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom); 
    }
   
  }
  return;
}


double Configuration_Energy (struct SuperCell *FullSys, struct Element *ElementList) 
{
  int i;
  double ene;
  
  ene = ECI_Vec[0];
  if (Loop_Parallelization) 
    for (i=0; i<nClusters; i++) ene +=  Average_ClusterEnergy_LP (FullSys, &ClusterList[i], ElementList);
  else                      
    for (i=0; i<nClusters; i++) ene +=  Average_ClusterEnergy (FullSys, &ClusterList[i], ElementList);

  return ene;
}


double Average_ClusterEnergy (struct SuperCell *FullSys, struct ClusterFamily *inCluster, struct Element *ElementList) 
{
  int i, j, k, m, isite, ncls;
  double *correl=NULL, phi, ClsEne;
  int *Cluster_Atoms=NULL;

  correl = alloc_1Ddouble_array (inCluster->nECIs);    
  for (m=0; m<inCluster->nECIs; m++) correl[m]=0.0;
  Cluster_Atoms = alloc_1Dint_array (inCluster->Ksites);
  
  for (i=0, ncls=0; i<FullSys->nSites; i++) {                 
    for (j=0; j<inCluster->nChildGroups; j++) {                                                       
      if (inCluster->ChildGroup_List[j].BL_Index == FullSys->Sites[i].iBaseLat && inCluster->ChildGroup_List[j].nChild > 0) {   
               
        for (k=0; k<inCluster->ChildGroup_List[j].nChild; k++) {                                      
          ncls+=1;
            
          Cluster_Occupations (&(FullSys->Sites[i].R), &inCluster->ChildGroup_List[j].Child[k], inCluster->Ksites, FullSys->Sites, &FullSys->Ncell, Cluster_Atoms);
          
          for (m=0; m<inCluster->nECIs; m++) {                                                        
            for (isite=0, phi=1.0; isite<inCluster->Ksites; isite++) 
               phi *= ElementList[Cluster_Atoms[isite]].Theta[inCluster->ECI_List[m].Index[isite]];   
            correl[m] += phi;        
          }
        }
        break;     
      }  
    }
  }
  
  ncls *= inCluster->Ksites;  
  for (m=0; m<inCluster->nECIs; m++) correl[m]/=ncls;
  for (m=0, ClsEne=0.0; m<inCluster->nECIs; m++) ClsEne += inCluster->ECI_List[m].ECI * correl[m];
  free (Cluster_Atoms); free (correl);
  return ClsEne;
}




void Cluster_Occupations (struct Coord *R, struct Cluster *cls, int ksites, struct AtomLabeled_Coord *FullSys_Sites, struct Coord *ncell, int *Cls_Occup)
{
  int i;
  struct Cluster Displaced_Cluster;

  Displaced_Cluster.R = alloc_1DCoord_array (ksites);
  Displace_Cluster_PeriodicBoundary (R, cls, &Displaced_Cluster, ksites, ncell);

#ifdef OLD_HASH_KEY
  for (i=0; i<ksites; i++) Cls_Occup[i] = FullSys_Sites[Old_Fetch_Coordinate_Index(&Displaced_Cluster.R[i])].iAtom;
#else
  for (i=0; i<ksites; i++) Cls_Occup[i] = FullSys_Sites[New_Fetch_Coordinate_Index(&Displaced_Cluster.R[i])].iAtom;
#endif

  free (Displaced_Cluster.R);
  return;
}



double Average_ClusterEnergy_LP (struct SuperCell *FullSys, struct ClusterFamily *inCluster, struct Element *ElementList) 
{
  int i, j, k, pncls, ncls=0, pchunk;
  double *Correl=NULL, *tCorrel=NULL, ClsEne;
  int *Cluster_Atoms=NULL;
  tCorrel = alloc_1Ddouble_array (inCluster->nECIs);  
  for (i=0; i<inCluster->nECIs; i++) tCorrel[i]=0.0;
  pchunk=(int)(FullSys->nSites/((float)(Max_OpenMP_Threads)));
  if (pchunk == 0) pchunk=1;


#if _OPENMP
#pragma omp parallel default(shared) num_threads (Max_OpenMP_Threads) private (Correl, Cluster_Atoms, k)
{     
  Correl = alloc_1Ddouble_array (inCluster->nECIs);        
  
  for (k=0; k<inCluster->nECIs; k++) Correl[k]=0.0;
  Cluster_Atoms = alloc_1Dint_array (inCluster->Ksites);
  
  #pragma omp for schedule(static, pchunk) private(i, j, pncls) reduction(+:ncls)
  for (i=0; i<FullSys->nSites; i++) {                                   
    Site_Correlation_LP (&(FullSys->Sites[i].R), inCluster, ElementList, FullSys->Sites[i].iBaseLat, FullSys->Sites, &FullSys->Ncell, &pncls, Cluster_Atoms, Correl);
    ncls+=pncls;
  }

  #pragma omp critical 
    for (j=0; j<inCluster->nECIs; j++) tCorrel[j]+=Correl[j];

  free(Correl);
  free (Cluster_Atoms);  
}
#else
  Exit_Error ("Called Average_ClusterEnergy_LP without OpenMP.");
#endif
  
  ncls *= inCluster->Ksites;  
  for (i=0; i<inCluster->nECIs; i++) tCorrel[i]/=ncls;
  
  for (i=0, ClsEne=0.0; i<inCluster->nECIs; i++) ClsEne += inCluster->ECI_List[i].ECI * tCorrel[i];

  return ClsEne;
}



void Site_Correlation_LP (struct Coord *R, struct ClusterFamily *inCluster, struct Element *ElementList, int iBaseLat, 
                          struct AtomLabeled_Coord *FullSys_Sites, struct Coord *ncell, int *ncls, int *cluster_atoms, double *correl)
{
  int j, k, m, isite;
  double phi; 
  (*ncls)=0;  
 for (j=0; j<inCluster->nChildGroups; j++) {                                
    if (inCluster->ChildGroup_List[j].BL_Index == iBaseLat && inCluster->ChildGroup_List[j].nChild > 0) {           
       for (k=0; k<inCluster->ChildGroup_List[j].nChild; k++) {          
        (*ncls)+=1;
          
        Cluster_Occupations (R, &inCluster->ChildGroup_List[j].Child[k], inCluster->Ksites, FullSys_Sites, ncell, cluster_atoms);
        for (m=0; m<inCluster->nECIs; m++) {                           
          for (isite=0, phi=1.0; isite<inCluster->Ksites; isite++)     
            phi *= ElementList[cluster_atoms[isite]].Theta[inCluster->ECI_List[m].Index[isite]];   
          correl[m] += phi;        
        }
      }
      break;     
    }  
  }
  
  return;
}


double Average_perSite_dE (struct AtomLabeled_Coord *iSite, struct SuperCell *FullSys, struct Element *ElementList, int iRef, int iExch, int Ref_Atom, int Exch_Atom) 
{
  int i, j, k, m, kk;
  int *Ref_Occup=NULL, *Exch_Occup=NULL;                    
  double *Ref_correl=NULL, *Exch_correl=NULL;
  double phi, Ref_Ene, Exch_Ene;
 
  Ref_Occup  = alloc_1Dint_array (LargestKsites);
  Exch_Occup = alloc_1Dint_array (LargestKsites);
  
  Ref_correl  = alloc_1Ddouble_array (maxClusterComponents);    
  Exch_correl = alloc_1Ddouble_array (maxClusterComponents);  
  Ref_Ene=Exch_Ene=0.0;
  for (i=0; i<nClusters; i++) {            
    for (j=0; j<maxClusterComponents; j++)  Ref_correl[j]=0.0;
    for (j=0; j<maxClusterComponents; j++)  Exch_correl[j]=0.0;
    
    for (j=0; j<ClusterList[i].nChildGroups; j++)                             
      if (ClusterList[i].ChildGroup_List[j].BL_Index == iSite->iBaseLat && ClusterList[i].ChildGroup_List[j].nChild > 0) {            
        
        for (k=0; k<ClusterList[i].ChildGroup_List[j].nChild; k++) {                           
          
          Xchngd_Cluster_Occupations (&(iSite->R), &ClusterList[i].ChildGroup_List[j].Child[k], ClusterList[i].Ksites, FullSys->Sites, &FullSys->Ncell, Ref_Occup, Exch_Occup, iRef, iExch, Ref_Atom, Exch_Atom);
          
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
  free (Ref_Occup);  free (Exch_Occup);
  free (Ref_correl); free (Exch_correl);
  return (Exch_Ene-Ref_Ene);
}



unsigned int Site_Perturbed (struct Coord *R0, struct Coord *R1, struct Coord *ncell, float *p)
{
  if (Periodic_Distance (R0, R1, ncell) < (*p)) return 1;
  else return 0;
}



void Xchngd_Cluster_Occupations (struct Coord *R, struct Cluster *cls, int ksites, struct AtomLabeled_Coord *FullSys_Sites, struct Coord *ncell, \
                                 int *Ref_Occup, int *Exch_Occup, int iRef, int iExch, int Ref_Atom, int Exch_Atom)
{
  int i;
  struct Cluster Displaced_Cluster;

  Displaced_Cluster.R = alloc_1DCoord_array (ksites);
  Displace_Cluster_PeriodicBoundary (R, cls, &Displaced_Cluster, ksites, ncell);
  
#ifdef OLD_HASH_KEY
  
  for (i=0; i<ksites; i++)
  {
 
  Exch_Occup[i] = Old_Fetch_Coordinate_Index(&Displaced_Cluster.R[i]);     
  }
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
  

void MC_Step_wLocal_dE (double *Ecurrent, struct SuperCell *FullSys, struct Element *ElementList, float TkB, float *pdist)  
{
  int i, site_key, jump_site_key;
  pcg32_random_t *rnglocal;                                                                  
  double dE, X;
  rnglocal = &rng1;                                                                          
  for (i=0; i<FullSys->nActiveSites; i++) {
    site_key = ActiveSites[pcg32_boundedrand(FullSys->nActiveSites)];
    if (! FullSys->Sites[site_key].Active) continue;
   
    jump_site_key = FullSys->Sites[site_key].Site_Neighbors[pcg32_boundedrand(FullSys->Sites[site_key].nNeighbors)];  
    if (FullSys->Sites[jump_site_key].iAtom == FullSys->Sites[site_key].iAtom) continue;
    
    if (Loop_Parallelization) dE = Local_dE_LP (FullSys, ElementList, site_key, jump_site_key, pdist);
    else                      dE = Local_dE (FullSys, ElementList, site_key, jump_site_key, pdist);

    if (dE < -ENE_CUTOFF) {                                                                            
      (*Ecurrent)+=dE;                                                           
      int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);}    
    else if (dE > ENE_CUTOFF) {                                                                                  
      X = pcg32_random_float (rnglocal);                                               
      if (X < exp ( -(FullSys->nSites * dE)/TkB) ) {                                         
        (*Ecurrent)+=dE;  
        int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);   
      }                                                                                      
    }
    
  }

  return;
}



double Local_dE (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, float *pdist)
{
  int i, Ref_Atom, Exch_Atom;
  double dE;

  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;

  for (i=0, dE=0.0; i<FullSys->nSites; i++) 
    if (Site_Perturbed (&(FullSys->Sites[Ref_Site].R), &(FullSys->Sites[i].R), &(FullSys->Ncell), pdist)) 
      dE += Average_perSite_dE (&(FullSys->Sites[i]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);

  return (dE/FullSys->nSites);
}




double Local_dE_LP (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, float *pdist)
{
  int i, Ref_Atom, Exch_Atom, pchunk;
  double dE=0.0;

  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;

  pchunk=(int)(FullSys->nSites/((float)(Max_OpenMP_Threads)));
  if (pchunk == 0) pchunk=1;
  
#if _OPENMP
  #pragma omp parallel default(shared) num_threads (Max_OpenMP_Threads)
{   
  #pragma omp for schedule (static, pchunk) private(i) reduction(+:dE)
  for (i=0; i<FullSys->nSites; i++) 
    if (Site_Perturbed (&(FullSys->Sites[Ref_Site].R), &(FullSys->Sites[i].R), &(FullSys->Ncell), pdist)) 
      dE += Average_perSite_dE (&(FullSys->Sites[i]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
}
#else
  Exit_Error ("Called Local_dE_LP without OpenMP");
#endif

  return (dE/FullSys->nSites);
}



void MC_Step_wLocal_dE_wMemory (double *Ecurrent, struct SuperCell *FullSys, struct Element *ElementList, float TkB)  
{
  int i, site_key, jump_site_key;
  pcg32_random_t *rnglocal;                                                                  
  double dE, X;
  rnglocal = &rng1;                                                                         
  for (i=0; i<FullSys->nActiveSites; i++) {
    site_key = ActiveSites[pcg32_boundedrand(FullSys->nActiveSites)];
    if (! FullSys->Sites[site_key].Active) continue;
    
    jump_site_key = FullSys->Sites[site_key].Site_Neighbors[pcg32_boundedrand(FullSys->Sites[site_key].nNeighbors)]; 
    if (FullSys->Sites[jump_site_key].iAtom == FullSys->Sites[site_key].iAtom) continue;

#ifdef CTIME    
    custom_swaps += 1;
#endif

    if (Loop_Parallelization) {
#ifdef OPT_PSITE
      dE = Local_dE_wMemory_LP_Psite (FullSys, ElementList, site_key, jump_site_key);
#else
      dE = Local_dE_wMemory_LP (FullSys, ElementList, site_key, jump_site_key, FullSys->Sites[site_key].Perturbed_Sites, FullSys->Sites[site_key].nPerturbed_Sites);
#endif
    }
    else {
      dE = Local_dE_wMemory (FullSys, ElementList, site_key, jump_site_key, FullSys->Sites[site_key].Perturbed_Sites, FullSys->Sites[site_key].nPerturbed_Sites);
    }

    if (dE < -ENE_CUTOFF) {                                                                            
      (*Ecurrent)+=dE;                                                           
      int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);}     
    else if (dE > ENE_CUTOFF) {                                                                                   
      X = pcg32_random_float (rnglocal);                                               
      if (X < exp (-(FullSys->nSites * dE)/TkB) ) {                                          
        (*Ecurrent)+=dE;  
        int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);    
      }                                                                                      
    }
    
  }
 
  return;
}



double Local_dE_wMemory (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, int *psites, int np_sites)
{
  int i, Ref_Atom, Exch_Atom;
  double dE;

  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
  for (i=0, dE=0.0; i<np_sites; i++) 
    dE += Average_perSite_dE (&(FullSys->Sites[psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
  
  return (dE/FullSys->nSites);
}



double Local_dE_wMemory_LP (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, int *psites, int np_sites)
{
  int i, Ref_Atom, Exch_Atom, pchunk;
  double dE=0.0;
  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
  pchunk = (int)(np_sites/((float)(Max_OpenMP_Threads)));

  if (pchunk == 0) pchunk=1;
  
#if _OPENMP
  #pragma omp parallel default(shared) num_threads (Max_OpenMP_Threads) 
{   
  #pragma omp for schedule (static, pchunk) private(i) reduction(+:dE)
  for (i=0; i<np_sites; i++) 
    dE += Average_perSite_dE (&(FullSys->Sites[psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
}

#else
  Exit_Error ("Called Local_wMemory_LP without OpenMP");
#endif
  
  return (dE/FullSys->nSites);
}

double Local_dE_wMemory_LP_Psite (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site)
{
  int np_sites = 0, psiteL, psiteR;

  for(psiteL=0,psiteR=0;
    (psiteL<FullSys->Sites[Ref_Site].nPerturbed_Sites)&&
    (psiteR<FullSys->Sites[Exch_Site].nPerturbed_Sites);
    np_sites++){
      if(FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL]
          <FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR]){
            total_psites[np_sites] = 
              FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL++];
      }
      else if(FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL]
          >FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR]){
            total_psites[np_sites] = 
              FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR++];
      }
      else{
            total_psites[np_sites] = 
              FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL];
            
            psiteL++;
            psiteR++;
      }
  }

  for(;psiteL<FullSys->Sites[Ref_Site].nPerturbed_Sites;
        np_sites++,psiteL++){
    total_psites[np_sites] = 
      FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL];
  }

  for(;psiteR<FullSys->Sites[Exch_Site].nPerturbed_Sites;
        np_sites++,psiteR++){
    total_psites[np_sites] = 
      FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR];
  } 
  
  int i, j, Ref_Atom, Exch_Atom, pchunk;
  double dE=0.0;
  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
  pchunk = (int)(np_sites/((float)(Max_OpenMP_Threads)));

  if (pchunk == 0) pchunk=1;
  
#if _OPENMP
  #pragma omp parallel default(shared) num_threads (Max_OpenMP_Threads) 
{   
  #pragma omp for schedule (static, pchunk) private(i) reduction(+:dE)
  for (i=0; i<np_sites; i++) 
    dE += Average_perSite_dE (&(FullSys->Sites[total_psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
}

#else
  Exit_Error ("Called Local_wMemory_LP without OpenMP");
#endif
  
  return (dE/FullSys->nSites);
}

void MonteCarlo (float TkB, int NHISTORY, int minEQUIL_STEPS, float EQUIL_TOL, int EQUIL_STEPS, float *EQUIL_TOL_ACHIEVED, int *EQUIL_STEPS_ACHIEVED, int CALC_STEPS, int BLOCK_SIZE, struct Stat *EnergyStat)
{
  int  MCSteps_Reached, k, nHist_Max, nBlocks, Block_Counter=0;                    
  double Ecurrent, *Ehist=NULL, *Eblock=NULL, *Spblock=NULL;
  double MeanE, StdE, MedianE, SpHeat, dSpHeat, tmp, MeanStdE;
  Ecurrent = Configuration_Energy (&MCSystem, ElementList); 
  if (CALC_STEPS) {                                          
    nHist_Max = BLOCK_SIZE;
    Ehist = alloc_1Ddouble_array (BLOCK_SIZE);

    nBlocks = ((int)(CALC_STEPS/BLOCK_SIZE));
    if (nBlocks == 0) nBlocks=1;
    Eblock = alloc_1Ddouble_array (nBlocks);
    Spblock = alloc_1Ddouble_array (nBlocks);}
  else {                                                    
    nHist_Max = NHISTORY;
    MeanStdE=0.0;
    nBlocks=0;
    Ehist = alloc_1Ddouble_array (NHISTORY);
  }
  
  MCSteps_Reached=0;
#ifdef TIME
  time_t Start_time,End_time;
  time(&Start_time);
#endif 

  if (Use_Local_EnergyDiff_Algo)            
    while (1) {                    
      if (Perturbed_Sites_List_Allocated) 
        for (k=0; k<nHist_Max; k++) {
#ifdef CTIME
          struct timeval start;	/* starting time */
          struct timeval end;	/* ending time */

          gettimeofday(&start, 0);
#endif
          MC_Step_wLocal_dE_wMemory (&Ecurrent, &MCSystem, ElementList, TkB);  
#ifdef CTIME
          gettimeofday(&end, 0);		/* mark the end time */

          /* now we can do the math. timeval has two elements: seconds and microseconds */
          custom_time += (((end.tv_sec * 1000000) + end.tv_usec) - ((start.tv_sec * 1000000) + start.tv_usec));
#endif
          Ehist[k]=Ecurrent;
          MCSteps_Reached++;
        }
      else
        for (k=0; k<nHist_Max; k++) {  
          MC_Step_wLocal_dE (&Ecurrent, &MCSystem, ElementList, TkB, &Perturbation_Distance);    
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
      
      tmp = (MeanStdE*nBlocks)/(nBlocks+1.0) + StdE/(nBlocks+1.0);
      (*EQUIL_TOL_ACHIEVED) = fabs(tmp-MeanStdE);
      
      if (nBlocks > minEQUIL_STEPS && EQUIL_TOL > 0.0 && (*EQUIL_TOL_ACHIEVED) < EQUIL_TOL) break;    
      MeanStdE=tmp;
      nBlocks++;  
      fprintf (MC_LogFile, "Current Energy: %f   Mean: %f   Median: %f    StdE: %f    MeanStdE: %f\n", Ecurrent, MeanE, MedianE, StdE, MeanStdE);
      if (MCSteps_Reached == EQUIL_STEPS) break;
      }
    }                  
    
  else                                       

    while (1) {                      
      for (k=0; k<nHist_Max; k++) {
        MC_Step_wGlobal_dE (&Ecurrent, &MCSystem, ElementList, TkB);     
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
        
        tmp = (MeanStdE*nBlocks)/(nBlocks+1.0) + StdE/(nBlocks+1.0);
        (*EQUIL_TOL_ACHIEVED) = fabs(tmp-MeanStdE);

        if (nBlocks > minEQUIL_STEPS && EQUIL_TOL > 0.0 && (*EQUIL_TOL_ACHIEVED) < EQUIL_TOL) break;  
        MeanStdE=tmp;
        nBlocks++;

        fprintf (MC_LogFile, "Current Energy: %f   Mean: %f   Median: %f    StdE: %f    MeanStdE: %f\n", Ecurrent, MeanE, MedianE, StdE, MeanStdE);
        if (MCSteps_Reached == EQUIL_STEPS) break;
      }
    }
  
  if (CALC_STEPS) {
    Get_Simple_Statistics_double (Eblock, nBlocks, &MeanE, &MedianE, &StdE);
    Get_Simple_Statistics_double (Spblock, nBlocks, &SpHeat, &tmp, &dSpHeat);
    (*EnergyStat) = (struct Stat) {.Mean=MeanE, .Median=MedianE, .StdDev=StdE, .SpHeat=SpHeat, .dSpHeat=dSpHeat};
    free (Eblock);
    free (Spblock);
  }

#ifdef TIME
  time(&End_time);
  fprintf (MC_LogFile, "\nTime taken for %d Metropolis equilibriation steps with single processor: %f hrs\n", MCSteps_Reached, ((float)(End_time-Start_time))/3600.00);
#endif 

  free (Ehist); 
  return;
}




