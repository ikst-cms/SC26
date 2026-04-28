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

int *ActiveSites=NULL;           

void Initialize_MonteCarlo (void)
{
  int i, j;

  ActiveSites = alloc_1Dint_array (MCSystem.nActiveSites);
  for (i=0, j=0; i<MCSystem.nSites; i++) if (MCSystem.Sites[i].Active) ActiveSites[j++]=i;
  return;
}



void Free_Global_Memory (void)
{
  free (ActiveSites);
}


double Configuration_Energy (struct SuperCell *FullSys, struct Element *ElementList) 
{
  int i;
  double ene;
  
  ene = ECI_Vec[0];
  for (i=0; i<nClusters; i++) ene +=  Average_ClusterEnergy_LP (FullSys, &ClusterList[i], ElementList);

  return ene;
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
    
#ifdef OPT_PRE_CLUSTER
    j = iSite->correspondingCluster[i];
      if(j!=-1){
#else
    for (j=0; j<ClusterList[i].nChildGroups; j++)                             
      if (ClusterList[i].ChildGroup_List[j].BL_Index == iSite->iBaseLat && ClusterList[i].ChildGroup_List[j].nChild > 0) {            
#endif        
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

#ifndef OPT_PRE_CLUSTER
        break;
#endif
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

int getR(int rank, int size, int tot){
  int x = (tot%size);
  int y = tot/size;
  int z = (rank+1)*y;
  if(rank<x){
    return z+rank+1;
  }
  else{
    return z+x;
  }
}

int getL(int rank, int size, int tot){
  if(rank==0)return 0;
  return getR(rank-1, size, tot);
}

void executeList(struct SuperCell *FullSys, struct Element *ElementList, double *Ecurrent, float TkB){
  
  double start = MPI_Wtime();
  
  int curList = psp_curList%psp_maxList,i,j,psiteL,psiteR,total_psites=0;

#ifdef CLIST
  if(grank==0){
    FILE *CUSTOM = fopen("custom_list_size","a");
    fprintf(CUSTOM,"%d\n",psp_curCap[curList]);
    fclose(CUSTOM);
  }
#endif

  psp_actCap = 0;
  //int curList = psp_curList%psp_maxList,i,j,psiteL,psiteR,total_psites=0;
  for(i=0;i<psp_curCap[curList];i++){
    int site_key = psp_list_jump[curList][i];
    int jump_site_key = psp_list_neigh[curList][i];
    if (FullSys->Sites[jump_site_key].iAtom == FullSys->Sites[site_key].iAtom) continue;
    psp_act_jump[psp_actCap] = site_key;
    psp_act_neigh[psp_actCap] = jump_site_key;

    // Sites Per Swap
    int sites_per_swap = 0;

#ifdef OPT_PSITE 
    for(psiteL=0,psiteR=0;
      (psiteL<FullSys->Sites[site_key].nPerturbed_Sites)&&
      (psiteR<FullSys->Sites[jump_site_key].nPerturbed_Sites);
      total_psites++){
        if(FullSys->Sites[site_key].Perturbed_Sites[psiteL]
            <FullSys->Sites[jump_site_key].Perturbed_Sites[psiteR]){
              psp_act_psites[total_psites] = 
                FullSys->Sites[site_key].Perturbed_Sites[psiteL++];
        }
        else if(FullSys->Sites[site_key].Perturbed_Sites[psiteL]
            >FullSys->Sites[jump_site_key].Perturbed_Sites[psiteR]){
              psp_act_psites[total_psites] = 
                FullSys->Sites[jump_site_key].Perturbed_Sites[psiteR++];
        }
        else{
              psp_act_psites[total_psites] = 
                FullSys->Sites[site_key].Perturbed_Sites[psiteL];
              
              psiteL++;
              psiteR++;
        }
        psp_swap_ind[total_psites] = psp_actCap;
        sites_per_swap++;
    }

    for(;psiteL<FullSys->Sites[site_key].nPerturbed_Sites;
          total_psites++,psiteL++){
      psp_act_psites[total_psites] = 
        FullSys->Sites[site_key].Perturbed_Sites[psiteL];
      psp_swap_ind[total_psites] = psp_actCap;
      sites_per_swap++;
    }

    for(;psiteR<FullSys->Sites[jump_site_key].nPerturbed_Sites;
          total_psites++,psiteR++){
      psp_act_psites[total_psites] = 
        FullSys->Sites[jump_site_key].Perturbed_Sites[psiteR];
      psp_swap_ind[total_psites] = psp_actCap;
      sites_per_swap++;
    } 
#else
    total_psites += FullSys->Sites[site_key].nPerturbed_Sites;  
    sites_per_swap = FullSys->Sites[site_key].nPerturbed_Sites;
#endif 
    psp_local_energy[psp_actCap] = psp_final_energy[psp_actCap] = 0;
    ++psp_actCap;
  }

#ifdef CTIME
  if(grank==0){
    custom_swaps += psp_actCap;
  }
#endif

  if(psp_actCap!=0){
    int L = getL(grank,gsize,total_psites);
    int R = getR(grank,gsize,total_psites);
#ifdef OPT_PSITE
    for(j=L;j<R;j++){
      i = psp_swap_ind[j];
      int Ref_Site = psp_act_jump[i];
      int Exch_Site = psp_act_neigh[i];
      int Ref_Atom = FullSys->Sites[Ref_Site].iAtom;
      int Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
      psp_local_energy[i] += Average_perSite_dE (
          &(FullSys->Sites[psp_act_psites[j]]), FullSys, ElementList, 
          Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);
    }  
#else
    for(i=0;(i<psp_actCap)&&(L<R);i++){
      int nt = FullSys->Sites[psp_act_jump[i]].nPerturbed_Sites;
      int *psites = FullSys->Sites[psp_act_jump[i]].Perturbed_Sites;
      int Ref_Site = psp_act_jump[i];
      int Exch_Site = psp_act_neigh[i];
      int Ref_Atom = FullSys->Sites[Ref_Site].iAtom;
      int Exch_Atom = FullSys->Sites[Exch_Site].iAtom;
      while((L<R)&&(L<nt)){
        psp_local_energy[i] += Average_perSite_dE (
          &(FullSys->Sites[psites[L]]), FullSys, ElementList, 
          Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);

        ++L;
      }
      L-=nt;
      R-=nt;
    } 
#endif  

    // To get the maximum time on a rank for executing swaps
    double mid = MPI_Wtime();
    double local_elapsed = mid - start;    
    double global_elapsed = 0;
    MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    

    MPI_Allreduce(psp_local_energy, psp_final_energy, psp_actCap, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double end = MPI_Wtime();
    double local_combine = end - mid;
      
    #ifdef CLIST
    if(grank==0){
      FILE *TEMP = fopen("list_profile","a");
      fprintf(TEMP,"%d, %f, %f\n",psp_curCap[curList],global_elapsed, local_combine);
      fclose(TEMP);
    }
    #endif
  


    for(i=0;i<psp_actCap;i++){
      int should_swap;
      double dE = psp_final_energy[i]/FullSys->nSites;
      if(dE < -ENE_CUTOFF) {
        should_swap = 1;
      }
      else{
        double X = pcg32_random_float (&rng_select);                                               
        if (X < exp (-(FullSys->nSites * dE)/TkB) ) {     
          should_swap = 1;                                
        }  
        else {
          should_swap = 0;
        }
      } 
      if(should_swap){
        int site_key = psp_act_jump[i];
        int jump_site_key = psp_act_neigh[i];
        int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom); 
        (*Ecurrent)+=dE;
      }
    }
  }

  psp_curCap[curList] = 0;
  ++psp_curList;
}

void MC_Step_wLocal_dE_wMemory_Hybrid (double *Ecurrent, struct SuperCell *FullSys, struct Element *ElementList, float TkB)  
{
  int i,j, site_key, jump_site_key;                                                             
  double dE, X;

  for(i=0; i<FullSys->nSites; i++) {
    psp_labels[i] = 0;
  }                   

  for(i=0; i<psp_maxList; i++) {
    psp_curCap[i] = 0;
  }              

  psp_curList = psp_totalList = 0;

  for (i=0; i<FullSys->nActiveSites; i++) {
    site_key = pcg32_boundedrand_r(&rng,FullSys->nActiveSites);
    if (! FullSys->Sites[site_key].Active) continue;

    jump_site_key = FullSys->Sites[site_key].Site_Neighbors[pcg32_boundedrand_r(&rng,FullSys->Sites[site_key].nNeighbors)];

    int curIndex = psp_labels[site_key];
    int jumpIndex = psp_labels[jump_site_key];

    if(jumpIndex>curIndex){
      curIndex = jumpIndex;
    }

    while((psp_curCap[curIndex%psp_maxList]==psp_maxCap)&&(curIndex<psp_totalList)){
      ++curIndex;
    }

    if((curIndex-psp_curList)==psp_maxList){
      executeList(FullSys, ElementList, Ecurrent, TkB);
    }

    int redIndex = curIndex%psp_maxList;
    psp_list_jump[redIndex][psp_curCap[redIndex]] = site_key;
    psp_list_neigh[redIndex][psp_curCap[redIndex]] = jump_site_key;

    ++psp_curCap[redIndex];
    ++curIndex;
    if(curIndex>psp_totalList){
      psp_totalList = curIndex;
    }

    for(j=0;j<FullSys->Sites[site_key].nPerturbed_Sites; j++) {
      int *cur = &psp_labels[FullSys->Sites[site_key].Perturbed_Sites[j]];
      if(curIndex>(*cur)){
        *cur = curIndex;
      }
    }
#ifdef OPT_PSITE 
    for(j=0;j<FullSys->Sites[jump_site_key].nPerturbed_Sites; j++) {
      int *cur = &psp_labels[FullSys->Sites[jump_site_key].Perturbed_Sites[j]];
      if(curIndex>(*cur)){
        *cur = curIndex;
      }
    }
#endif
  }

  // completing remaining lists
  while(psp_curList<psp_totalList){
    executeList(FullSys, ElementList, Ecurrent, TkB);
  }

  return;
}

void MC_Step_wLocal_dE_wMemory (double *Ecurrent, struct SuperCell *FullSys, struct Element *ElementList, float TkB)  
{
  int i,j, site_key, jump_site_key;
  pcg32_random_t *rnglocal;                                                                  
  double dE, X;                                                                
  for (i=0; i<FullSys->nActiveSites; i++) {
    site_key = ActiveSites[pcg32_boundedrand_r(&rng,FullSys->nActiveSites)];
    if (! FullSys->Sites[site_key].Active) continue;
    jump_site_key = FullSys->Sites[site_key].Site_Neighbors[pcg32_boundedrand_r(&rng,FullSys->Sites[site_key].nNeighbors)]; 
    if (FullSys->Sites[jump_site_key].iAtom == FullSys->Sites[site_key].iAtom) continue;

#ifdef CTIME
    custom_swaps++;
#endif

#ifdef OPT_PSITE
    double lDE = Local_dE_wMemory_LP_Psite (FullSys, ElementList, site_key, jump_site_key);
#else
    double lDE = Local_dE_wMemory_LP (FullSys, ElementList, site_key, jump_site_key, FullSys->Sites[site_key].Perturbed_Sites, FullSys->Sites[site_key].nPerturbed_Sites);
#endif

    dE = 0;
    MPI_Allreduce(&lDE, &dE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(dE < -ENE_CUTOFF) {
      (*Ecurrent)+=dE;                                                           
      int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom); 
    }
    else{
      X = pcg32_random_float(&rng_select);                                               
      if (X < exp (-(FullSys->nSites * dE)/TkB) ) {     
        (*Ecurrent)+=dE;                                                           
        int_swap (&FullSys->Sites[site_key].iAtom, &FullSys->Sites[jump_site_key].iAtom);                                    
      }  
    }
  }
 
  return;
}

double Local_dE_wMemory_LP (struct SuperCell *FullSys, struct Element *ElementList, int Ref_Site, int Exch_Site, int *psites, int np_sites)
{
  int i, j, Ref_Atom, Exch_Atom, pchunk;
  double dE=0.0;
  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;

  j = getR(grank, gsize, np_sites);

  for (i=getL(grank, gsize, np_sites); i<j; i++) 
    dE += Average_perSite_dE (&(FullSys->Sites[psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);

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
            mpi_total_psites[np_sites] = 
              FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL++];
      }
      else if(FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL]
          >FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR]){
            mpi_total_psites[np_sites] = 
              FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR++];
      }
      else{
            mpi_total_psites[np_sites] = 
              FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL];
            
            psiteL++;
            psiteR++;
      }
  }

  for(;psiteL<FullSys->Sites[Ref_Site].nPerturbed_Sites;
        np_sites++,psiteL++){
    mpi_total_psites[np_sites] = 
      FullSys->Sites[Ref_Site].Perturbed_Sites[psiteL];
  }

  for(;psiteR<FullSys->Sites[Exch_Site].nPerturbed_Sites;
        np_sites++,psiteR++){
    mpi_total_psites[np_sites] = 
      FullSys->Sites[Exch_Site].Perturbed_Sites[psiteR];
  } 
  
  int i, j, Ref_Atom, Exch_Atom, pchunk;
  double dE=0.0;
  Ref_Atom  = FullSys->Sites[Ref_Site].iAtom;     
  Exch_Atom = FullSys->Sites[Exch_Site].iAtom;

  j = getR(grank, gsize, np_sites);

  for (i=getL(grank, gsize, np_sites); i<j; i++) 
    dE += Average_perSite_dE (&(FullSys->Sites[mpi_total_psites[i]]), FullSys, ElementList, Ref_Site, Exch_Site, Ref_Atom, Exch_Atom);

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
  if(grank==0) time(&Start_time);
#endif 
  

  while (1) {                    
    for (k=0; k<nHist_Max; k++) {

#ifdef CTIME
      struct timeval start;	/* starting time */
      struct timeval end;	/* ending time */

      if(grank==0){
        gettimeofday(&start, 0);
      }
#endif

#ifdef Hybrid
      MC_Step_wLocal_dE_wMemory_Hybrid (&Ecurrent, &MCSystem, ElementList, TkB);
#else
      MC_Step_wLocal_dE_wMemory (&Ecurrent, &MCSystem, ElementList, TkB);
#endif

#ifdef CTIME
      if(grank==0){
        gettimeofday(&end, 0);		/* mark the end time */

        /* now we can do the math. timeval has two elements: seconds and microseconds */
        custom_time += (((end.tv_sec * 1000000) + end.tv_usec) - ((start.tv_sec * 1000000) + start.tv_usec));
      }
#endif

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
      
      if(grank==0){
        fprintf (MC_LogFile, "Current Energy: %f   Mean: %f   Median: %f    StdE: %f    MeanStdE: %f\n", Ecurrent, MeanE, MedianE, StdE, MeanStdE);
      }
  
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
  if(grank==0){
    time(&End_time);
    fprintf (MC_LogFile, "\nTime taken for %d Metropolis equilibriation steps with single processor: %f hrs\n", MCSteps_Reached, ((float)(End_time-Start_time))/3600.00);
  }
#endif   
  
  free (Ehist);
  return;
}




