//---------------------------------------------------------------------------------------- 
// A simple hash-function to look-up the coordinate of atoms in the system
//----------------------------------------------------------------------------------------


void Generate_Old_HashTable (struct SuperCell *, struct UnitCell *);
void Generate_BaseLattice_HashTable (struct UnitCell *);

int Minimum_GridSize_for_Key (struct Labeled_Coord *, int);
int Get_UnitCell_HashKey (struct Coord *);
void Get_SuperCell_HashKey (struct Coord *, int *);
int Old_Fetch_Coordinate_Index (struct Coord *);


int KEYS1[3], KEYS2[3];
int *UnitCell_HashTable=NULL, **SuperCell_HashTable=NULL;
int UnitCell_GridSize;


void Generate_Old_HashTable (struct SuperCell *inSys, struct UnitCell *BaseSystem)
{
  int i, j, out[2];
  int Nmain_list, Ncell_list;
  struct AtomLabeled_Coord *atmp=NULL;
  Generate_BaseLattice_HashTable (BaseSystem);
  KEYS1[0]=1;   KEYS1[1]=(int)(inSys->Ncell.x);  KEYS1[2]=((int)(inSys->Ncell.x)) * ((int)(inSys->Ncell.y));
  Nmain_list = ((int)(inSys->Ncell.x)) * ((int)(inSys->Ncell.y)) * ((int) (inSys->Ncell.z));
  Ncell_list = BaseSystem->nSites;
  SuperCell_HashTable = alloc_2Dint_array (Nmain_list, Ncell_list);   
  for (i=0; i<Nmain_list; i++) for (j=0; j<Ncell_list; j++) SuperCell_HashTable[i][j]=-1;
  
  atmp=inSys->Sites;
  for (i=0; i<inSys->nSites; i++){
    Get_SuperCell_HashKey (&(atmp[i].R), out); 
    SuperCell_HashTable [out[0]][out[1]] = i;
  }

  return;
}



void Generate_BaseLattice_HashTable (struct UnitCell *BaseSystem)
{
  int i, nLEN, itmp;
  struct Labeled_Coord *atmp=NULL;

  UnitCell_GridSize = Minimum_GridSize_for_Key (BaseSystem->Sites, BaseSystem->nSites);    

  if (UnitCell_GridSize == 1) {KEYS2[0]=1; KEYS2[1]=1; KEYS2[2]=1;}         
  else { 
    KEYS2[0]=1; KEYS2[1]=UnitCell_GridSize; KEYS2[2]=UnitCell_GridSize * UnitCell_GridSize;  
  }

  nLEN = UnitCell_GridSize * UnitCell_GridSize * UnitCell_GridSize;
  UnitCell_HashTable = alloc_1Dint_array (nLEN);  
  for (i=0; i<nLEN; i++) UnitCell_HashTable[i]=-1;
  
  atmp=BaseSystem->Sites;
  for (i=0; i<BaseSystem->nSites; i++) {
    itmp = Get_UnitCell_HashKey (&atmp[i].R);

    if (UnitCell_HashTable[itmp] != -1) Exit_Error ("\n Error in Generate_Atom_KeyList_Base_Lattice. Possible collision --> change the KEYS Exit\n"); 
    else UnitCell_HashTable[itmp]=i;
  }

  return;
}


int Get_UnitCell_HashKey (struct Coord *r)
{
  int ix, iy, iz;

  ix = (int)(round(r->x * UnitCell_GridSize));   
  iy = (int)(round(r->y * UnitCell_GridSize)); 
  iz = (int)(round(r->z * UnitCell_GridSize)); 
  
  return ix*KEYS2[0] + iy*KEYS2[1] + iz*KEYS2[2];
}


void Get_SuperCell_HashKey (struct Coord *r, int *out)
{
  int ix, iy, iz;
  struct Coord rtmp;
   
  rtmp.x=r->x - (ix=(int)(r->x)); if (fabs(rtmp.x-1.0) < COORD_TOL) {ix+=1; rtmp.x=0.0;}  
  rtmp.y=r->y - (iy=(int)(r->y)); if (fabs(rtmp.y-1.0) < COORD_TOL) {iy+=1; rtmp.y=0.0;}
  rtmp.z=r->z - (iz=(int)(r->z)); if (fabs(rtmp.z-1.0) < COORD_TOL) {iz+=1; rtmp.z=0.0;}
  out[0] = iz*KEYS1[2] + iy*KEYS1[1] + ix;                      
  out[1] = UnitCell_HashTable[Get_UnitCell_HashKey (&rtmp)];    

  return;
 }  


int Minimum_GridSize_for_Key (struct Labeled_Coord *inList, int nList)         
{
  int i, j, n;
  float d;
  struct Coord dR;

  d=LARGE;
  for (i=0; i<nList; i++) 
    for (j=i+1; j<nList; j++) {
      dR=Subtract_Coords (&inList[i].R, &inList[j].R);
      if (fabs(dR.x) < d && fabs(dR.x) > COORD_TOL) d=fabs(dR.x); 
      if (fabs(dR.y) < d && fabs(dR.y) > COORD_TOL) d=fabs(dR.y); 
      if (fabs(dR.z) < d && fabs(dR.z) > COORD_TOL) d=fabs(dR.z);
    }

  n=1;
  while (1) {
    if ((d*n) > 1) break; else n*=10;
  };

  return n;
}


int Old_Fetch_Coordinate_Index (struct Coord *r)
{
  int out[2];

  Get_SuperCell_HashKey (r, out); 
  return (SuperCell_HashTable [out[0]][out[1]]);
}



