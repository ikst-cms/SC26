//---------------------------------------------------------------------------------------- 
// A simple hash-function to look-up the coordinate of atoms in the system
//----------------------------------------------------------------------------------------

#define ROWPRIME 12289   
#define COLPRIME 769     

void Generate_New_HashTable (struct SuperCell *);
int *HashKey (struct Coord *);
void pHashKey (struct Coord *, int *);
struct IntCoord Coord2IntCoord (struct Coord *);
int New_Fetch_Coordinate_Index (struct Coord *);


float HASHCOORD_TOL=1.0e-4;
int COMPRESS_KEYSa[3], COMPRESS_KEYSb[3];
int HASH_TABLE[ROWPRIME][COLPRIME];


void Generate_New_HashTable (struct SuperCell *inSys)
{
  int i, j, *key;
  struct AtomLabeled_Coord *atmp=NULL;
  
  for (i=0; i<ROWPRIME; i++)
    for (j=0; j<COLPRIME; j++) HASH_TABLE[i][j]=-1;

  COMPRESS_KEYSa[0]=1;   COMPRESS_KEYSa[1]=(int)(inSys->Ncell.x);   COMPRESS_KEYSa[2]=((int)(inSys->Ncell.x)) * ((int)(inSys->Ncell.y));
  COMPRESS_KEYSb[0]=1;   COMPRESS_KEYSb[1]=100;                     COMPRESS_KEYSb[2]=10000;
  
  atmp=inSys->Sites;
  for (i=0; i<inSys->nSites; i++) {
    key=HashKey (&(atmp[i].R));
    if (HASH_TABLE[key[0]][key[1]] == -1) HASH_TABLE [key[0]][key[1]] = i;
    else Exit_Error ("Collision in hash-table. Change the table size.");
  }

  return;
}



int *HashKey (struct Coord *a)
{
  static int p[2];
  struct IntCoord r, r2;

  r=Coord2IntCoord (a);
  p[0] = (r.x + r.y*COMPRESS_KEYSa[1] + r.z*COMPRESS_KEYSa[2])%ROWPRIME;

  r2=(struct IntCoord ){.x=(int)(roundf((a->x-r.x)/HASHCOORD_TOL)), .y=(int)(roundf((a->y-r.y)/HASHCOORD_TOL)), .z=(int)(roundf((a->z-r.z)/HASHCOORD_TOL))};
  p[1] = (abs(r2.x) + abs(r2.y*COMPRESS_KEYSb[1]) + abs(r2.z*COMPRESS_KEYSb[2]))%COLPRIME;

  return p;
}



void pHashKey (struct Coord *a, int *p)
{
  struct IntCoord r, r2;

  r=Coord2IntCoord (a);
  p[0] = (r.x + r.y*COMPRESS_KEYSa[1] + r.z*COMPRESS_KEYSa[2])%ROWPRIME;

  r2=(struct IntCoord ){.x=(int)(roundf((a->x-r.x)/HASHCOORD_TOL)), .y=(int)(roundf((a->y-r.y)/HASHCOORD_TOL)), .z=(int)(roundf((a->z-r.z)/HASHCOORD_TOL))};
  p[1] = (abs(r2.x) + abs(r2.y*COMPRESS_KEYSb[1]) + abs(r2.z*COMPRESS_KEYSb[2]))%COLPRIME;

  return;
}



struct IntCoord Coord2IntCoord (struct Coord *a)
{
  return (struct IntCoord ){.x=(int)(roundf(a->x/HASHCOORD_TOL)*HASHCOORD_TOL), .y=(int)(roundf(a->y/HASHCOORD_TOL)*HASHCOORD_TOL), .z=(int)(roundf(a->z/HASHCOORD_TOL)*HASHCOORD_TOL)};
}



int New_Fetch_Coordinate_Index (struct Coord *r)
{
  int out[2], p;
  pHashKey(r, out);
  if ((p=HASH_TABLE[out[0]][out[1]]) == -1) Exit_Error ("Unable to get hash-key. Exit.");
  return p;
}



