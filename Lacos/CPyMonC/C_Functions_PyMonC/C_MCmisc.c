//---------------------------------------------------------------------------------------- 
// Miscellaneous functions used in Monte Carlo routine attached to Python.
//----------------------------------------------------------------------------------------

struct Coord *alloc_1DCoord_array (int);

float Periodic_Distance (struct Coord *, struct Coord *, struct Coord *);
float Norm3d (struct Coord *);
struct Coord Add_Coords (struct Coord *, struct Coord *); 
struct Coord Subtract_Coords (struct Coord *, struct Coord *);
struct Coord ScalarMultiply_Coords (struct Coord *, float);
struct Coord Wrap_Coordinates_PeriodicBoundary (struct Coord *, struct Coord *);
float Distance3D (struct Coord *, struct Coord *);

void Copy_Cluster (struct Cluster *, struct Cluster *, int);
void Displace_Cluster (struct Coord *, struct Cluster *, struct Cluster *, int);
void Displace_Cluster_PeriodicBoundary (struct Coord *, struct Cluster *, struct Cluster *, int, struct Coord *);

void Get_Specific_Heat (double *, int, float, double *);


struct Coord *alloc_1DCoord_array (int k)
{
  struct Coord *a=NULL;

  if (k > 0) 
    if (!(a = (struct Coord *) malloc (k*sizeof(struct Coord)))) Exit_Error ("Unable to allocate 1D Coord array.");
  return (a); 
}



float Periodic_Distance (struct Coord *R1, struct Coord *R2, struct Coord *ncell)
{
  struct Coord dR;

  dR=Subtract_Coords (R1, R2);
  dR.x=fabsf(dR.x); dR.y=fabsf(dR.y); dR.z=fabsf(dR.z);

  dR.x = (dR.x-ncell->x/2.0)>COORD_TOL ? dR.x - (float)(ncell->x) : dR.x;
  dR.y = (dR.y-ncell->y/2.0)>COORD_TOL ? dR.y - (float)(ncell->y) : dR.y;
  dR.z = (dR.z-ncell->z/2.0)>COORD_TOL ? dR.z - (float)(ncell->z) : dR.z;
  
  return Norm3d(&dR);
}


float Norm3d (struct Coord *in)
{
  return (sqrt(in->x*in->x + in->y*in->y + in->z*in->z));
}


struct Coord Subtract_Coords (struct Coord *in1, struct Coord *in2)
{
  return ((struct Coord) {.x=in1->x-in2->x, .y=in1->y-in2->y, .z=in1->z-in2->z});
}


struct Coord Add_Coords (struct Coord *in1, struct Coord *in2)
{
  return ((struct Coord) {.x=in1->x+in2->x, .y=in1->y+in2->y, .z=in1->z+in2->z}); 
}


struct Coord ScalarMultiply_Coords (struct Coord *in1, float r)
{
  return ((struct Coord) {.x=in1->x*r, .y=in1->y*r, .z=in1->z*r}); 
}


float Distance3D (struct Coord *R1, struct Coord *R2)
{
  struct Coord a=Subtract_Coords (R1, R2);
  return Norm3d (&a);
}


struct Coord Wrap_Coordinates_PeriodicBoundary (struct Coord *in, struct Coord *ncell)
{
  struct Coord R;
  float q;

  R.x = (q=fmodf (in->x, ncell->x)) < 0 ? ncell->x + q : q;  R.x = fabs(R.x-ncell->x) < COORD_TOL ? 0.0 : R.x;
  R.y = (q=fmodf (in->y, ncell->y)) < 0 ? ncell->y + q : q;  R.y = fabs(R.y-ncell->y) < COORD_TOL ? 0.0 : R.y;
  R.z = (q=fmodf (in->z, ncell->z)) < 0 ? ncell->z + q : q;  R.z = fabs(R.z-ncell->z) < COORD_TOL ? 0.0 : R.z;

  return R;   
}


void Copy_Cluster (struct Cluster *dest, struct Cluster *in, int Nk)
{
  int j;

  for (j=0; j<Nk; j++) {
    dest->R[j].x = in->R[j].x;
    dest->R[j].y = in->R[j].y;
    dest->R[j].z = in->R[j].z;
  }

  return;
}


void Displace_Cluster (struct Coord *inR, struct Cluster *incls, struct Cluster *outcls, int Nk)
{
  int i;

  for (i=0; i<Nk; i++) outcls->R[i] = Add_Coords (&(incls->R[i]), inR);
  return;
}


void Displace_Cluster_PeriodicBoundary (struct Coord *inR, struct Cluster *incls, struct Cluster *outcls, int Nk, struct Coord *ncell)
{
  int i;
  struct Coord R;

  for (i=0; i<Nk; i++) {
    R = Add_Coords (&(incls->R[i]), inR);
    outcls->R[i] = Wrap_Coordinates_PeriodicBoundary (&R, ncell);
  }

  return;
}


void Get_Specific_Heat (double *Elist, int nlist, float beta, double *cp)
{
  int i;
  double Eavg, E2avg;

  if (nlist == 0) {(*cp)=0.0; return;}

  for (i=0, Eavg=0.0; i<nlist; i++)  Eavg+=Elist[i];             Eavg/=nlist;
  for (i=0, E2avg=0.0; i<nlist; i++) E2avg+=Elist[i]*Elist[i];   E2avg/=nlist;

  (*cp) = (E2avg-Eavg*Eavg)*(beta)*(beta)*kBOLTZMANN;
  return;
}
