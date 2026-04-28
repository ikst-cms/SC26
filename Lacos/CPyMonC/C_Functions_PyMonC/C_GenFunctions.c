//------------------------------------------------------------------------------------------------ 
// These are set of miscellaneous functions which are routinely required in programming.
//------------------------------------------------------------------------------------------------


char *alloc_1Dchar_array (int);             
int *alloc_1Dint_array (int);               
unsigned int *alloc_1Duint_array (int);     
long int *alloc_1Dlint_array (int);         

double *alloc_1Ddouble_array (int);         
float *alloc_1Dfloat_array (int);           

double **alloc_2Ddouble_array (int, int);   
void free_2Ddouble_array (double **, int);
float **alloc_2Dfloat_array (int, int);     
void free_2Dfloat_array (float **, int);    
int **alloc_2Dint_array (int, int);         
void free_2Dint_array (int **, int);        
int ***alloc_3Dint_array(int, int, int);  
void free_3Dint_array (int ***, int, int);
void digestchar (FILE *);                   
void digestline (FILE *);                    
int find_line (char *, FILE *, char *);     
void find_greater (float *, float *);       
void find_lesser (float *, float *);        
int min_of_array (float *, int);            
void double_swap (double *, double *);
void float_swap (float *, float *);
void int_swap (int *, int *);
void char_swap (char *, char *);
void long_int_swap (long int *, long int *);
int long_int_compare (long int *, long int *);
int exist_in_list (int *, int, int);             
int fexist_in_list (float *, int, float *);      
void print_header (char *, FILE *);
void print_float_vector (float *, int, FILE *); 
void print_double_vector (double *, int, FILE *); 
void print_float_matrix (float **, int, int, FILE *); 
void print_double_matrix (double **, int, int, FILE *); 
void Exit_Error (char *);
void Get_Simple_Statistics_float (float *, int, float *, float *, float *);
void Get_Simple_Statistics_double (double *, int, double *, double *, double *);
void hsortf (int , float []);
void siftdownf (float [], int, int);
void hsortd (int , double []);
void siftdownd (double [], int, int);


char *alloc_1Dchar_array (int nrow)
{
  char *z=NULL;

  if (nrow > 0) 
   if (!(z = (char *) malloc (nrow*sizeof(char)))) Exit_Error ("Unable to allocate 1D character array."); 
  
  return (z);
}


int *alloc_1Dint_array (int nrow)
{
  int *z=NULL;
  
  if (nrow > 0) 
    if (!(z = (int *) malloc (nrow*sizeof(int)))) Exit_Error ("Unable to allocate 1D int array."); 
 
  return (z);
}


unsigned int *alloc_1Duint_array (int nrow)
{
  unsigned int *z=NULL;
  
  if (nrow > 0) 
    if (!(z = (unsigned int *) malloc (nrow*sizeof(unsigned int)))) Exit_Error ("Unable to allocate 1D unsigned int array."); 

  return (z);
}


long int *alloc_1Dlint_array (int nrow)
{
  long int *z=NULL;
  
  if (nrow > 0) 
    if (!(z = (long int *) malloc (nrow*sizeof(long int)))) Exit_Error ("Unable to allocate 1D long int array."); 
  return (z);
}


double *alloc_1Ddouble_array (int nrow)
{
  double *z=NULL;

  if (nrow > 0) 
    if (!(z = (double *) malloc (nrow*sizeof(double)))) Exit_Error ("Unable to allocate 1D double array."); 
     
  return (z);
}


float *alloc_1Dfloat_array (int nrow)
{
  float *z=NULL;

  if (nrow > 0) 
    if (!(z = (float *) malloc (nrow*sizeof(float)))) Exit_Error  ("Unable to allocate 1D float array."); 
    
  return (z);
}


double **alloc_2Ddouble_array (int nrow, int ncol)
{
  int i;
  double **z=NULL;

  if (!(z = (double **) malloc (nrow*sizeof(double *)))) Exit_Error  ("Unable to allocate 2D double array."); 
 
  for (i=0; i<nrow; i++) z[i] = (double *) malloc (ncol*sizeof(double));
  return (z);
}


void free_2Ddouble_array (double **a, int nrow)
{
  int i;

  for (i=0; i<nrow; i++) free (a[i]);
  free (a); 
}


float **alloc_2Dfloat_array (int nrow, int ncol)
{
  int i;
  float **z=NULL;

  if (!(z = (float **) malloc (nrow*sizeof(float *)))) Exit_Error ("Unable to allocate 2D float array."); 
    
  for (i=0; i<nrow; i++) z[i] = (float *) malloc (ncol*sizeof(float));
  return (z);
}


void free_2Dfloat_array (float **a, int nrow)
{
  int i;
  
  for (i=0; i<nrow; i++) free (a[i]);
  free (a); 
}


int **alloc_2Dint_array (int nrow, int ncol)
{
  int i;
  int **z=NULL;

  if (!(z = (int **) malloc (nrow*sizeof(int *)))) Exit_Error ("Unable to allocate 2D int array."); 
  for (i=0; i<nrow; i++) z[i] = (int *) malloc (ncol*sizeof(int));
  return (z);
}


void free_2Dint_array (int **a, int nrow)
{
  int i;
  for (i=0; i<nrow; i++) free (a[i]);
  free (a); 
}


int ***alloc_3Dint_array(int nx, int ny, int nz)       
{
  int i, j;
  int ***z=NULL;

  z=(int ***) malloc (nx*sizeof(int **));
  for (i=0; i<nx; i++) z[i]=(int **) malloc (ny*sizeof(int *));
  for (i=0; i<nx; i++) for (j=0; j<ny; j++) z[i][j]=(int *) malloc (nz*sizeof(int));

  return (z);;
}


void free_3Dint_array (int ***z, int nx, int ny)
{
  int i, j;

  for (i=0; i<nx; i++) for (j=0; j<ny; j++) free (z[i][j]);
  for (i=0; i<nx; i++) free (z[i]);
  free (z); 

  return;
}




void digestchar (FILE *fp)          
{
  int k;
  k=fscanf (fp, "%*[^\n]\n");
  if (k == EOF) printf ("\nError in reading file (in digestchar); reached EOF.\n");    
}


void digestline (FILE *fp)          
{
  int k;
  k=fscanf (fp, "%*[^\n]\n");
  if (k == EOF) printf ("\nError in reading file (in digestline); reached EOF.\n");    
}


int find_line (char *in, FILE *fpx, char *fname)
{
  char t[200];
  char *r=NULL;

  rewind (fpx);
  while ( fgets (t, sizeof(t), fpx) != NULL) {
    if ((r=strstr(t, in)) != NULL) return (1);
  }
  printf ("Unable to find string: %s in file %s\nExit\n", in, fname);
  exit (EXIT_FAILURE);
}


void find_greater (float *p, float *q)
{
  if ((*p)>(*q)) {return;}
  else {(*p)=(*q); return;}
}


void find_lesser (float *p, float *q)
{
  if ((*p)<(*q)) {return;}
  else {(*p)=(*q); return;}
}


int min_of_array (float *a, int asize)
{
  int i, k;

  for (i=1, k=0; i<asize; i++) if (a[i]<a[k]) k=i;

  return (k);
}


void float_swap (float *a, float *b)
{
  float temp;
  temp=*a; *a=*b; *b=temp;
}


void int_swap (int *a, int *b)
{
  int temp;
  temp=*a; *a=*b; *b=temp;
}


void long_int_swap (long int *a, long int *b)
{
  long int temp;
  temp=*a; *a=*b; *b=temp;
}


void double_swap (double *a, double *b)
{
  double temp;
  temp=*a; *a=*b; *b=temp;
}


int long_int_compare (long int *a, long int *b)
{
  if (*a < *b) return -1;
  if (*a > *b) return 1;
  return 0;
} 


int exist_in_list (int *inlist, int list_len, int k)
{
  int i;

  for (i=0; i<list_len; i++) if (inlist[i] == k) return (i);
  return (-1);
}



int fexist_in_list (float *inlist, int list_len, float *k)
{
  int i;

  for (i=0; i<list_len; i++) if (fabs(inlist[i]-(*k))<1.0e-5) return (i);
  return (-1);
}



void print_header (char *str, FILE *fin)
{
  int i;

  for (i=0; i<(60+strlen(str)); i++) fprintf (fin, "-"); fprintf (fin, "\n");
  fprintf (fin, "                              %s                              \n", str);    
  for (i=0; i<(60+strlen(str)); i++) fprintf (fin, "-"); fprintf (fin, "\n");
  return;
}


void print_float_matrix (float **a, int nrow, int ncol, FILE *fin)
{
  int i, j;

  if (fin == NULL) {
    printf ("\n");
    for (i=0; i<nrow; i++) {
      for (j=0; j<ncol; j++) printf ("%7.4f ", a[i][j]);
      printf ("\n");
    }}
  else {
    fprintf (fin, "\n");
    for (i=0; i<nrow; i++) {
      for (j=0; j<ncol; j++) fprintf (fin, "%7.4f ", a[i][j]);
      fprintf (fin, "\n");
    }
  }

  return;
}


void print_double_matrix (double **a, int nrow, int ncol, FILE *fin)
{
  int i, j;

  if (fin == NULL) {
    printf ("\n");
    for (i=0; i<nrow; i++) {
      for (j=0; j<ncol; j++) printf ("%7.4f ", a[i][j]);
      printf ("\n");
    }}
  else {
    fprintf (fin, "\n");
    for (i=0; i<nrow; i++) {
      for (j=0; j<ncol; j++) fprintf (fin, "%7.4f ", a[i][j]);
      fprintf (fin, "\n");
    }
  }

  return;
}


void print_double_vector (double *a, int n, FILE *fin)
{
  int j;

  if (fin == NULL) {
    printf ("\n");
    for (j=0; j<n; j++) printf ("%7.4f ", a[j]); printf ("\n"); }
  else {
    fprintf (fin, "\n");
    for (j=0; j<n; j++) fprintf (fin, "%7.4f ", a[j]); fprintf (fin, "\n"); 
  }

  return;
}


void print_float_vector (float *a, int n, FILE *fin)
{
  int j;

  if (fin == NULL) {
    printf ("\n");
    for (j=0; j<n; j++) printf ("%7.4f ", a[j]); printf ("\n"); }
  else {
    fprintf (fin, "\n");
    for (j=0; j<n; j++) fprintf (fin, "%7.4f ", a[j]); fprintf (fin, "\n"); 
  }

  return;
}



void Exit_Error (char *p)
{
  printf ("\n\n >>>>>>>  Error :: %s   Exit.  <<<<<<<<\n", p);
  exit (EXIT_FAILURE);
}



void Get_Simple_Statistics_float (float *list, int nlist, float *mean, float *median, float *std)
{
  int i;
  float tmp, atmp;

  if (nlist == 0) { (*mean)=0.0; (*median)=0.0; (*std)=0.0; return;}

  for (i=0, atmp=0.0; i<nlist; i++) atmp+=list[i];
  (*mean)=(atmp/=nlist);
  
  for (i=0, tmp=0.0; i<nlist; i++) tmp+=powf(list[i]-atmp, 2.0);
  (*std)=sqrtf(tmp/nlist);
  
  if (nlist < 3) (*median) = (*mean);
  else {
    hsortf (nlist, list);
    if (nlist%2 == 1) (*median) = list[nlist/2];
    else (*median) = 0.5*(list[nlist/2]+list[(nlist/2)-1]);
  }

  return;
}



void Get_Simple_Statistics_double (double *list, int nlist, double *mean, double *median, double *std)
{
  int i;
  double tmp, atmp;

  if (nlist == 0) { (*mean)=0.0; (*median)=0.0; (*std)=0.0; return;}
  
  for (i=0, atmp=0.0; i<nlist; i++) atmp+=list[i];
  (*mean)=(atmp/=nlist);
  
  for (i=0, tmp=0.0; i<nlist; i++) tmp+=pow(list[i]-atmp, 2.0);
  (*std)=sqrt(tmp/nlist);

  if (nlist < 3) (*median) = (*mean);
  else {
    hsortd (nlist, list);
    if (nlist%2 == 1) (*median) = list[nlist/2];
    else (*median) = 0.5*(list[nlist/2]+list[(nlist/2)-1]);
  }

  return;
}



void hsortf (int n, float ra[])
{
  int start, end;

  for (start=(n-2)/2; start>=0; start--) siftdownf (ra, start, n);
  for (end=n-1; end>0; end--) {float_swap (&ra[0], &ra[end]); siftdownf (ra, 0, end);}
  return;
}

void siftdownf (float ra[], int start, int end)
{
  int root=start, child;

  while (2*root+1 < end) {
    child=2*root+1;
    if ((child+1 < end) && (ra[child]<ra[child+1])) child+=1;
    if(ra[root]<ra[child]) {float_swap (&ra[root], &ra[child]); root=child;} else return;
  }
}


void hsortd (int n, double ra[])
{
  int start, end;

  for (start=(n-2)/2; start>=0; start--) siftdownd (ra, start, n);
  for (end=n-1; end>0; end--) {double_swap (&ra[0], &ra[end]); siftdownd (ra, 0, end);}
  return;
}

void siftdownd (double ra[], int start, int end)
{
  int root=start, child;

  while (2*root+1 < end) {
    child=2*root+1;
    if ((child+1 < end) && (ra[child]<ra[child+1])) child+=1;
    if(ra[root]<ra[child]) {double_swap (&ra[root], &ra[child]); root=child;} else return;
  }
}





