/* prototypes utilisés dans main.h */
double mytimer_cpu();
double mytimer_wall();

int prob(int m, int *n, int **ia, int **ja, double **a);
int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs);

typedef long long fint;
extern void dpjd_(int *N,double *A,fint *JA,fint *IA,double *EIGS,double *RES,double *X,fint   *LX,fint   *NEIG,double *SIGMA,fint  *ISEARCH,fint  *NINIT,fint  *MADSPACE,fint  *ITER,double *TOL,double *SHIFT,double *DROPTOL,double *MEM,fint  *ICNTL, fint  *IPRINT,fint  *INFO, double *GAP);
static double maxd(const double *v, int len);
static double mind(const double *v, int len);

