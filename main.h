/* prototypes utilisés dans main.h */
double mytimer_cpu();
double mytimer_wall();

int prob(int m, int *n, int **ia, int **ja, double **a);
int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs);

extern void dpjd_(int *N,double *A,int *JA,int *IA,double *EIGS,double *RES,double *X,int   *LX,int   *NEIG,double *SIGMA,int  *ISEARCH,int  *NINIT,int  *MADSPACE,int  *ITER,double *TOL,double *SHIFT,double *DROPTOL,double *MEM,int  *ICNTL, int  *IPRINT,int  *INFO, double *GAP);

