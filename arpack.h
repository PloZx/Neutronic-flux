int arpack_sym_csr_eigs(
    int n, const int *ia, const int *ja, const double *a,
    int nev, const char *which, int maxit, double tol,
    double *evals, double *evecs
);

extern void dsaupd_(
    int *ido, char *bmat, int *n, char *which,
    int *nev, double *tol, double *resid,
    int *ncv, double *v, int *ldv, int *iparam,
    int *ipntr, double *workd, double *workl,
    int *lworkl, int *info
);

extern void dseupd_(
    int *rvec, char *howmny, int *select, double *d,
    double *z, int *ldz, double *sigma,
    char *bmat, int *n, char *which, int *nev,
    double *tol, double *resid, int *ncv, double *v,
    int *ldv, int *iparam, int *ipntr, double *workd,
    double *workl, int *lworkl, int *info
);
