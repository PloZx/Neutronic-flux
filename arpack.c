#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "arpack.h"


/* Produit matrice-vecteur CSR : y = A x
   - n  : taille
   - ia : taille n+1
   - ja : taille nnz
   - a  : taille nnz
*/
static void csr_matvec(int n,
                       const int *ia, const int *ja, const double *a,
                       const double *x, double *y)
{
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = ia[i]; k < ia[i+1]; ++k) {
            s += a[k] * x[ja[k]];
        }
        y[i] = s;
    }
}

/*
   - calcule nev valeurs propres + vecteurs propres d’une matrice symétrique A
     stockée en CSR.
   - nev valeurs propres sont retournées dans evals[0..nev-1]
   - les vecteurs propres sont colonnes de evecs (taille n*nev)
   Retour :
     0  : succès, nev valeurs propres convergées
     >0 : seulement "ret" valeurs propres convergées
     <0 : erreur (allocation ou erreur ARPACK sérieuse)


    Cette partie ne vient pas de moi-même. Ce code est basé sur des exemples trouvés sur le github de Arpack-ng et certains éléments ont été fournis par une IA (d'où le fait que j'ai implementé par la suite Jadamilu, par convictions personnelles).
*/
int arpack_sym_csr_eigs(
    int n, const int *ia, const int *ja, const double *a,
    int nev, const char *which, int maxit, double tol,
    double *evals, double *evecs
) {
    if (n <= 0 || nev <= 0 || nev >= n) {
        fprintf(stderr,
                "arpack_sym_csr_eigs: arguments invalides (n=%d, nev=%d)\n",
                n, nev);
        return -1;
    }

    /* Dimension du sous-espace de Krylov */
    int ncv = 2*nev + 20;
    if (ncv > n)      ncv = n;
    if (ncv < nev+2)  ncv = nev + 2;

    int ldv    = n;
    int ldz    = n;
    int lworkl = ncv*(ncv + 8);   /* doc ARPACK : réel sym. */

    double *resid = (double*) malloc(n * sizeof(double));
    double *V     = (double*) malloc(ldv * ncv * sizeof(double));
    double *Z     = evecs;        /* vecteurs propres de sortie */
    double *D     = evals;        /* valeurs propres de sortie */
    double *workd = (double*) malloc(3 * n * sizeof(double));
    double *workl = (double*) malloc(lworkl * sizeof(double));
    int    *select= (int*)    malloc(ncv * sizeof(int));

    if (!resid || !V || !workd || !workl || !select) {
        fprintf(stderr, "arpack_sym_csr_eigs: erreur malloc\n");
        free(resid); free(V); free(workd); free(workl); free(select);
        return -2;
    }

    int iparam[11], ipntr[11];
    memset(iparam, 0, sizeof(iparam));
    memset(ipntr,  0, sizeof(ipntr));

    iparam[0] = 1;         /* ishift = 1 */
    iparam[2] = maxit;     /* maxit */
    iparam[3] = 1;         /* NB (block size) = 1 */
    iparam[6] = 1;         /* mode standard */

    int ido   = 0;
    int info  = 0;
    char bmat = 'I';

    char which_[3] = "SA";
    if (which && strlen(which) >= 2) {
        which_[0] = which[0];
        which_[1] = which[1];
        which_[2] = '\0';
    }

    if (tol < 0.0)
        tol = 0.0;         /* tol=0 -> tolérance machine (ARPACK) */

    /* Boucle de reverse communication */
    do {
        dsaupd_(&ido, &bmat, &n, which_, &nev, &tol,
                resid, &ncv, V, &ldv, iparam, ipntr,
                workd, workl, &lworkl, &info);

        if (ido == -1 || ido == 1) {
            /* x = workd[ipntr[0]-1], y = workd[ipntr[1]-1] */
            csr_matvec(n, ia, ja, a,
                       &workd[ipntr[0]-1],
                       &workd[ipntr[1]-1]);
        }
    } while (ido == -1 || ido == 1);

    if (info < 0) {
        fprintf(stderr, "arpack_sym_csr_eigs: dsaupd_ a échoué, info=%d\n", info);
        free(resid); free(V); free(workd); free(workl); free(select);
        return info;  /* négatif */
    }

    int nconv = iparam[4];   /* nombre de valeurs propres convergées */

    if (nconv <= 0) {
        fprintf(stderr,
                "arpack_sym_csr_eigs: aucune valeur propre n'a convergé (nconv=%d)\n",
                nconv);
        free(resid); free(V); free(workd); free(workl); free(select);
        return 0;   /* 0 convergées -> rien d’utile dans evals */
    }

    int rvec   = 1;         /* on veut les vecteurs propres */
    char howmny = 'A';
    double sigma = 0.0;
    int ierr = 0;

    dseupd_(&rvec, &howmny, select, D, Z, &ldz, &sigma,
            &bmat, &n, which_, &nev, &tol,
            resid, &ncv, V, &ldv, iparam, ipntr,
            workd, workl, &lworkl, &ierr);

    free(resid); free(V); free(workd); free(workl); free(select);

    if (ierr != 0) {
        fprintf(stderr, "arpack_sym_csr_eigs: dseupd_ a échoué, ierr=%d\n", ierr);
        return -10 + ierr;
    }

    if (nconv < nev) {
        fprintf(stderr,
                "arpack_sym_csr_eigs: seulement %d valeurs propres ont convergé sur %d.\n",
                nconv, nev);
        /* tu peux remplir le reste avec NaN si tu veux être sûr de ne pas lire 0 par erreur */
        for (int i = nconv; i < nev; ++i) {
            D[i] = NAN;
        }
        return nconv;     /* nombre de valeurs propres fiables */
    }

    return 0;  /* succès total */
}
