#include <stdlib.h>
#include <stdint.h>

typedef long fint;

/*
 * Convertit une matrice CSR C 0-based (n, ia, ja, a)
 * en CSR Fortran 1-based (IA_F, JA_F, A_F) mais en version triangulaire supérieur pour Jadamilu,
 *
 *   n        : taille de la matrice
 *   ia, ja, a: CSR 0-based
 *   *IA_F_out : tableau de taille n+1, de type fint, 1-based fortran
 *   *JA_F_out : tableau de taille nnz_upper, de type fint, colonnes 1-based
 *   *A_F_out  : tableau de taille nnz_upper, double
 *   *nnz_F_out: nombre d’entrées gardées (triang sup)
 *
 *   Remarque : en réalité, on fait comme dans prob.c mais il faut juste retirer certaines parties.
 */
int csrC_to_upperCSRFortran(int n,int *ia,int *ja,double *a,fint **IA_F_out,fint **JA_F_out,double **A_F_out,fint *nnz_F_out) {
    // compter les nnz dans la partie triangulaire sup
    fint nnz_upper = 0;
    for (int i = 0; i < n; ++i) {
        for (int k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) {
                nnz_upper++;
            }
        }
    }

    // Allocation
    fint *IA_F = (fint*)malloc((n + 1) * sizeof(fint));
    fint *JA_F = (fint*)malloc((size_t)nnz_upper * sizeof(fint));
    double *A_F  = (double*)malloc((size_t)nnz_upper * sizeof(double));
    if (!IA_F || !JA_F || !A_F) {
        free(IA_F); free(JA_F); free(A_F);
        return -1;
    }

    // remplir IA_F (version 0-based)
    IA_F[0] = 0;
    fint count = 0;
    for (int i = 0; i < n; ++i) {
        for (int k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) {
                count++;
            }
        }
        IA_F[i+1] = count;
    }

    // remplir JA_F et A_F (toujours en 0-based)
    for (int i = 0; i < n; ++i) {
        fint write_pos = IA_F[i];  // début de la ligne i en 0-based
        for (int k = ia[i]; k < ia[i+1]; ++k) {
            int j = ja[k];
            if (j >= i) {
                JA_F[write_pos] = (fint)j;   // 0-based
                A_F[write_pos]  = a[k];
                write_pos++;
            }
        }
    }

    //  conversion en indices 1-based
    for (int i = 0; i <= n; ++i) {
        IA_F[i] = IA_F[i] + 1;
    }
    for (fint k = 0; k < nnz_upper; ++k) {
        JA_F[k] = JA_F[k] + 1;
    }

    *IA_F_out  = IA_F;
    *JA_F_out  = JA_F;
    *A_F_out   = A_F;
    *nnz_F_out = nnz_upper;

    return 0;
}
