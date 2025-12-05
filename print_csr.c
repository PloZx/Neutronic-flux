#include <stdio.h>
/*
 *
 *
 * Ce code ne sert qu'à des fins de verifications pour la première question, il n'est pas à considérer. Je le laisse simplement comme trace.
 *
 *
 */
void print_csr(int n, const int *ia, const int *ja, const double *a) {
    int nnz = ia[n];

    printf("n = %d, nnz = %d\n", n, nnz);

    printf("ia (%d+1 = %d elems):\n", n, n+1);
    for (int i = 0; i <= n; i++) {
        printf("%d%s", ia[i], (i==n? "\n" : " "));
    }

    printf("ja (%d elems):\n", nnz);
    for (int k = 0; k < nnz; k++) {
        printf("%d%s", ja[k], (k==nnz-1? "\n" : " "));
    }

    printf("a (%d elems):\n", nnz);
    for (int k = 0; k < nnz; k++) {
        printf("%.15g%s", a[k], (k==nnz-1? "\n" : " "));
    }
}

