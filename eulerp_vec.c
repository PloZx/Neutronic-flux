#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "norme_residu.h"
#include "eulerp_vec.h"

/* Euler progressive : phi(t+dt) = phi(t) - dt * tau * [(A - beta2 I) phi(t)]
   - n, ia, ja, a : matrice A en format CSR
   - beta2, tau   : paramètres
   - dt           : pas de temps, choisi avec le théorème de Greshgorin
   - nsteps       : nombre de pas
   - tol          : tolérance sur la norme du résidu
   - phi          : vecteur état (IN/OUT) ; doit contenir l’état initial à l’appel donc il faut le faire dans main.c
   - log_every    : sampling
   - pipe         : pipe vers Python
*/
void eulerp_vec(int n,
                int *ia, int *ja, double *a,
                double beta2, double tau,
                double dt, int nsteps,
                double tol,
                double *phi,
                int log_every,
                FILE *pipe)
{
	/* initialisation vecteur */
    double *y = (double*)malloc((size_t)n * sizeof(double));
    if (!y) {
        fprintf(stderr, "eulerp_vec: OOM\n");
        return;
    }

    /* Ceci sert à faire de légères optimisations en garantissant au compilateur qu'il n'y aura pas d'aliasing "restrict"*/
    int * restrict ia_r  = ia;
    int * restrict ja_r  = ja;
    double * restrict a_r = a;
    double * restrict phi_r = phi;
    double * restrict y_r = y;
    double c = dt * tau;


    /*Envoi initial à Python*/
        if (pipe) {
        	fprintf(pipe, "# step %d t=%e\n", 0, 0.0);
        	for (int i=0; i<n; ++i)
        		fprintf(pipe, "%.17g ", phi_r[i]); // À ce stade, on a la condition initale d'envoyée
        	fprintf(pipe, "\n"); fflush(pipe);
        }

    /* Boucle Euler prograssive vectorielle */
    for (int k = 1; k <= nsteps; ++k) {
        /* y = A * phi */
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int p = ia_r[i]; p < ia_r[i+1]; ++p)
                s += a_r[p] * phi_r[ ja_r[p] ];
            y_r[i] = s;
            /* phi <- phi - dt * tau * (y - beta2*phi) */
            double ri = s - beta2 * phi_r[i];
            phi_r[i] -= c * ri;

        }

        residu_stats s = norme_residu(n, (int*)ia_r, (int*)ja_r, (double*)a_r, beta2, phi_r);

        /* Envoie des données vers python avec un sampling sur les vecteurs */
        if (log_every > 0 && pipe && (k % log_every == 0)) {
        	fprintf(pipe, "# step %d t=%e\n", k, k*dt);
        	for (int i=0; i<n; ++i)
        		fprintf(pipe, "%.17g ", phi_r[i]);
        	fprintf(pipe, "\n");
        	fflush(pipe);
        }

        /* Au cas où, en pratique n'est jamais encore arrivé. Il est possible de retirer ce if ainsi que l'appel à norme_residu pour un code plus effiace. Je ne le fais pas ici */
        if (s.norme < tol) {
           fprintf(stderr, "[euler] convergence atteinte au pas %d\n", k);
           break;
        }
    }

    free(y);
}
