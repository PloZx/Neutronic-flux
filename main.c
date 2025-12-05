#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "print_csr.h"
#include "norme_residu.h"
#include "eulerp_vec.h"
#include <math.h>
#include "arpack.h"
#include "primme.h"
#include <stdint.h>
#include "CtoFortran.h"

typedef long long fint;


int main(int argc, char *argv[])
{
  /* déclarer les variables */

  int m = 505, nev = 1;
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs, *resNorms, *phi_t, *evals_ar, *evecs_ar;
  double tc1, tc2, tw1, tw2;
  double alpha = sqrt(5.85/2.00);
  double L = 4.0;
  double h = (L / (m - 1)) * alpha;

  /* générer le problème */
  if (prob(m, &n, &ia, &ja, &a))
     return 1;

  printf("\nPROBLÈME: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  
  /* allouer la mémoire pour vecteurs & valeurs propres de PRIMME et ARPACK */
  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));
  resNorms = malloc(nev * sizeof(double));
  phi_t = malloc(nev * n * sizeof(double));
  evals_ar = malloc(nev * sizeof(double));
  evecs_ar = malloc(n * nev * sizeof(double));

  /* Condition initiale pour le problème temporel */
  for (int i=0;i<n;i++) phi_t[i] = 3.0;


  if (evals == NULL || evecs == NULL || resNorms == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour les vecteurs, valeurs propres et normes des résidus\n\n");
      return 1;
  }

  /* primme - résolution */
  tc1 = mytimer_cpu(); tw1 = mytimer_wall();
  if(primme(n, ia, ja, a, nev, evals, evecs, resNorms,primme_smallest))
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall();

  /* ARPACK - résolution, désactivée pour les cas où m est trop grand car trop lent*/
  //arpack_sym_csr_eigs(n, ia, ja, a, nev, "SA", 10000, 0, evals_ar, evecs_ar);

  /* Conversion -> triang sup Fortran */

  fint *IA_J = NULL;
  fint *JA_J = NULL;
  double *A_J = NULL;
  fint NNZ_F = 0;
  if (csrC_to_upperCSRFortran(n, ia, ja, a,&IA_J, &JA_J, &A_J, &NNZ_F) != 0) {
	  fprintf(stderr, "Erreur conversion CSR -> Fortran\n");
      return 1;
}
  fint N = (fint)n;
  fint NEIG = 1;          // on veut 1 valeur propre
  fint MAXEIG = NEIG;     // comme dans l’exemple Fortran
  fint MAXSP  = 20;       // espace de Krylov

  // LX = N*(3*MAXSP + MAXEIG + 1) + 4*MAXSP*MAXSP
  fint LX = N * (3*MAXSP + MAXEIG + 1) + 4*MAXSP*MAXSP;

  double *EIGS_J = malloc((size_t)MAXEIG * sizeof(double));
  double *RES_J  = malloc((size_t)MAXEIG * sizeof(double));
  double *X    = malloc((size_t)LX     * sizeof(double));


  if (!EIGS_J || !RES_J || !X) {
       fprintf(stderr, "Erreur alloc EIGS/RES/X\n");
       return 1;
  }

  // paramètres de contrôle (d’après la doc JADAMILU)
  double SIGMA   = 2.0;
  fint   ISEARCH = 1; // REMARQUE IMPORTANTE : Ici j'utilise le fait que j'ai une idée de la valeur propre la plus petite (à peu près 2) si on s'éloigne pas trop non plus des dimensions stationnaires du réacteur.
  fint   NINIT   = 0;
  fint   MADSPACE = MAXSP;
  fint   ITER    = 1000;
  double TOL     = 1e-6;
  double SHIFT   = 2.0;
  double DROPTOL = 1e-3;
  double MEM     = 20.0;

  fint ICNTL[5] = {0,0,1,0,0};
  fint IPRINT   = 6;
  fint INFO     = 0;
  double GAP    = 0.0;

  /* Appel JADAMILU*/
  dpjd_(&N,A_J,JA_J,IA_J,EIGS_J,RES_J,X,&LX,&NEIG,&SIGMA,&ISEARCH,&NINIT,&MADSPACE,&ITER,&TOL,&SHIFT,&DROPTOL,&MEM,ICNTL,&IPRINT,&INFO,&GAP);


  // Visualisations avec matplotlib (python)
  // Création de pipe pour communiquer m et les vecteurs phi à python pour les visualisations
  char cmd[512];
  snprintf(cmd,sizeof(cmd),"python3 plotfinal.py --L %f --m %d",L,m); // ceci adapte la commande executée en fonction de la valeur de m et L
  FILE *pipe = popen(cmd, "w");
  // Cas stationnaire
  for (int i = 0; i < n; i++)
	  fprintf(pipe, "%.17g\n", evecs[i]);
  fprintf(pipe, "\n");
  fflush(pipe);
  // Cas temporel, pour m grand on ne le fait pas (choix personnel)
  if(m<=141)
	  eulerp_vec(n, ia, ja, a, 2.0, 100, (h*h)/(4*100), 3000, 1e-12, phi_t, 100,pipe);
  else
	  fprintf(stderr, "Evolution temporelle désactivée dû à m=%d trop grand",m);

  pclose(pipe);
  /* Calcul de la norme du résidu et affichage dans la console*/
  residu_stats nres = norme_residu(n, ia, ja, a, evals[0], evecs); //norme_residu renvoie une structure
  printf("\nNorme du résidu calculée : %5.10e\n",nres.norme);
  printf("\nNorme du résidu PRIMME: %5.10e\n",resNorms[0]);

  /* temps de solution */
  printf("\nTemps de solution (CPU): %5.1f sec",tc2-tc1);
  printf("\nTemps de solution (horloge): %5.1f sec \n",tw2-tw1);
  printf("\nValeur propre minimale calculée: %5.5f\n",evals[0]);
  printf("Valeurs propre minimale ARPACK calculée: %5.5f\n", evals_ar[0]);


  /* libérer la mémoire */
  free(ia); free(ja); free(a); free(evals); free(evecs), free(resNorms), free(phi_t), free(evals_ar); free(evecs_ar);free(IA_J); free(JA_J); free(A_J);free(EIGS_J); free(RES_J); free(X);
  return 0;
}

