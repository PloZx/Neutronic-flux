#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interface_primme.h"
#include "time.h"
#include "norme_residu.h"
#include "primme.h"

/* Calcul de la norme du résidu
 *
 * 	- n : taille de la matrice
 * 	- ia, ja et a : matrice en format CSR
 * 	- beta2 : valeur propre
 * 	- phi : flux
 *
 * Remarque : Au départ j'ai simplement utilisé matvec de primme mais je me suis rendu compte qu'on pouvait optimiser la chose suivante :
 * dans le calcul du produit matrice-vecteur on peut réutiliser le "sum" directement après au lieu d'aller le relire plus tard.
 */
residu_stats norme_residu(int n,int *ia,int *ja, double *a, double beta2, double *phi){
	double num2 = 0.0;
	double den2 = 0.0;
	/* Matvec*/
	for (int i = 0; i < n; ++i) {
	    double sum = 0.0;
	    for (int p = ia[i]; p < ia[i+1]; ++p)
	        sum += a[p] * phi[ ja[p] ];

	    double phii = phi[i];
	    double ri = sum - beta2 * phii;

	    num2 += ri * ri;
	    den2 += phii * phii;
	}
	/* Structure */
    residu_stats out;
    out.num = sqrt(num2);
    out.den = sqrt(den2);
    out.norme = (den2 > 0.0) ? sqrt(num2) / sqrt(den2) : 0.0;
    return out;
}
