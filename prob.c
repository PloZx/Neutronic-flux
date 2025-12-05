#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int prob(int m, int *n, int **ia, int **ja, double **a)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la discrétisation sur une grille 
   cartésienne régulière m x m de l'opérateur de Laplace à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )        sur [0,2] x [0,2]
           dx   dx       dy   dy

  avec la fonction u qui satisfait les conditions aux limites de Dirichlet
         
         u = 0  sur (0,y), (2,y), (x,0) et (x,2), avec 0 <= x,y <= 2 .
  
  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice est retournée dans le format CRS
  qui est défini par le scalaire 'n' et les trois tableaux 'ia, 'ja' et 'a'.

  Arguments
  =========
  m (input)   - nombre de points par direction dans la grille 
  n  (output) - pointeur vers le nombre d'inconnues dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
 
  Sortie
  ======
  0 - exécution avec succès
  1 - erreurs
*/
{
    int  nnz = 0, ix, iy, nx, ind, compteur = 0;
    double invh2, h = 0, alpha = sqrt(5.85/2.00), L = 4.0;
    nx = m - 2; /* nœuds de Dirichlet ne sont pas pris en compte */
    h = (L/(m-1))*alpha; /* Pas de discretisation */
    invh2 = (m-1)*(m-1)/(L*L*alpha*alpha); /* h^-2 pour L=2 */

    /* Valeur -1 accordée aux points ne faisant pas parti du domaine, les autres sont normalement numérotés */
    int *map = malloc(nx*nx * sizeof(int)); // Allocation dans la mémoire

    for (iy= 0; iy<nx; iy++) {
        for (ix=0; ix<nx; ix++) {
            double x = (ix+1) * h;
            double y = (iy+1) * h;
            if (x >= (1.5/2)*L*alpha && x <= L*alpha  && y >= 0.0 && y <= (0.75/2)*L*alpha) {
                map[ix+nx*iy] = -1; // On marque les éléments hors du domaine

            } else {
                map[ix+nx*iy] = compteur++;
          }
       }
    }
    *n = compteur; // nombre réel d'inconnues

    /* Compte le nombre d'éléments non nuls réel */
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            int ind = ix + nx * iy;
            int row = map[ind];
            if (row < 0) continue;   /* nœud exclus skippés du domaine puisqu'ils ont la valeur -1*/

            /* diagonale toujours présente */
            nnz += 1;

            /* voisins valides uniquement */
            if (iy > 0     && map[ind - nx] >= 0) nnz++;   /* sud */
            if (ix > 0     && map[ind - 1]  >= 0) nnz++;   /* ouest */
            if (ix < nx-1  && map[ind + 1]  >= 0) nnz++;   /* est */
            if (iy < nx-1  && map[ind + nx] >= 0) nnz++;   /* nord */
        }
    }

    /* allocation des tableaux */

    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));

    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */

    nnz = 0;
    for (iy = 0; iy < nx; iy++) {
        for (ix = 0; ix < nx; ix++) {
            /* index de la grille */
            ind = ix + nx * iy;
            int row = map[ind];

            /* ligne de la matrice (=-1 si nœud exclu) */
            row = map[ind];
            if (row < 0) continue;  /* nœud retiré du domaine si indice négatif*/

            /* marquer le début de la ligne 'row' */
            (*ia)[row] = nnz;

            /* voisin sud */
            if (iy > 0) {
                int col = map[ind - nx];
                if (col >= 0) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = col;
                    nnz++;
                }
            }

            /* voisin ouest */
            if (ix > 0) {
                int col = map[ind - 1];
                if (col >= 0) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = col;
                    nnz++;
                }
            }

            /* diagonale */
            (*a)[nnz]  = 4.0 * invh2;
            (*ja)[nnz] = row;
            nnz++;

            /* voisin est */
            if (ix < nx - 1) {
                int col = map[ind + 1];
                if (col >= 0) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = col;
                    nnz++;
                }
            }

            /* voisin nord */
            if (iy < nx - 1) {
                int col = map[ind + nx];
                if (col >= 0) {
                    (*a)[nnz]  = -invh2;
                    (*ja)[nnz] = col;
                    nnz++;
                }
            }
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[*n] = nnz;

    /* retour habituel de fonction */
    free(map);
    return 0;
}
