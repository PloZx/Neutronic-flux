## Voici quelques indications importantes pour modifier certains paramètres du projet
  - Pour modifier la longueur du réacteur : modifier la variable "L" dans main.c & prob.c
  - Pour modifier la discrétisation spatiale : modifier la variable "m" dans main.c en respectant l'équation pour les m tels que mod8 = 1 (e.g:9,17,81,505)
## m = 505
![Img](projetfinalm505.png)

## Licenses
Ce projet utilise plusieurs bibliothèques numériques externes. Les informations de licence suivantes sont fournies afin d'assurer la conformité et d'attribuer correctement chaque composant.

### JADAMILU
JADAMILU est un solveur distribué sous licence GNU LGPL (Lesser General Public License).
Répertoire officiel de JADAMILU : [JADAMILU](https://github.com/mbollhoefer/JADAMILU)
### MC64 & MC21
Ce projet inclut des fichiers objets compilés provenant de MC64 & MC21 (voir jadamilu-main/lib/INT64YGNU), une routine de la HSL Mathematical Software Library, développée par le Science and Technology Facilities Council (STFC). <br>
Informations & licenses officielles :
[HSL](https://www.hsl.rl.ac.uk/)
### PRIMME
PRIMME (PReconditioned Iterative MultiMethod Eigensolver) est distribué sous la licence GNU LGPL.<br>
Répertoire github officiel : [PRIMME](https://github.com/primme/primme)


### Remarques pour la dernière question
Il faut considérer Jadamilu. J'ai d'abord mis Arpack puis j'ai décidé de changer pour la raison expliquée en commentaire dans arpack.c et je me suis dit autant laisser le code là au lieu de supprimer.

