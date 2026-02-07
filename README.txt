ATTENTION : 
SI JAMAIS IL Y A UN PROBLEME AVEC LES NOTES DE CALCUL (ça m'est arrivé), RELANCER LE MAIN AVEC LE BON CAS_TEST PUIS LA NOTE DE CALCUL.
(PAR EXEMPLE SI JE VEUX VOIR LA NOTE DE CALCUL DU CAS TEST 5, JE RELANCE D'ABORD LE MAIN AVEC LE CAS_TEST5 PUIS LA NOTE DE CALCUL)

README BE Méthode des Déplacements 

Ce programme MATLAB permet l’analyse d’une structure plane en poutres 2D (Euler–Bernoulli) à l’aide de la méthode des déplacements.
Il calcule :

- les déplacements nodaux,
- les réactions d’appui,
- les efforts internes (N, V, M),
- la déformée
------------------------------------------------------------
0) AUTEUR 
------------------------------------------------------------

Matthieu MISKOWSKI - M1 MAISES - ENS Paris Saclay
BE Calcul des structures - Méthode des déplacements 

------------------------------------------------------------
1) CONTENU DU FICHIER
------------------------------------------------------------

main.m : script principal
CAS_TEST/                
   Cas_test_1.m : cas test de base du BE, 3 barres, 3 appuis, 4 nœuds 
   Cas_test_2.m : cas test sur une structure hyperstatique d'application, double portique contreventé, 14 barres, 3 appuis, 9 noeuds 
   Cas_test_3.m : cas test de base du BE avec un appui élastique, 3 barres, 3 appuis, 4 nœuds  
   Cas_test_4.m : cas test de base du BE avec un effort sur barre centré, 3 barres, 3 appuis, 4 nœuds 
   Cas_test_5.m : cas test de base du BE avec une barre multifibre, 3 barres, 3 appuis, 4 nœuds 

FONCTIONS/
   f_rigidite.m : calcul de la matrice de rigidité locale (6x6) d'un élément de poutre 
   f_passage.m : calcul de la matrice de passage (rotation) 6x6 entre le repère local et global.
   f_forme.m : calcul de la matrice des fonctions de forme d’un élément de poutre.
   f_charge_equivalente_point_centre.m : calcul du vecteur des efforts nodaux équivalents (6x1) dans le repère local
   f_section_multifibre.m : calcul EA_eq, EI_eq, yE pour une section multifibre discrétisée
   f_ampli.m : calcul automatique des coefficients d’amplification

Une note de calcul pour le cas_test_1 plus épaisse qui détaille comment le programme fonctionne
Puis des note de calculs pour les cas tests 2,3,4 et 5 plus succinctes qui présentent les résultats et les hypothèses supplémentaires 

Deux vidéos faites avec manim pour expliquer le passage local -> global et l'assemblage de la matrice K (la deuxième est un peu raté sur la fin)

------------------------------------------------------------
2) UTILISATION
------------------------------------------------------------

Dans main.m (4e ligne) modifier :

    Cas_test_1;   % en remplaçant le 1 par le n° du cas test voulu

Puis exécuter le programme

Le programme affiche automatiquement :
- la déformée,
- N(x),
- V(x),
- M(x).

------------------------------------------------------------
3) CAS MULTIFIBRE (Cas_test_5)
------------------------------------------------------------

L’élément n°2 peut être modélisé comme une section composite discrétisée en un nombre quelconque de fibres.

Il faut donner l'indice de l'élément qui est multifibre et ses paramètres géométriques dans Cas_test_5, 

"MULTIFIBRE : paramètres pour l’élément"

Le cas multifibre calcule :

EA_eq : rigidité axiale équivalente,

EI_eq : rigidité en flexion équivalente,

yE : position de l’axe neutre équivalent.

La fonction utilisée est :

[EA_eq, EI_eq, yE] = f_section_multifibre(E_f, b_f, h_f, y_f);


Elle prend directement :

les modules E_f par fibre,

les dimensions b_f, h_f,

la position de la fibre y_f.

Cas_test_5 construit automatiquement les fibres donc il n'y a aucun réglage manuel à faire.
A_eq et I_eq sont directement injectés dans la matrice de rigidité de l’élément.
------------------------------------------------------------
4) PARAMÈTRES
------------------------------------------------------------

Unités : SI

Longueur : m
Force    : N
Module   : Pa
Section  : m²
Inertie  : m⁴


5) GESTION DES AMPLIFICATIONS (déformée + diagrammes)

Deux modes :

-Mode automatique (ampli_manuel = 0;)

Le script utilise f_ampli.m qui adapte l'amplification proportionellement à la longueur moyenne L_moy de la structure

-Mode manuel (ampli_manuel = 1;)

Les valeurs à fixer (en fractions de L_moy) :

Par défaut :
ampli_def_man = 0.2;   % 20% de L_moy
ampliN_man    = 0.05;  % 5% de L_moy
ampliV_man    = 0.05;  % 5% de L_moy
ampliM_man    = 0.05;  % 5% de L_moy


Le programme convertit automatiquement ces valeurs en amplifications cohérentes.
