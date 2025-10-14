%%%%%%%%%%%%%%%%%%%%%
%%Données d'entrées%%
%%%%%%%%%%%%%%%%%%%%%

n = % Nombre de noeuds 
n_e = % Nombre d'éléments (poutres)
N = % Coordonnées des noeuds (taille n   × 2)
B = % Numéros des nœuds de début et de fin par élément ((taille n_e × 2)
E = % Module d'Young par élément (taille n_e × 1)
A = % Aire par élément (taille n_e × 1)
I = % Moment quadratique par élément (taille n_e × 1)
DDL = % Degrés de liberté DDL (taille 3n × 1)
F_imp = % Forces et couples imposées aux nœuds (taille 3n × 1)
U_imp = % Déplacements et rotations imposés aux nœuds (taille 3n×1)