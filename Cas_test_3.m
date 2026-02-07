%%%%%%%%%%%%%%%%%%%%%
%%Données d'entrées%%
%%%%%%%%%%%%%%%%%%%%%

n = 4; % Nombre de noeuds 
n_e = 3; % Nombre d'éléments (poutres)
N = [[0,0];[10,10];[25,0];[10,0]]; % Coordonnées des noeuds (taille n   × 2)
B = [[1,2];[4,2];[2,3]]; % Numéros des nœuds de début et de fin par élément ((taille n_e × 2)
E = [210e9;210e9;210e9]; % Module d'Young par élément (taille n_e × 1)
A = [10.3e-4;40e-4;10.3e-4]; % Aire par élément (taille n_e × 1)
I = [171e-8;850e-8;171e-8]; % Moment quadratique par élément (taille n_e × 1)
DDL = [1;0;1;
    1;1;1;
    1;1;1;
    0;0;0];   % Degrés de liberté DDL (taille 3n × 1)
F_imp = [0;0;0;
    3e3;0;0;
    0;0;0;
    0;0;0];  % Forces et couples imposées aux nœuds (taille 3n × 1)
U_imp = [0;0;0;
    0;0;0;
    0;0;0;
    0;0;0]; % Déplacements et rotations imposés aux nœuds (taille 3n×1)

Ressort = [3,2,5e2] %[noeud,n°DDL(1=Ux,2=Uy,3=rotation), k]
%Donc ici j'ai mis un ressort au noeud 3 selon Uy de raideur 500 N/m