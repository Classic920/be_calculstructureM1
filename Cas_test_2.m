%%%%%%%%%%%%%%%%%%%%%
%%Données d'entrées%%
%%%%%%%%%%%%%%%%%%%%%

n = 9; % Nombre de noeuds 
n_e = 14; % Nombre d'éléments (poutres)
N = [[0,0];[6,0];[12,0];[0,4];[6,4];[12,4];[0,8];[6,8];[12,8]]; % Coordonnées des noeuds (taille n × 2)
B = [[1,4];[4,7];[2,5];[5,8];[3,6];[6,9];[4,5];[5,6];[7,8];[8,9];[4,8];[5,7];[5,9];[6,8]]; % Numéros des nœuds de début et de fin par élément (taille n_e × 2)
E = [210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9;210e9]; % Module d'Young par élément (taille n_e × 1)
A = [0.9e-3;0.9e-3;0.9e-3;0.9e-3;0.9e-3;0.9e-3;1.2e-3;1.2e-3;1.2e-3;1.2e-3;0.6e-3;0.6e-3;0.6e-3;0.6e-3]; % Aire par élément (taille n_e × 1)
I = [180e-8;180e-8;180e-8;180e-8;180e-8;180e-8;240e-8;240e-8;240e-8;240e-8;60e-8;60e-8;60e-8;60e-8]; % Moment quadratique par élément (taille n_e × 1)
DDL = [0;0;0;
       0;0;0;
       0;0;0;
       1;1;1;
       1;1;1;
       1;1;1;
       1;1;1;
       1;1;1;
       1;1;1];   % Degrés de liberté DDL (taille 3n × 1)
F_imp = [0;0;0;
         0;0;0;
         0;0;0;
         0;-2.0e4;0;
         5.0e4;-1.0e5;0;
         0;-1.5e5;0;
         0;-8.0e4;0;
         5.0e4;0;0;
         0;-1.2e5;0];  % Forces et couples imposées aux nœuds (taille 3n × 1)
U_imp = [0;0;0;
         0;0;0;
         0;0;0;
         0;0;0;
         0;0;0;
         0;0;0;
         0;0;0;
         0;0;0;
         0;0;0]; % Déplacements et rotations imposés aux nœuds (taille 3n×1)