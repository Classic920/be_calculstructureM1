%%%%%%%%%%%%%%%%%%%%%
%% Données d'entrées
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
    1;0;1;
    0;0;0];   % Degrés de liberté DDL (taille 3n × 1)
F_imp = [0;0;0;
    3e3;0;0;
    0;0;0;
    0;0;0];  % Forces et couples imposées aux nœuds (taille 3n × 1)
U_imp = [0;0;0;
    0;0;0;
    0;0;0;
    0;0;0]; % Déplacements et rotations imposés aux nœuds (taille 3n×1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MULTIFIBRE : paramètres pour l’élément 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_MF = 2;          % indice de l'élément traité en multifibre

% Paramètres géométriques et matériaux pour la section multifibre
h_MF    = 0.16;     % hauteur totale (m)
A_MF    = 40e-4;    % aire totale (m^2) 
E_basMF = 210e9;    % module de Young de la zone basse (Pa)
E_hautMF= 420e9;    % module de Young de la zone haute (Pa)
Ny_MF   = 40;       % nombre de couches en hauteur 

% Largeur totale équivalente
b_MF = A_MF / h_MF;

% Discrétisation en Ny_MF fibres sur la hauteur
y_bords = linspace(-h_MF/2, h_MF/2, Ny_MF+1)';   
h_f     = diff(y_bords);                         
b_f     = b_MF * ones(Ny_MF,1);                  
y_f     = (y_bords(1:end-1) + y_bords(2:end))/2; 

% Modules de Young par fibre :
E_f = [E_basMF  * ones(Ny_MF/2,1);
       E_hautMF * ones(Ny_MF/2,1)];

% Calcul EA_eq, EI_eq, yE 
[EA_eq, EI_eq, yE] = f_section_multifibre(E_f, b_f, h_f, y_f);

% Mise à jour de E, A, I pour l'élément multifibre
E_ref       = E_basMF;       
E(e_MF)     = E_ref;
A(e_MF)     = EA_eq / E_ref;
I(e_MF)     = EI_eq / E_ref;
