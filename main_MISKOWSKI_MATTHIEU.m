addpath("CAS_TEST","FONCTIONS")
Cas_test_1;

for k = 1:n_e
    % Noeuds de début et de fin
    i = B(k,1);
    j = B(k,2);

    % Coordonnées
    x_i = N(i,1);
    y_i = N(i,2);
    x_j = N(j,1);
    y_j = N(j,2);

    % Différences
    x = x_j - x_i;
    y = y_j - y_i;

    % Longueur et inclinaison
    L(k) = sqrt(x^2 + y^2);
    alpha(k) = atan2(y, x); % angle dans [-pi, pi]

end


K = zeros(3*n,3*n)

for e = 1:n_e
    % Noeuds de début et de fin
    i = B(e,1);
    j = B(e,2);
    indice = [(3*i-2):(3*i), (3*j-2):(3*j)];

    ke = f_rigidite(E(e),A(e),I(e),L(e));
    M = f_passage(alpha(e));
    Ke = M'*ke*M;
    K(indice,indice) = K(indice,indice) + Ke

end



%Calcul 
% --- Entrées attendues ---
% K      : matrice de raideur globale (3n x 3n)
% Fimp   : vecteur des forces imposées (3n x 1)
% Uimp   : vecteur des déplacements imposés (3n x 1). Valeur 0 si non imposé.
% DDL    : vecteur indicateur (3n x 1) : 1 = libre, 0 = bloqué

% Exemples: (déjà dans ton cas test)
% load ou définir K, Fimp, Uimp, DDL avant d'exécuter ce bloc

% --- indices libres et bloqués ---
indicesL = find(DDL == 1);   % indices des DDL libres (L)
indicesB = find(DDL == 0);   % indices des DDL bloqués (B)

% --- extraire les blocs de K ---
K_LL = K(indicesL, indicesL);
K_LB = K(indicesL, indicesB);
K_BL = K(indicesB, indicesL);
K_BB = K(indicesB, indicesB);

% --- extraire vecteurs F et U ---
F_L = F_imp(indicesL);
F_B = F_imp(indicesB);

U_B = U_imp(indicesB);      % déplacements imposés (connus) aux DDL bloqués

% --- second membre réduit pour les inconnues U_L ---
F_reduit_L = F_L - K_LB * U_B;


% Résolution  
% On contrôle le conditionnement de KLL, en 
if rcond(K_LL) < 1e-12
    warning('K_LL singulière, des appuis ont été mal définis : ');
    return;
else
    U_L = K_LL \ F_reduit_L;       
end

% --- Reconstruction du vecteur complet des déplacements ---
U = zeros(3*n, 1);   % initialisation

% On place les valeurs connues et calculées
U(indicesL) = U_L;   % déplacements libres calculés
U(indicesB) = U_B;   % déplacements imposés

% --- Question 9 : déplacements des nœuds par élément ---

Ue = zeros(6, n_e);  % déplacements par élément (repère global)
ue = zeros(6, n_e);  % déplacements par élément (repère local)

for e = 1:n_e
    % Noeuds de début et de fin
    i = B(e,1);
    j = B(e,2);

    % Indices des ddl globaux de ces nœuds
    indices = [(3*i-2):(3*i), (3*j-2):(3*j)];

    % Déplacements globaux des noeuds de l'élément e
    Ue(:,e) = U(indices);

    % Matrice de passage (rotation)
    M = f_passage(alpha(e));

    % Déplacements dans le repère local
    ue(:,e) = M * Ue(:,e);
end

% --- Question 10 : Calcul des réactions d'appuis ---

% Calcul des réactions globales
R = K * U - F_imp;       % vecteur de réactions global (3n x 1)

% Extraction des réactions aux DDL bloqués
R_B = R(indicesB);      % seules les réactions sur appuis bloqués

% Création d'un tableau lisible
reac_table = [];
for i = 1:n
    % indices des ddl du noeud i
    ind = (3*i-2):(3*i);
    if any(DDL(ind) == 0)   % s’il y a au moins un ddl bloqué
        reac_table = [reac_table; ...
            i, ...
            R(3*i-2), R(3*i-1), R(3*i)];  % stocker les 3 réactions
    end
end

% --- Question 12 corrigée : Déplacements le long de chaque élément ---

nint = 5;   % nombre d’intervalles
ut = zeros(3, nint+1, n_e);   % déplacements locaux
Ut = zeros(3, nint+1, n_e);   % déplacements globaux

for e = 1:n_e
    Le = L(e);
    M = f_passage(alpha(e));

    % Déplacements nodaux locaux et globaux
    ue_local = ue(:,e);        % 6x1 (local)
    Ue_global = Ue(:,e);       % 6x1 (global)

    % Discrétisation
    x_local = linspace(0, Le, nint+1);

    for k = 1:(nint+1)
        phi = f_forme(x_local(k), Le);

        % déplacement dans le repère local
        ut(:,k,e) = phi(:,:,1) * ue_local;

        % déplacement dans le repère global
        Ut(:,k,e) = phi(:,:,1) * Ue_global;
    end
end


% --- Question 13 : Tracé de la déformée de la structure ---

ampli = 10;  % facteur d'amplification visuelle (ajuster selon les unités)
figure; hold on; axis equal; grid on;
title('Déformée de la structure');
xlabel('X [m]');
ylabel('Y [m]');

for e = 1:n_e
    % --- Coordonnées initiales ---
    i = B(e,1);
    j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);

    % Points initiaux le long de l’élément
    Xini = linspace(xi, xj, nint+1);
    Yini = linspace(yi, yj, nint+1);

    % --- Direction locale de la barre ---
    Lx = (xj - xi);
    Ly = (yj - yi);
    Le = L(e);
    cosA = Lx / Le;
    sinA = Ly / Le;

    % --- Déplacements globaux Ut le long de l’élément ---
    % Ut(1,:,e) : déplacement X local projeté global
    % Ut(2,:,e) : déplacement Y local projeté global
    % On suppose Ut exprimé dans le repère global (déjà fait à Q12)
    Xdef = Xini + ampli * Ut(1,:,e);
    Ydef = Yini + ampli * Ut(2,:,e);

    % --- Tracés ---
    plot(Xini, Yini, 'b--', 'LineWidth', 1);     % forme initiale
    plot(Xdef, Ydef, 'k-', 'LineWidth', 1.8);    % déformée
end

legend('Structure initiale', 'Déformée amplifiée');

% --- Question 14 : Efforts nodaux par élément (repère local) ---

fe = zeros(6, n_e);   % vecteur d'efforts locaux par élément

for e = 1:n_e
    % Récupération des données de l’élément
    Ee = E(e);
    Ae = A(e);
    Ie = I(e);
    Le = L(e);

    % Matrice de rigidité locale
    ke = f_rigidite(Ee, Ae, Ie, Le);

    % Déplacements nodaux locaux
    ue_local = ue(:,e);   % (6x1)

    % Efforts nodaux locaux
    fe(:,e) = ke * ue_local;
end

