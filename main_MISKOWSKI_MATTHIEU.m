clearvars Ressort
clearvars Efforts_loc_T % pour pas que matlab garde en mémoire ressort quand j'utilise d'autres cas test
addpath("CAS_TEST","FONCTIONS")

ampli_manuel = 0;  % 0 = auto, 1 = manuel
% Valeurs manuelles d'amplitudes
ampli_def_man = 0.2;   % déformée : .% de L_moy
ampliN_man    = 0.05;  % N :  .% de L_moy
ampliV_man    = 0.05;  % V :  .% de L_moy
ampliM_man    = 0.05;  % M :  .% de L_moy

Cas_test_5;   % charge le jeu de données (noeuds, éléments, E, A, I, DDL, charges...)

% -------------------------------------------------------------------------
% Question 2 
% géométrie des éléments : longueurs et angles (repère global)
% -------------------------------------------------------------------------
for k = 1:n_e
    % nœuds de début et de fin de l’élément k
    i = B(k,1);
    j = B(k,2);

    % coordonnées des nœuds i et j
    x_i = N(i,1);
    y_i = N(i,2);
    x_j = N(j,1);
    y_j = N(j,2);

    % vecteur ij et caractéristiques
    x = x_j - x_i;
    y = y_j - y_i;

    L(k)     = sqrt(x^2 + y^2);   % longueur de l’élément k
    alpha(k) = atan2(y, x);       % orientation de l’élément k (en radians)
end

% -------------------------------------------------------------------------
% Question 5,6 : assemblage de la matrice de raideur globale
% -------------------------------------------------------------------------
K = zeros(3*n,3*n);
for e = 1:n_e
    % nœuds de début et fin
    i = B(e,1);
    j = B(e,2);

    % indices des 6 ddl globaux associés à l’élément e
    indice = [(3*i-2):(3*i), (3*j-2):(3*j)];

    % rigidité locale (repère élément) puis passage en global
    ke = f_rigidite(E(e),A(e),I(e),L(e))
    M  = f_passage(alpha(e));
    Ke = M'*ke*M;

    % ajout dans la matrice globale (assemblage)
    K(indice,indice) = K(indice,indice) + Ke
end

% -- Seulement si on a des appuis élastiques --
if exist('Ressort','var') && ~isempty(Ressort)
    for r = 1:size(Ressort,1)
        noeud    = Ressort(r,1);
        ddl      = Ressort(r,2);     % 1=Ux, 2=Uy, 3=rotation
        kressort = Ressort(r,3);     % N/m ou N·m/rad pour dof=3
        g = 3*(noeud-1) + ddl;       % indice global du ddl
        K(g,g) = K(g,g) + kressort;
        % un ressort = ddl libre
        DDL(g)   = 1;
        U_imp(g) = 0;
    end
end

% -- Seulement si on a des efforts centrés sur barre --
% -------------------------------------------------------------------------
% Ajout des efforts équivalents d'éléments (effort ponctuel centré)
% Conventions :
%   - Efforts_loc_T est un tableau (nb_lignes x 2) : [ e, P ]
%     * e = indice d'élément (1..n_e)
%     * P = intensité (P > 0 agit vers -y_local)
% -------------------------------------------------------------------------
F_eq = zeros(3*n,1);  % second membre équivalent global

if exist('Efforts_loc_T','var') && ~isempty(Efforts_loc_T)
    for e = 1:n_e
        % filtre : toutes les charges déclarées sur l'élément e
        mask = (Efforts_loc_T(:,1) == e);
        if any(mask)
            Le = L(e);
            a  = alpha(e);
            M  = f_passage(a);  % 6x6 (local -> global)
            feq_loc_sum = zeros(6,1);

            % somme des contributions si plusieurs charges sur le même élément
            for r = find(mask).'
                P = Efforts_loc_T(r,2);
                feq_loc_sum = feq_loc_sum + f_charge_equivalente_point_centre(Le, P);
            end

            % passage en global pour l'assemblage
            feq_glob = M.' * feq_loc_sum;

            % indices des DDL globaux de l'élément
            i = B(e,1); j = B(e,2);
            indice = [(3*i-2):(3*i), (3*j-2):(3*j)];

            % assemblage dans F_eq
            F_eq(indice) = F_eq(indice) + feq_glob;
        end
    end
end

% Second membre total : forces nodales + équivalents d’éléments
F_tot = F_imp + F_eq;


% =========================================================================
% Question 7  : on sépare K en quatre blocs K_LL, K_LB, K_BL, K_BB
% =========================================================================

% répartition des ddl libres/bloqués
indicesL = find(DDL == 1);   % ddl libres
indicesB = find(DDL == 0);   % ddl bloqués

% sous-blocs de K selon L/B
K_LL = K(indicesL, indicesL);
K_LB = K(indicesL, indicesB);
K_BL = K(indicesB, indicesL);
K_BB = K(indicesB, indicesB);

% vecteurs de charges et de déplacements imposés filtrés
F_L = F_tot(indicesL);
F_B = F_tot(indicesB);
U_B = U_imp(indicesB);       % déplacements imposés sur les ddl bloqués

% second membre réduit (équilibre sur ddl libres)
F_reduit_L = F_L - K_LB * U_B;

U_L = K_LL \ F_reduit_L;

% -------------------------------------------------------------------------
% Question 8 : recomposition du vecteur de déplacements complet
% -------------------------------------------------------------------------
U = zeros(3*n, 1);
U(indicesL) = U_L;   % inconnues résolues
U(indicesB) = U_B;   % valeurs imposées


% -------------------------------------------------------------------------
% Question 9 : déplacements nodaux par élément
% -------------------------------------------------------------------------
Ue = zeros(6, n_e);  % global
ue = zeros(6, n_e);  % local

for e = 1:n_e
    % nœuds de l’élément
    i = B(e,1);
    j = B(e,2);

    % indices des 6 ddl globaux de l’élément
    indices = [(3*i-2):(3*i), (3*j-2):(3*j)];

    % déplacements aux nœuds (global)
    Ue(:,e) = U(indices);

    % passage en local via la matrice de rotation
    M = f_passage(alpha(e));
    ue(:,e) = M * Ue(:,e);
end


% -------------------------------------------------------------------------
% Question 10 : réactions d’appui
% -------------------------------------------------------------------------
R = K * U - F_tot;   % vecteur global des réactions (3n x 1)
R_B = R(indicesB);   % réactions sur ddl bloqués (si besoin)

% tableau 
reac_table = [];
for i = 1:n
    ind = (3*i-2):(3*i);
    if any(DDL(ind) == 0)
        reac_table = [reac_table; ...
            i, R(3*i-2), R(3*i-1), R(3*i)];
    end
end

% -------------------------------------------------------------------------
% Question 12  : déplacements le long de chaque élément
% -------------------------------------------------------------------------
nint = 100;                         % nombre de segments pour l’affichage
ut = zeros(3, nint+1, n_e);         % local : [u; v; th]
Ut = zeros(3, nint+1, n_e);         % global : [Ux; Uy; th]

for e = 1:n_e
    Le = L(e);
    a  = alpha(e);

    % déplacements nodaux LOCAUX [u1; v1; th1; u2; v2; th2]
    ue_local  = ue(:,e);
    x_local   = linspace(0, Le, nint+1);

    % matrice de passage GLOBAL -> LOCAL (6x6)
    M6 = f_passage(a);
    % sous-bloc rotations/déplacements (global -> local) pour (Ux,Uy,th)
    T3_glob2loc = M6(1:3,1:3);
    % pour passer de local -> global : on prend la transposée
    T3_loc2glob = T3_glob2loc.';  

    for k = 1:(nint+1)
        phi = f_forme(x_local(k), Le);   % taille 3x6x1
        ut(:,k,e) = phi(:,:,1) * ue_local;
        Ut(:,k,e) = T3_loc2glob * ut(:,k,e);
    end 
end


% -------------------------------------------------------------------------
% Question 14 : efforts nodaux par élément en repère local
% (ke * ue_local - charges équivalentes d’élément éventuelles)
% -------------------------------------------------------------------------
fe = zeros(6, n_e);   % efforts nodaux LOCAUX [N1; V1; M1; N2; V2; M2]

for e = 1:n_e
    Ee = E(e); Ae = A(e); Ie = I(e); Le = L(e);

    % matrice de rigidité LOCALE (ta fonction)
    ke = f_rigidite(Ee, Ae, Ie, Le);

    % charges équivalentes locales (si chargement ponctuel au centre)
    fe_eq_loc = zeros(6,1);
    if exist('Efforts_loc_T','var') && ~isempty(Efforts_loc_T)
        mask = (Efforts_loc_T(:,1) == e);
        if any(mask)
            Ptot = 0;
            for r = find(mask).'
                Ptot = Ptot + Efforts_loc_T(r,2);
            end
            fe_eq_loc = f_charge_equivalente_point_centre(Le, Ptot);
        end
    end

    % efforts locaux internes résultants
    fe(:,e) = ke * ue(:,e) - fe_eq_loc;
end


% -------------------------------------------------------------------------
% Question 16 : diagrammes N(x), V(x), M(x) + déformée
% -------------------------------------------------------------------------
npoints = 200;
s = linspace(0,1,npoints);   % abscisse réduite 0→1
Xlim = [min(N(:,1)) max(N(:,1))];
Ylim = [min(N(:,2)) max(N(:,2))];
L_moy = mean(L);

% ===========================
% 1) Pré-calcul N, V, M + max
% ===========================

maxN = eps; clear Nplot
for e = 1:n_e
    Le = L(e); 
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);

    % efforts nodaux LOCAUX : N1, N2
    N1 = fe(1,e);
    N2 = -fe(4,e);       % signe opposé à droite
    Nbar = 0.5*(N1+N2);  % poutre de Bernoulli → N ~ constant

    % inversion pour coller au schéma de l’énoncé
    Nloc = -Nbar * ones(1,npoints);
    maxN = max(maxN, max(abs(Nloc)));

    % géométrie de l’axe
    X = xi + (xj - xi).*s;
    Y = yi + (yj - yi).*s;

    % vecteur normal géométrique (global)
    t = [xj - xi, yj - yi]; t = t / norm(t);
    nvec = [-t(2), t(1)];

    Nplot(e).X = X; Nplot(e).Y = Y; Nplot(e).N = Nloc; Nplot(e).n = nvec;
end

maxV = eps; clear Vplot
for e = 1:n_e
    Le = L(e);
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);

    X = xi + (xj - xi).*s;
    Y = yi + (yj - yi).*s;

    % efforts nodaux Locaux
    M1 = fe(3,e);       % moment à gauche
    M2 = -fe(6,e);      % moment à droite

    Vbar = (M2 - M1)/Le;

    % charge ponctuelle au milieu ?
    Ptot = 0;
    if exist('Efforts_loc_T','var') && ~isempty(Efforts_loc_T)
        mask = (Efforts_loc_T(:,1) == e);
        if any(mask), Ptot = sum(Efforts_loc_T(mask,2)); end
    end

    if abs(Ptot) > 0
        Vleft  = Vbar + Ptot/2;
        Vright = Vbar - Ptot/2;
        Vloc = Vleft*ones(1,npoints);
        Vloc(round(npoints/2):end) = Vright; % saut au milieu
    else
        Vloc = Vbar*ones(1,npoints);
    end

    maxV = max(maxV, max(abs(Vloc)));

    % normale géométrique pour tracé
    t = [xj - xi, yj - yi]; t = t/norm(t);
    nvec = [-t(2), t(1)];

    Vplot(e).X = X; Vplot(e).Y = Y; Vplot(e).V = Vloc; Vplot(e).n = nvec;
end

maxM = eps; clear Mplot
for e = 1:n_e
    Le = L(e);
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);

    X = xi + (xj - xi).*s;
    Y = yi + (yj - yi).*s;

    % normale géométrique pour décalage
    t = [xj - xi, yj - yi]; t = t/norm(t);
    nvec = [-t(2), t(1)];

    % efforts nodaux LOCAUX
    M1 = fe(3,e);
    M2 = -fe(6,e);
    Vbar = (M2 - M1)/Le;

    % charge ponctuelle au milieu ?
    Ptot = 0;
    if exist('Efforts_loc_T','var') && ~isempty(Efforts_loc_T)
        mask = (Efforts_loc_T(:,1) == e);
        if any(mask), Ptot = sum(Efforts_loc_T(mask,2)); end
    end

    xloc = s*Le;
    Mloc = zeros(1,npoints);

    if abs(Ptot) > 0
        Vleft  = Vbar + Ptot/2;
        Vright = Vbar - Ptot/2;

        for k = 1:npoints
            xk = xloc(k);
            if xk <= Le/2
                Mloc(k) = M1 + Vleft * xk;
            else
                Mmid = M1 + Vleft*(Le/2);
                Mloc(k) = Mmid + Vright*(xk - Le/2);
            end
        end
    else
        % pas de charge interne : M linéaire continu
        Mloc = M1 + Vbar * xloc;
    end

    maxM = max(maxM, max(abs(Mloc)));

    Mplot(e).X = X; Mplot(e).Y = Y; Mplot(e).M = Mloc; Mplot(e).n = nvec;
end

% -------------------------------------------------------------------------
% Choix des amplitudes (automatiques ou manuelles)
% -------------------------------------------------------------------------
L_moy = mean(L);   % si pas déjà calculé plus haut

% norme max des déplacements pour la déformée
Ux = Ut(1,:,:);
Uy = Ut(2,:,:);
U_norm = sqrt(Ux.^2 + Uy.^2);
maxU = max(U_norm(:));

if ampli_manuel == 1
    % Mode manuel pour ampli

    if maxU > 0
        ampli_def = (ampli_def_man * L_moy) / maxU;
    else
        ampli_def = 1.0;
    end

    if maxN > 0
        ampliN = (ampliN_man * L_moy) / maxN;
    else
        ampliN = 0;
    end

    if maxV > 0
        ampliV = (ampliV_man * L_moy) / maxV;
    else
        ampliV = 0;
    end

    if maxM > 0
        ampliM = (ampliM_man * L_moy) / maxM;
    else
        ampliM = 0;
    end

else
    % ----- Géré automatiqement par la fonction ----
    [ampli_def, ampliN, ampliV, ampliM] = f_ampli(Ut, maxN, maxV, maxM, L_moy);
end


% ===========================
% 2) Tracé avec ces amplitudes
% ===========================
figure; tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

%% (1) Déformée
nexttile; hold on; axis equal; grid on;
title('Déformée'); xlabel('X [m]'); ylabel('Y [m]');
ampli = ampli_def;

for e = 1:n_e
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);

    Xini = linspace(xi, xj, nint+1);
    Yini = linspace(yi, yj, nint+1);
    Xdef = Xini + ampli * Ut(1,:,e);
    Ydef = Yini + ampli * Ut(2,:,e);

    plot(Xini, Yini, 'b:', 'LineWidth', 1);
    plot(Xdef, Ydef, 'k-', 'LineWidth', 1.8);
end
xlim([Xlim(1)-0.1*L_moy, Xlim(2)+0.1*L_moy]);
ylim([Ylim(1)-0.1*L_moy, Ylim(2)+0.1*L_moy]);

%% (2) Efforts normaux N(x)
nexttile; hold on; axis equal; grid on;
title('Efforts normaux'); xlabel('X [m]'); ylabel('Y [m]');
for e = 1:n_e
    % axe géométrique
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);
    plot([xi,xj],[yi,yj],'b:','LineWidth',1);

    % diagramme N
    plot(Nplot(e).X + ampliN*Nplot(e).N.*Nplot(e).n(1), ...
         Nplot(e).Y + ampliN*Nplot(e).N.*Nplot(e).n(2), ...
         'b','LineWidth',1.5);
end
xlim([Xlim(1)-0.1*L_moy, Xlim(2)+0.1*L_moy]);
ylim([Ylim(1)-0.1*L_moy, Ylim(2)+0.1*L_moy]);

%% (3) Efforts tranchants V(x)
nexttile; hold on; grid on;
title('Efforts tranchants'); xlabel('X [m]'); ylabel('Y [m]');
for e = 1:n_e
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);
    plot([xi,xj],[yi,yj],'b:','LineWidth',1);  % axe

    plot(Vplot(e).X + ampliV*Vplot(e).V.*Vplot(e).n(1), ...
         Vplot(e).Y + ampliV*Vplot(e).V.*Vplot(e).n(2), ...
         'r','LineWidth',1.5);
end
xlim([Xlim(1)-0.1*L_moy, Xlim(2)+0.1*L_moy]);
ylim([Ylim(1)-0.1*L_moy, Ylim(2)+0.1*L_moy]);

%% (4) Moments fléchissants M(x)
nexttile; hold on; grid on;
title('Moments fléchissants'); xlabel('X [m]'); ylabel('Y [m]');
for e = 1:n_e
    i = B(e,1); j = B(e,2);
    xi = N(i,1); yi = N(i,2);
    xj = N(j,1); yj = N(j,2);
    plot([xi,xj],[yi,yj],'b:','LineWidth',1); % axe géométrique

    plot(Mplot(e).X - ampliM*Mplot(e).M.*Mplot(e).n(1), ...
          Mplot(e).Y - ampliM*Mplot(e).M.*Mplot(e).n(2), ...
          'g','LineWidth',1.5);
end
xlim([Xlim(1)-0.1*L_moy, Xlim(2)+0.1*L_moy]);
ylim([Ylim(1)-0.1*L_moy, Ylim(2)+0.1*L_moy]);

%% ================================
%  Note de calcul : tableaux de synthèse
% ================================

% -------------------------------
% 1) Table des déplacements aux nœuds
% -------------------------------
% U est de taille 3n x 1 : [Ux1; Uy1; Th1; Ux2; Uy2; Th2; ...]
U_node = reshape(U, 3, n)';   % n x 3 : [Ux Uy Th]

% norme du déplacement au nœud (sans la rotation)
Unorm = sqrt(U_node(:,1).^2 + U_node(:,2).^2);

TableDeplNoeuds = table((1:n)', ...
                        U_node(:,1), U_node(:,2), U_node(:,3), Unorm, ...
    'VariableNames', {'Noeud','Ux_m','Uy_m','Theta_rad','Unorm_m'})

% -------------------------------
% 2) Table des déplacements max
% -------------------------------
[Umax, indMax] = max( Unorm );

TableDeplMax = table(indMax, Umax, ...
    'VariableNames', {'Noeud_max','Deplacement_max_m'})

% -------------------------------
% 3) Efforts de cohésion min et max (N, V, M)
% -------------------------------
N_tout = [];
V_tout = [];
M_tout = [];

for e = 1:n_e
    N_tout = [N_tout, Nplot(e).N];
    V_tout = [V_tout, Vplot(e).V];
    M_tout = [M_tout, Mplot(e).M];
end

N_min = min(N_all);  N_max = max(N_all);
V_min = min(V_all);  V_max = max(V_all);
M_min = min(M_all);  M_max = max(M_all);

TableEffortsExtremes = table(N_min, N_max, V_min, V_max, M_min, M_max, ...
    'VariableNames', {'N_min','N_max','V_min','V_max','M_min','M_max'})

% -------------------------------
% 4) Table des réactions d’appui
% -------------------------------
%   [Noeud, Rx, Ry, Mz] pour chaque nœud avec au moins un ddl bloqué

TableReactions = array2table(reac_table, ...
    'VariableNames', {'Noeud','Rx_N','Ry_N','Mz_Nm'})
