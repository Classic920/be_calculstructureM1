function phi = f_forme(x, L)
% Question 11
% Calcule la matrice des fonctions de forme d’un élément de poutre.
% Entrées :
%   x : position(s) locale(s) le long de l’élément [m] (scalaire ou vecteur, 0 <= x <= L)
%   L : longueur de l’élément [m]
% Sortie :
%   phi : tableau [3 x 6 x nx]
%         → pour chaque position x(k), phi(:,:,k) donne la matrice reliant
%           les déplacements nodaux locaux [u_i, v_i, theta_i, u_j, v_j, theta_j]
%           aux déplacements internes [u, v, theta] au point x(k)
%
%

% mise en forme des données d’entrée
x = x(:)';             
nx = numel(x);
s  = x / L;            % position adimensionnée
s2 = s.^2; s3 = s.^3;

% on calcule les fonctions de forme
phi1 = 1 - s;
phi4 = s;
phi8  = 1 - 3*s2 + 2*s3;
phi9  = x - 2*(x.^2)/L + (x.^3)/(L^2);      
phi11 = 3*s2 - 2*s3;
phi12 = - (x.^2)/L + (x.^3)/(L^2);         

phi8p  = (-6*s + 6*s2) / L;                
phi9p  = 1 - 4*x/L + 3*(x.^2)/(L^2);       
phi11p = (6*s - 6*s2) / L;                  
phi12p = -2*x/L + 3*(x.^2)/(L^2);           

%on introduit un tableau de zéros [3x6xnx]
phi = zeros(3,6,nx);

% Assemblage de la matrice des fonctions de forme pour chaque x(k)
for k = 1:nx
    phi(:,:,k) = [ ...
        phi1(k),     0,        0,      phi4(k),    0,         0;      % déplacement axial
            0,    phi8(k),  phi9(k),     0,     phi11(k),  phi12(k);  % déplacement transversal
            0,    phi8p(k), phi9p(k),    0,     phi11p(k), phi12p(k)];% rotation
end
end
