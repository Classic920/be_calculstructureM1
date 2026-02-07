function ke = f_rigidite(E, A, I, L)
% Question 4
% Calcule la matrice de rigidité locale (6x6) d’un élément de poutre.
% Entrées :
%   E : module d’Young [Pa = N/m²]
%   A : aire de la section [m²]
%   I : moment quadratique [m⁴]
%   L : longueur de l’élément [m]
% Sortie :
%   ke : matrice de rigidité locale 
%
% cette matrice relie les déplacements et rotations locaux (u, v, theta)
% aux efforts internes correspondants (N, V, M) de l’élément.

KN = E*A/L;       
KM = E*I/L^3;     

ke = [ KN,     0,        0,    -KN,     0,        0;
        0,  12*KM,   6*L*KM,     0, -12*KM,   6*L*KM;
        0,   6*L*KM, 4*L^2*KM,   0,  -6*L*KM, 2*L^2*KM;
      -KN,     0,        0,     KN,     0,        0;
        0, -12*KM,  -6*L*KM,     0,  12*KM,  -6*L*KM;
        0,   6*L*KM, 2*L^2*KM,   0,  -6*L*KM, 4*L^2*KM];
end
