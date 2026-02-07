function M = f_passage(alpha)
% Question 3 
% Calcule la matrice de passage (rotation) 6x6 entre le repère local et global.
% Entrée :
%   alpha : angle d’inclinaison de l’élément [rad]
% Sortie :
%   M : matrice de rotation 6x6 
%
% elle permet de transformer les déplacements ou efforts :
%   - du repère local de la poutre (axe de la barre)
%   - vers le repère global de la structure.

M = [ cos(alpha),  sin(alpha),  0,   0,  0,  0;
     -sin(alpha),  cos(alpha),  0,   0,  0,  0;
      0,           0,           1,   0,  0,  0;
      0,           0,           0,   cos(alpha),  sin(alpha),  0;
      0,           0,           0,  -sin(alpha),  cos(alpha),  0;
      0,           0,           0,   0,           0,           1];
end
