function [EA_eq, EI_eq, yE] = f_section_multifibre(E_f, b_f, h_f, y_f)
% f_section_multifibre
%
%   Entrée :
%       E_f :  vecteur (nf × 1) des modules d'Young des fibres
%       b_f :  vecteur (nf × 1) des largeurs des fibres
%       h_f :  vecteur (nf × 1) des hauteurs des fibres
%       y_f :  vecteur (nf × 1) des positions des centres des fibres
%
%   Sorties :
%       EA_eq : rigidité axiale équivalente
%       EI_eq : rigidité en flexion équivalente (référencée à l’axe neutre équivalent)
%       yE    : position de l’axe neutre équivalent

    % Aire de chaque fibre
    A_f = b_f .* h_f;

    % Rigidité axiale équivalente
    EA_eq = sum(E_f .* A_f);

    % Axe neutre équivalent pondéré par EA
    yE = sum(E_f .* A_f .* y_f) / EA_eq;

    % Inertie propre de chaque fibre
    I_f = b_f .* (h_f.^3) / 12;

    % Rigidité équivalente selon Huygens
    EI_eq = sum(E_f .* (I_f + A_f .* (y_f - yE).^2));

end
