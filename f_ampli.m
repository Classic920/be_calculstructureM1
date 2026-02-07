function [ampli_def, ampliN, ampliV, ampliM] = f_ampli(Ut, maxN, maxV, maxM, L_moy)
% f_ampli
% Calcule automatiquement :
%   - le facteur d'amplification de la déformée (ampli_def)
%   - les facteurs d'échelle pour les diagrammes N, V, M
%
% Entrées :
%   Ut    : tableau 3 x (nint+1) x n_e des déplacements en repère GLOBAL
%           Ut(1,:,:) = Ux, Ut(2,:,:) = Uy, Ut(3,:,:) = th
%   maxN  : valeur max absolue des efforts normaux (déjà calculée)
%   maxV  : valeur max absolue des efforts tranchants (déjà calculée)
%   maxM  : valeur max absolue des moments fléchissants (déjà calculée)
%   L_moy : longueur moyenne des éléments
%
% Sorties :
%   ampli_def : facteur d'amplification pour la déformée
%   ampliN    : facteur d'échelle pour tracer N(x)
%   ampliV    : facteur d'échelle pour tracer V(x)
%   ampliM    : facteur d'échelle pour tracer M(x)

    % ---------------------------
    % 1) Ampli pour la déformée
    % ---------------------------
    % Norme max des déplacements (Ux,Uy)
    Ux = Ut(1,:,:);
    Uy = Ut(2,:,:);
    U_norm = sqrt(Ux.^2 + Uy.^2);
    maxU = max(U_norm(:));

    if maxU > 0
        % Déformée visible mais pas énorme
        cible = 0.5 * L_moy;   % ~ moitié de la longueur moyenne
        ampli_def = cible / maxU;
    else
        ampli_def = 1.0;
    end

    % --------------------------------------
    % 2) Échelles pour N, V, M
    %    Hauteur des diagrammes plus grande
    % --------------------------------------
    if maxN > 0
        ampliN = 0.1 * L_moy / maxN;   
    else
        ampliN = 0;
    end

    if maxV > 0
        ampliV = 0.15 * L_moy / maxV;   
    else
        ampliV = 0;
    end

    % Moments un peu moins énormes que N,V
    if maxM > 0
        ampliM = 0.2 * L_moy / maxM;
    else
        ampliM = 0;
    end
end
