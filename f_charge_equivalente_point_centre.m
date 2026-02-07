function feq_loc = f_charge_equivalente_point_centre(L, P)
% Effort ponctuel transversal P au milieu, dirigé vers -y_local
% Retourne le vecteur 6x1 des efforts nodaux équivalents (repère local)
% Ordre des ddl locaux : [u1; v1; th1; u2; v2; th2]
feq_loc = [ 0;
           -P/2;
           -P*L/8;
            0;
           -P/2;
            P*L/8 ];
end
