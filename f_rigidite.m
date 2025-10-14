function ke = f_rigidite(E, A, I, L)

% Termes de raideur
KN = E*A/L;       % raideur axiale
KM = E*I/(L^3);   % raideur en flexion

ke = [ KN     ,0           ,0      ,-KN      ,0           ,0;
        0   ,12*KM     ,6*L*KM      ,0   ,-12*KM     ,6*L*KM;
        0   ,6*L*KM   ,4*L^2*KM     ,0   ,-6*L*KM   ,2*L^2*KM;
      -KN     ,0           ,0       ,KN      ,0           ,0;
        0  ,-12*KM    ,-6*L*KM      ,0    ,12*KM    ,-6*L*KM;
        0   ,6*L*KM   ,2*L^2*KM     ,0   ,-6*L*KM   ,4*L^2*KM ];

end
