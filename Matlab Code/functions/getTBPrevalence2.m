

%% calculate the TB Incidence Rate
function TB_Prevalence = getTBPrevalence2(y)
% also grabs  from R
    Xi = y(:,1);
    Ei = y(:,2);
    Li = y(:,3);
    Ti = y(:,4);
    Ri = y(:,5);

    TFPi = Xi + Ei + Li + Ti + Ri; %total Foreign Population
    
    
    TB_Prevalence = Ti * 100000./TFPi;

end 

