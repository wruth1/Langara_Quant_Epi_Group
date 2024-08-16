function err = IncidenceError6(x, allParams, ImmigrationRate, ReportedIncidence)
% input: x is in R2, q1, q2

localBP = allParams(1:11);
localBP(9)=x(1); % q1
localBP(10)=x(2); % q2
localBP(11) = x(3); %q3
localBP(8) = x(4); % u

localIC(1) = x(5); %X0
localIC(2) = x(6); %E0
localIC(3) = x(7); %L0
localIC(4) = x(8); %T0
localIC(5) = x(9); %R0


[XELTR, EstimatedIncidence] = solveGuoWu4(localBP, localIC, ImmigrationRate);

   


err_Incidence = sum(((EstimatedIncidence(1:end-1)'-ReportedIncidence)./ReportedIncidence).^2);
% err = err_Incidence;

    power_relapse = allParams(18);

    err_Relapse = error_PR(XELTR, localBP, power_relapse);
    
    err_Relapse0 = error_PR(localIC, localBP, power_relapse);
    % err_Relapse0 = 0;


    % JordanEoverTB_mean = 1/488;
    % err_E2021ratio = error_E2021(XELTR, JordanEoverTB_mean);
    err_E2021ratio = 0;

    JordanX0Ratio = 0.76;
    err_JordanX0 = error_JordanX0(XELTR,  JordanX0Ratio);
    
err = err_Incidence + err_Relapse + err_Relapse0 + err_E2021ratio + err_JordanX0 ;






end
