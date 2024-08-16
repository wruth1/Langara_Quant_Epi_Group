function err = IncidenceError5(x, allParams, ImmigrationRate, ReportedIncidence)
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

   
err_Incidence = norm((EstimatedIncidence(1:end-1)'-ReportedIncidence)./ReportedIncidence);
% err = err_Incidence;

    power_relapse = allParams(18);

% % relapse should be account for NgRelapseFraction~=10%
% NgRelapseFraction = 110/1120;
%     p = localBP(2); % ~probability someone in E goes straight into T; pi in Guo-Wu
%     w = localBP(3); % period of time new infectee considered E rather than L
%     v = localBP(4); % rate people in L develop TB
%     u = localBP(8); % relapse rate
% 
%     proportionRelapse =  sum(u*XELTR(:,5))/sum(p*w*XELTR(:,2)+v*XELTR(:,3)+u*XELTR(:,5)) ;
% 
% 
%     % vector
% err_Relapse = ((proportionRelapse-NgRelapseFraction)/NgRelapseFraction).^power_relapse;
% % scalar


    err_Relapse = error_PR(XELTR, localBP, power_relapse);
    
    % err_Relapse0 = error_PR(localIC, localBP, power_relapse);
    % err_Relapse0 = 0;

err = err_Incidence + err_Relapse ;






end
