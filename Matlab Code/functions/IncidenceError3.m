function err = IncidenceError3(x, bioParameters, initialConditions, ImmigrationRate, ReportedIncidence,ReportedPrevalence)
% input: x is in R2, q1, q2

localBP = bioParameters;
localBP(9)=x(1); % q1
localBP(10)=x(2); % q2
localBP(11) = x(3); %q3
% localBP(8) = x(4); % u


%compute initialConditions
ysteady = findSteadyState2(localBP, initialConditions, ImmigrationRate(1));
% rescale to correct total population.  localIC(1) has correct TP0
ysteady = ysteady*sum(initialConditions)/sum(ysteady);
localIC = ysteady;
    


[XELTR, EstimatedIncidence, EstimatedPrevalence] = solveGuoWu4(localBP, localIC, ImmigrationRate);

   
err_Incidence = norm((EstimatedIncidence(1:end-1)'-ReportedIncidence)./ReportedIncidence);
err = err_Incidence;

% relapsde error
% relapse should be account for NgRelapseFraction~=10%
% NgRelapseFraction = 110/1120;
%     p = localBP(2); % ~probability someone in E goes straight into T; pi in Guo-Wu
%     w = localBP(3); % period of time new infectee considered E rather than L
%     v = localBP(4); % rate people in L develop TB
%     % u = localBP(8); % relapse rate
% 
%     proportionRelapse =  sum(u*XELTR(:,5))/sum(p*w*XELTR(:,2)+v*XELTR(:,3)+u*XELTR(:,5)) ;
% 
% 
%     power_relapse = 100;
%     % vector
% err_Relapse = ((proportionRelapse-NgRelapseFraction)/NgRelapseFraction).^power_relapse;
% % scalar
% err = err_Incidence + err_Relapse;






end
