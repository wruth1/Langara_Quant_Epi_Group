function err = error_Prevalence(XELTR,  ReportedPrevalence)

    % numt = length(ReportedPrevalence);

    % compute "actual" = reported prevalence ratio p
    TotalPrevalence_report = sum(ReportedPrevalence);
    p = ReportedPrevalence/TotalPrevalence_report; % function of time
    
    % now compute our estimated p
    That = XELTR(1:end-1,4)'; %estimated prevalence
    TotalPrevalence_hat = sum(That);
    phat = That/TotalPrevalence_hat;


    %   weight all terms equally
    err = mean( ( (phat-p)./p).^2);

    % special weight on first term
    % err = mean( ( (phat(2:end)-p(2:end))./p(2:end)).^2);
    % err = err + ((phat(1)-p(1))/p(1))^2;
    

    
   
   

end