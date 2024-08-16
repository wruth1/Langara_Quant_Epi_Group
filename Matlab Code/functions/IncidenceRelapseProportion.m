function proportionRelapse = IncidenceRelapseProportion(XELTR, BP)


    % NgRelapseFraction = 110/1120;
    p = BP(2); % ~probability someone in E goes straight into T; pi in Guo-Wu
    w = BP(3); % period of time new infectee considered E rather than L
    v = BP(4); % rate people in L develop TB
    u = BP(8); % relapse rate

    proportionRelapse =  sum(u*XELTR(:,5))/sum(p*w*XELTR(:,2)+v*XELTR(:,3)+u*XELTR(:,5)) ;




end