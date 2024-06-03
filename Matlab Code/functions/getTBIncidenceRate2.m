

%% calculate the TB Incidence Rate
function Tri = getTBIncidenceRate2(y, p, w, v, u)
% also grabs  from R
    Xi = y(:,1);
    Ei = y(:,2);
    Li = y(:,3);
    Ti = y(:,4);
    Ri = y(:,5);

    TFPi = Xi + Ei + Li + Ti + Ri; %total Foreign Population
    
    Tri = p*w*Ei + v*Li + u*Ri; %Page 10(704) of the York paper
    Tri = Tri * 100000./TFPi;

end 

