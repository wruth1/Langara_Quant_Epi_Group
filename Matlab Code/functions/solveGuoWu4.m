function [XELTR, TBIncidence, TBPrevalence] = solveGuoWu4(bioParameters, initialConditions, ImmigrationRate)
%{
Updated to include sigma, the reactivation rate of R to T
  Input: 
    - bioParameters vector in R9, hard coded so position [beta p w v a d n q1 q2]
    - initialConditions vector in R5, hard coded so position [X0 E0 L0 T0
    R0]. 
    - ImmigrationRate, which gives annual immigration rate, vector in R^N

  Output:
    - XELTR, an N+1 (N is size of ImmigrationRate) by 5 matrix where each row is a time
    - TBIncidence, vector in R^N+1 with calculated incidence
%}

%% Unpack parameters

 beta = bioParameters(1); %TB infectivity
 p = bioParameters(2); % ~probability someone in E goes straight into E; pi in Guo-Wu
 w = bioParameters(3); % period of time new infectee considered E rather than L
 v = bioParameters(4); % rate people in L develop TB
 a = bioParameters(5); % TB death rate for people in T
 d = bioParameters(6); % non-TB death rate
 n = bioParameters(7); % recovery rate; delta in Guo-Wu
 u = bioParameters(8); % percentage people in R develop active TB 
 q1 = bioParameters(9); % percent of immigrants into E
 q2 = bioParameters(10); % percent of immigrants into L
 q3 = bioParameters(11); % percent of immigrants into R

 
 %% Initialize outputs
 Nt = length(ImmigrationRate);
 
 XELTR = zeros(Nt+1, 5);
  X0 = initialConditions(1);
  E0 = initialConditions(2);
  L0 = initialConditions(3);
  T0 = initialConditions(4);
  R0 = initialConditions(5);
  
 XELTR0 = [X0 E0 L0 T0 R0]; % this will be replaced soon
 XELTR(1,:) = XELTR0;
 
 
 TBIncidence = zeros(Nt+1,1);
 TBPrevalence = zeros(Nt+1,1);
 
 TBIncidence(1) = getTBIncidenceRate3(XELTR0, p, w, v, u);
 TBPrevalence(1) = XELTR0(4);

%% Calculate TB Incididence 

    XELTRlocal = XELTR0;
    for i = 1:Nt

        pi = ImmigrationRate(i);
        
        % y = (x,e,l,t,r) is solution to system of odes.
        [t,y] = ode23(@(t,y) odefcn(t, y), [0,1], XELTRlocal); % [0,1] means 1 year forward

        % update outputs
        XELTRlocal = y(end,:);
        XELTR(i+1,:) = XELTRlocal;

        TBIncidence(i+1) = getTBIncidenceRate3(XELTRlocal, p, w, v, u);

        TBPrevalence(i+1) = XELTRlocal(4);
        % getTBPrevalence2(XELTR0local);




    end
%% 
    function Tri = getTBIncidenceRate3(y, p, w, v, u)
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




%% %% Setup the ODE System
    function dydt = odefcn(t, y  )

    %pi = ImmigrationRate(t);
 
dydt = zeros(5,1);
dydt(1) = (1-q1-q2-q3)*pi - beta*y(1)*y(4) - n*y(1);
dydt(2) = q1*pi + beta*y(1)*y(4)-(n+w)*y(2);
dydt(3) = q2*pi + (1-p)*w*y(2) - (n+v)*y(3);
dydt(4) = p*w*y(2) + v*y(3) - (n+a+d)*y(4) + u*y(5);
dydt(5) = q3*pi + d*y(4) - n*y(5) - u*y(5);
end

end 