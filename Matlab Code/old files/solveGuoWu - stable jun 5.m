function [XELTR, TBestimate] = solveGuoWu(params1, params2, params3, ImmigrationRate)
%{
  Input: 
    - params1 vector in R7, hard coded so position [beta p w v a d n]
    - params2 vector in R2, hard coded so position [X0 T0]
    - params3 vector in R5, hard coded so position [q1, q2, E0, L0, R0]
    - ImmigrationRate, which gives annual immigration rate, vector in R^N

  Output:
    - XELTR, an N+1 by 5 matrix where each row is a time
    - TBestimate, vector in R^N+1 with calculated incidence
%}

%% Unpack parameters

 beta = params1(1);
 p = params1(2);
 w = params1(3);
 v = params1(4);
 a = params1(5);
 d = params1(6);
 n= params1(7);
 
 q1 = params3(1);
 q2 = params3(2);
 
 %% Initialize outputs
 Nt = length(ImmigrationRate);
 
 XELTR = zeros(Nt+1, 5);
  X0 = params2(1);
  E0 = params3(3);
  L0 = params3(4);
  T0 = params2(2);
  R0 = params3(5);
  
 XELTR0 = [X0 E0 L0 T0 R0]; % this will be replaced soon
 XELTR(1,:) = XELTR0;
 
 
 TBestimate = zeros(Nt+1,1);
 TBestimate(1) = getTBIncidenceRate(XELTR0, p, w, v);
 

%% Calculate TB Incididence 

    
   
    XELTR0local = XELTR0;
    for i = 1:Nt

        pi = ImmigrationRate(i);
        
        % y = (x,e,l,t,r) is solution to system of odes.
        [t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n), [0,1], XELTR0local); % [0,1] means 1 year forward

        % update outputs
        XELTR0local = y(end,:);
        XELTR(i+1,:) = XELTR0local

        Tri = getTBIncidenceRate(XELTR0local, p, w, v);
        TBestimate(i+1) = Tri(1,1);
    end


%% 
function dydt = odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n)

dydt = zeros(4,1);
dydt(1) = (1-q1-q2)*pi - beta*y(1)*y(4) - n*y(1);
dydt(2) = q1*pi + beta*y(1)*y(4)-(n+w)*y(2);
dydt(3) = q2*pi + (1-p)*w*y(2) - (n+v)*y(3);
dydt(4) = p*w*y(2) + v*y(3) - (n+a+d)*y(4);
dydt(5) = d*y(4) - n*y(5);
end

%% calculate the TB Incidence Rate
function Tri = getTBIncidenceRate(y, p, w, v)
    Xi = y(:,1);
    Ei = y(:,2);
    Li = y(:,3);
    Ti = y(:,4);
    Ri = y(:,5);

    FPi = Xi + Ei + Li + Ti + Ri; %Foreign Population
    
    Tri = p*w*Ei + v*Li; %Page 10(704) of the York paper
    Tri = Tri * 100000./FPi;
end
end 