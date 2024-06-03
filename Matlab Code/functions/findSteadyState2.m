%% find steady state 
function [y, k] = findSteadyState2(bioParameters, initialConditions, pi)
% input: initialConditions are [TP0, E0, L0, T0, R0]


 beta = bioParameters(1); %TB infectivity
 p = bioParameters(2); % ~probability someone in E goes straight into E; pi in Guo-Wu
 w = bioParameters(3); % period of time new infectee considered E rather than L
 v = bioParameters(4); % rate people in L develop TB
 a = bioParameters(5); % TB death rate for people in T
 d = bioParameters(6); % non-TB death rate
 n = bioParameters(7); % recovery rate; delta in Guo-Wu
 u = bioParameters(8); % percentage people in R develop TB relative to L
 q1 = bioParameters(9); % percent of immigrants into E
 q2 = bioParameters(10); % percent of immigrants into L
 q3 = bioParameters(11); % percent of immigrants into R
 
 %% Initialize outputs
 % Nt = length(ImmigrationRate);
 
 % XELTR = zeros(Nt+1, 5);
 
 X0 = initialConditions(1);
  E0 = initialConditions(2);
  L0 = initialConditions(3);
  T0 = initialConditions(4);
  R0 = initialConditions(5);
  

  y0 =[X0 E0 L0 T0 R0];


%while loop parameters
tstep = 50;
TOL = 1e-4;
stopflag = 0;
maxiters = 10000;
k = 1;

y = y0;
% ybar = norm(y0);

while k < maxiters && stopflag == 0
    % compute error
    k = k+1;
   
    % solve ode
    [t,yv] = ode23(@odefcn, [0, tstep],y);

    ynew = yv(end,:);

    if norm(ynew-y)/norm(y)<TOL
        stopflag = 1;
    end
    y = ynew;
    % yfinal = y;
end

if k == maxiters
    display('after max steps, did not terminate')
end
% [t,yv] = ode23(@odefcn, tspan, y0);

% yfinal = ynew;

%% Guo-Wu
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