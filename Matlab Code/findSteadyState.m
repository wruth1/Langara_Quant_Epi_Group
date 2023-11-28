%% find steady state 
function yfinal = findSteadyState(bioParameters, initialConditions, pi)

 beta = bioParameters(1); %TB infectivity
 p = bioParameters(2); % ~probability someone in E goes straight into E; pi in Guo-Wu
 w = bioParameters(3); % period of time new infectee considered E rather than L
 v = bioParameters(4); % rate people in L develop TB
 a = bioParameters(5); % TB death rate for people in T
 d = bioParameters(6); % non-TB death rate
 n = bioParameters(7); % recovery rate; delta in Guo-Wu
 q1 = bioParameters(8); % percent of immigrants with early latent TB
 q2 = bioParameters(9); % percent of immigrants with early latent TB
 
 %% Initialize outputs
 % Nt = length(ImmigrationRate);
 
 % XELTR = zeros(Nt+1, 5);
 
  E0 = initialConditions(2);
  L0 = initialConditions(3);
  T0 = initialConditions(4);
  R0 = initialConditions(5);
  X0 = initialConditions(1) - (E0+L0+T0+R0);

  y0 =[X0 E0 L0 T0 R0];

tfinal = 2000;

tspan = 0:tfinal;

[t,yv] = ode23(@odefcn, tspan, y0);

yfinal = yv(end,:);



function dydt = odefcn(t, y  )

    %pi = ImmigrationRate(t);
 
dydt = zeros(5,1);
dydt(1) = (1-q1-q2)*pi - beta*y(1)*y(4) - n*y(1);
dydt(2) = q1*pi + beta*y(1)*y(4)-(n+w)*y(2);
dydt(3) = q2*pi + (1-p)*w*y(2) - (n+v)*y(3);
dydt(4) = p*w*y(2) + v*y(3) - (n+a+d)*y(4);
dydt(5) = d*y(4) - n*y(5);
end

end