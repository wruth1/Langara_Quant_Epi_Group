%% Parameters Setup
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada. 
pi = 223840;           % Average number of annual new immigrants to Canada     
w = 0.4;               % (q1+q2) Proportion of all LTBI immigrants pass through latent stage in the first 2.5 years
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
v = 0.0002;            % The rate of slow progression to active TB due to reactivation 
a = 0.06;              % TB-caused death rate 
d = 0.8;               % Constant rate of recovery by nature or treatment 
q1 = 0.03;             % Percentages of early latent (high risk) new immigrants to develop TB
q2 = 0.37;             % Percentages of late latent latent (low risk) new immigrants to develop TB
tspan = [0:1:100];     % Time span

n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039

y0 = [4431746, 9784, 1196551, 1094, 0]; %[X0, E0, L0, T0, R0] Initial Value
% y0(1) = 3000000;


%% TB incidence rate varying with q1
% Note that q1 > 0, q2 = 0 (fixed) in this simulation

Tr = [];
Q1_init = 0;
Q1_final = 0.25;
Q1_step = 0.01;
Q1 = [Q1_init:Q1_step:Q1_final];

% Version 1
% for q1 = Q1
%     [t,y] = ode23(@(t,y) odefcn(t, y, beta,pi, p, w, v, a ,d, q1, 0, n), tspan, y0);
% 
%     Xi = y(:,1);
%     Ei = y(:,2);
%     Li = y(:,3);
%     Ti = y(:,4);
% 
%     %Foreign Population, FP
%     FPi = Xi + Ei + Li + Ti;
% %     FPi = sum(y0);
% 
%     %Tri = Ti*100000/(y0(1));
%     Tri = p*w*Ei + v*Li; %Page 10(704) of the York paper
%     Tri = 100000 * (Tri./(FPi));
% 
%     if q1 == Q1_init
%         Tr = [Tri];  %Tr stands for TB incidence rate
%     else
%         Tr = [Tr Tri];
%     end
% end

%Version2
for q1 = Q1
    [t,y] = ode23(@(t,y) odefcn(t, y, beta,pi, p, w, v, a ,d, q1, 0, n), tspan, y0);
    Tri = getTBIncidenceRate(y, p, w, v);
    Tr = [Tr Tri];
end



%% Plot TB incidence rate varying with q1
figure(2)
[t1, q1] = meshgrid(t,Q1);
surf(t1,q1,Tr')
colormap("parula")
title('Fig 3AL. TB incidence rate varying with q_1')
xlabel('time (year)')
ylabel('q_1')
zlabel('TB incidence rate per 100,000 population')
%% TB incidence rate varying with q2
% Note that q1 = 0(fixed), q2 > 0 in this simulation
Tr = [];
Q2_init = 0;
Q2_step = 0.01;
Q2_final = 0.4;

Q2 = [Q2_init:Q2_step:Q2_final];
for q2 = Q2
    [t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, 0, q2, n), tspan, y0);
    Tri = getTBIncidenceRate(y, p, w, v);
    Tr = [Tr Tri];

end
%% Plot TB incidence rate varying with q2
[t1, q2] = meshgrid(t,Q2);
surf(t1,q2,Tr')
title('Fig 3AR. TB incidence rate varying with q_2')
xlabel('time (year)')
ylabel('q_2')
zlabel('TB incidence rate per 100,000 population')

%% TB incidence rate varying with pi


Tr = [];
Pi_init = 111920;
Pi_final = 447680;
Pi_step = 10000;
Pi = [Pi_init:Pi_step:Pi_final];
for pi = Pi
    [t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n), tspan, y0);
    Tri = getTBIncidenceRate(y, p, w, v);
    Tr = [Tr Tri];
end
%% Plot TB incidence rate varying with pi
[t1, pi1] = meshgrid(t,Pi);
surf(t1,pi1,Tr')
colormap("winter")
title('TB incidence rate varying with new immigration rate \pi')
xlabel('time (year)')
ylabel('\pi')
zlabel('TB incidence rate per 100,000 population')

%% TB incidence rate varying with \beta

Tr = [];
Beta_init = 0*10^-7;
Beta_final = 2*10^-7;
Beta_step = 0.1*10^-7;
Beta = [Beta_init:Beta_step:Beta_final];
q1 = 0.03;
q2 = 0.37;
pi = 1.3*10^5; %Altered
for beta = Beta
    [t,y] = ode23(@(t,y) odefcn(t, y, beta, pi, p, w, v, a ,d, q1, q2, n), tspan, y0);
    Tri = getTBIncidenceRate(y, p, w, v);
    Tr = [Tr Tri];
end
%% Plot TB incidence rate varying with \beta
[t1, beta1] = meshgrid(t,Beta);
surf(t1,beta1,Tr')
colormap("default")
title('Fig 4b TB Incidence rate with varying transmission coefficient \beta')
xlabel('time (year)')
ylabel('\beta')
zlabel('TB incidence rate per 100,000 population')


%% TB incidence rate varying with q1 and q2
% Note that t = 101 (fixed) in this simulation


Q1_init = 0;
Q1_final = 0.25;
Q1_step = 0.01;
Q1 = [Q1_init:Q1_step:Q1_final];

Q2_init = 0;
Q2_final = 0.4;
Q2_step = 0.0125;
Q2 = [Q2_init:Q2_step:Q2_final];

Tr = zeros(length(Q2),length(Q1));

p = 0.05;

for q2 = Q2
    for q1 = Q1
        [t,y] = ode23(@(t,y) odefcn(t, y, beta,pi, p, w, v, a ,d, q1, q2, n), tspan, y0);
        Tri = getTBIncidenceRate(y, p, w, v);
        Tr(int32(q2/Q2_step+1),int32(q1/Q1_step+1)) = Tri(length(tspan),1);
    end
end




%% Plot TB incidence rate varying with q1 and q2
[q1, q2] = meshgrid(Q1,Q2);
surf(q1,q2,Tr)
colormap(hsv)
title('Fig 3BL. TB incidence rate varying with q_1 and q_2')
xlabel('q_1')
ylabel('q_2')
zlabel('TB incidence rate per 100,000 population')

%% Plot TB incidence rate varying with q1 and q2 isocline
[q1, q2] = meshgrid(Q1*100,Q2*100);
contour(q1, q2, Tr, 10)
xlabel('q_1')
ylabel('q_2')
colormap(turbo)
title('Fig 3BR. TB incidence rate relative to q_1 and q_2')
% caxis([0 52])
cb = colorbar('eastoutside');
% surf(q1,q2,Tr)


% zlabel('TB incidence rate per 100,000 population')


%% Setup the ODE System
%P.4 of York Paper

% X <-> y(1)
% E <-> y(2)
% L <-> y(3)
% T <-> y(4)
% R <-> y(5)


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
