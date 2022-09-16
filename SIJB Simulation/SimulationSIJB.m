% SIJB Simulation using ODE45
% 
% Simulation of two diseases interacting.
%  Pre-conditions:
% Post-conditions:

t0 = 0;
tfinal = 3000;
N0 = 1000; % initial population size in unit persons

lambda = 0.01*N0;   % Birth per unit time
tau_I = 0.005/N0;   % Infectivity of I ?
tau_J = 0.0025/N0;  % Infectivity of J ?
gamma_I = 0.003;    % Recovered rate of I
gamma_J = 0.0015;   % Recovered rate of J
delta_S = 0.01;     % Proportion of death in S per unit time
delta_I = delta_S;  % Proportion of death in I per unit time
delta_J = delta_S;  % Proportion of death in J per unit time
delta_B = delta_S;  % Proportion of death in B per unit time

y0 = N0 * [0.7; 0.2; 0.1; 0]; % initial conditions (proportions of S, I, J, B)

SIJB = @(t,y) [lambda - tau_I*y(1).*(y(2) + y(4)) - tau_J * y(1) .* (y(3) + y(4)) + gamma_I*y(2) + gamma_J*y(3) - delta_S * y(1); 
    tau_I * y(1) .* (y(2) + y(4)) - gamma_I*y(2) - tau_J*y(2).*(y(3) + y(4)) + gamma_J * y(4) - delta_I * y(2);
    tau_J * y(1) .* (y(3) + y(4)) - gamma_J*y(3) - tau_I*y(3).*(y(2) + y(4)) + gamma_I * y(4) - delta_J * y(3);
    tau_J * y(2) .* (y(3) + y(4)) + tau_I * y(3) .* (y(2) + y(4)) - gamma_I * y(4) - gamma_J * y(4) - delta_B * y(4)];

[t,y] = ode45(SIJB, [t0 tfinal], y0);
N = sum(y,2); % population size

baseR_I = tau_I * lambda / delta_S / (gamma_I + delta_I);
baseR_J = tau_J * lambda / delta_S / (gamma_J + delta_J);

fprintf('Base Reproduction Number: (R0_I, R0_J < 1 for stability)\n');
fprintf(['R0_I: ',num2str(baseR_I),'\n']);
fprintf(['R0_J: ',num2str(baseR_J),'\n']);

figure
hold on
plot(t,y,'lineWidth',2)
%plot(t,N,'LineWidth',2)
title('Disease-Free Equilibrium')
xlabel('t')
ylabel('Population')
legend('S','I','J','B','location','eastoutside')
hold off
% print('ODE Disease-Free Equilibrium', '-dpng');


