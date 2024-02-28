
%% Plot populations vs time 

addpath('functions\')
addpath('figure creation\')


%% Parameters Setup.  Hard code parameter values.

% Observed Data

%Actual Pi, annual new immigration to Canada every year 2010-2020
ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];
%Actual Pi, annual new immigration to Canada every year 2001-2020
%Actual_Pi = [223840 199170 239083 244578 254374 238125 249622 245289 270581 259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020

ReportedT = [1054 1108 1112 1153 1110 1178 1231 1319 1315 1427 1303];



% bioParameters

%Parameters from York Paper
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada.
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
w = 0.4;               % Proportion of all LTBI immigrants pass through latent stage in the first 2.5 years
v = 0.0002;            % The rate of slow progression to active TB due to reactivation
% v = p*w/15;    %15x more likely to develop
a = 0.06;              % TB-caused death rate
d = 0.8;               % Constant rate of recovery by nature or treatment
n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039
%q1 =     0.0674;
%q2 =     0.1069;
q1 = 0.03;
q2 = 0.37;

bioParameters = [beta p w v a d n q1 q2];
initialConditions = [TP0, E0, L0, T0, R0 ];

%% Guo-Wu Parameters naked

[XELTR_naked, EstimatedIncidence_naked] = solveGuoWu(bioParameters, initialConditions, ReportedImmigration);


% plot popluation vs time


t0 = 2010;
numYear = length(ReportedTB);
yearGrid = t0:t0+numYear;

figure (1)
title('Guo-Wu Parameters naked')
plotPopulationVsTime(yearGrid,XELTR_naked,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence_naked, ReportedTB )
% saveas(gcf,'PopVsTime_unoptimized.png')

%% Minimize using fmincon

f5=@(x)IncidenceError5(x,bioParameters, initialConditions, ReportedImmigration, ReportedTB);

% minimizer starting value
x0 = [q1 q2 E0 L0 R0];

% intialize fmincon matrix constraints

% row1 is q1 + q2 <= 1
A  = [1 1 0 0 0];
b = [1];
% define equality constraints
Aeq = [];
beq = [];
% bounds
lb = zeros(size(x0));
ub = lb + Inf;

%
%xmin = fmincon(f, x0) ; % about 10 seconds to run
[xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f5, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run


% calculate 
[BP2, IC2] = updateParameters5(bioParameters, initialConditions, xmin5);
[XELTR_opt1, EstimatedIncidence_opt1] = solveGuoWu(BP2, IC2, ReportedImmigration);

 

figure(2)
plotPopulationVsTime(yearGrid,XELTR_opt1,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence_opt1, ReportedTB )
% compareIncidence(xmin, bioParameters, initialConditions, ReportedImmigration, ReportedTB, yearGrid );
saveas(gcf,'PopVsTime_optimized.png')

%% Steady state

[ysteady,k] = findSteadyState(bioParameters, initialConditions, ReportedImmigration(1));

% rescale 
ysteady = ysteady*TP0/sum(ysteady);

% calculate XELTR

IC3 = ysteady;
IC3(1) = TP0;

[XELTR_steady, EstimatedIncidence_steady] = solveGuoWu(bioParameters, IC3, ReportedImmigration);

figure(3)
plotPopulationVsTime(yearGrid,XELTR_steady,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence_steady, ReportedTB )
% compareIncidence(xmin, bioParameters, initialConditions, ReportedImmigration, ReportedTB, yearGrid );
% saveas(gcf,'PopVsTime_optimized.png')

%% Minimize using fmincon, with steadystate


% IC3 from steadystate
f5b=@(x)IncidenceError5(x,bioParameters, IC3, ReportedImmigration, ReportedTB);

% minimizer starting value
x1 = [q1 q2 IC3(2) IC3(3) IC3(5)]; %q1 q2 E0 L0 R0

% intialize fmincon matrix constraints

% row1 is q1 + q2 <= 1
A  = [1 1 0 0 0];
b = [1];
% define equality constraints
Aeq = [];
beq = [];
% bounds
lb = zeros(size(x0));
ub = lb + Inf;

%
%xmin = fmincon(f, x0) ; % about 10 seconds to run
[xmin5b,fval5b,exitflag,output,lambda,grad,hessian] = fmincon(f5b, x1 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run


% calculate 
[BP2b, IC2b] = updateParameters5(bioParameters, IC3, xmin5b);
[XELTR_opt2, EstimatedIncidence_opt2] = solveGuoWu(BP2b, IC2b, ReportedImmigration);

 

figure(4)
plotPopulationVsTime(yearGrid,XELTR_opt2,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence_opt2, ReportedTB )
% compareIncidence(xmin, bioParameters, initialConditions, ReportedImmigration, ReportedTB, yearGrid );
saveas(gcf,'PopVsTime_optimized.png')

%% 

names = {'X','E','L','T','R'};

% figure(5)
figure('units','normalized','outerposition',[0 0 1 1]); clf;
for k=1:5
    subplot(2,3,k)
    plot(yearGrid,XELTR_naked(:,k),':','DisplayName','Guo-Wu'); grid on;
    hold on
    plot(yearGrid,XELTR_opt1(:,k),'.-','DisplayName','fmincon with Guo-Wu');
    plot(yearGrid,XELTR_steady(:,k),'*-','DisplayName','steady state');
    plot(yearGrid,XELTR_opt2(:,k),'o-','DisplayName','fmincon with steady state');
    title([names{k},' vs Year'])
    xlabel('Year')
    ylabel(names{k})
    if k==4
        plot(yearGrid(1:end-1), ReportedT,'DisplayName','Reported','LineWidth',3)
    end
    % legend('Location','West')
    % hold off
end
subplot(2,3,6)
    plot(yearGrid,EstimatedIncidence_naked,':','DisplayName','Guo-Wu'); grid on;
    hold on
    plot(yearGrid,EstimatedIncidence_opt1,'.-','DisplayName','fmincon with Guo-Wu');
    plot(yearGrid,EstimatedIncidence_steady,'*-','DisplayName','steady state');
    plot(yearGrid,EstimatedIncidence_opt2,'o-','DisplayName','fmincon with steady state');
    plot(yearGrid(1:end-1), ReportedTB,'DisplayName','Reported','LineWidth',3)
    title(['Incidence vs Year'])
    xlabel('Year')
    ylabel('Incidence')
    legend('Location','SouthEast')
    hold off
saveas(gcf,'PopVsTime_various.png')
% exportgraphics(gca,'PopVsTime_various.png','Resolution',300)
%%

%% Functions

function [BP, IC] = updateParameters5(BP0, IC0, x5)
    BP = BP0;
    IC = IC0;

    BP(8) = x5(1); %q1
    BP(9) = x5(2); %q2
    IC(2) = x5(3); %E0
    IC(3) = x5(4); %L0
    IC(5) = x5(5); %R0

    

end
%% IncidenceError, function to optmmize put into fmincon
function err = IncidenceError5(x, bioParameters, initialConditions, ImmigrationRate, ReportedTB)
% input: x is in R5, q1, q2, E0, L0, R0

% load bioParameters
localBP = bioParameters;
localBP(8)=x(1); % q1
localBP(9)=x(2); % q2

%load initialConditions
localIC = initialConditions;
localIC(2)=x(3); % E0
localIC(3)=x(4); % L0
localIC(5)=x(5);

[~, EstimatedTB] = solveGuoWu(localBP, localIC, ImmigrationRate);

   
err = norm((EstimatedTB-ReportedTB));
%err = norm((EstimatedTB-ReportedTB)./ReportedTB);
end


