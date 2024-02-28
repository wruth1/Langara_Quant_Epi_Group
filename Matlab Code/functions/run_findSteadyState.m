
%% Demonstrate steady state

addpath('functions\')

%% Parameters Setup.  Hard code parameter values.

% Observed Data

%Actual Pi, annual new immigration to Canada every year 2010-2020
ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];
%Actual Pi, annual new immigration to Canada every year 2001-2020
%Actual_Pi = [223840 199170 239083 244578 254374 238125 249622 245289 270581 259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020

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


% 2010 Initial Condition from TB Surveilance by Stat Can
TP0 = 6775765 - ReportedImmigration(1) ; % 6775765 2011 census  % https://en.wikipedia.org/wiki/Canada_immigration_statistics#2011_census
T0 = 1054; % reported in 2010
% R0 = 50000; 
R0 = 0;
% TP = total FB population

% initial % of populations from Guo Wu
E0 = 0.001735006* TP0;
L0 = 0.21218547 * TP0; 

% y0 = yfinal;

% bioParameters = [beta p w v a d n q1 q2];
% initialConditions = [TP0, E0, L0, T0, R0];

initialConditions = [TP0, E0, L0, T0, R0];

[yfinal,k] = findSteadyState(bioParameters, initialConditions, ReportedImmigration(1));

display(['Estimated steady state total population is: ', num2str(sum(yfinal)), '; reported is ', num2str(TP0)]);


estimated_incidence = getTBIncidenceRate(yfinal, p, w, v);

display(['Estimated steady state incidence is: ', num2str(estimated_incidence), '; reported is ',num2str(ReportedTB(1))]);

display(['Estimated steady state T is: ', num2str(yfinal(4)), '; reported is ',num2str(1054)]);