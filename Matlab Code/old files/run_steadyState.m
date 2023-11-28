 %optimization toolbox

%% Parameters Setup.  Hard code parameter values.

%Actual Pi, annual new immigration to Canada every year 2010-2020
ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

%Actual Pi, annual new immigration to Canada every year 2001-2020
%Actual_Pi = [223840 199170 239083 244578 254374 238125 249622 245289 270581 259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];


ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020


%% bioParameters

%Parameters from York Paper
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada.

p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
w = 0.4;               % (q1+q2) Pro portion of all LTBI immigrants pass through latent stage in the first 2.5 years
% v = 0.0002;            % The rate of slow progression to active TB due to reactivation
v = p*w/15;    %15x more likely to develop
% According to Standards Canada [cite \url{https://www.canada.ca/en/public-health/services/infectious-diseases/canadian-tuberculosis-standards-7th-edition/edition-18.html}], ``newly" infected patients [$\leq2$ years from infection] are $15\times$ more likely to develop active TB than ``regular" [people with no known risk factor]. 
a = 0.06;              % TB-caused death rate
d = 0.8;               % Constant rate of recovery by nature or treatment
n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039
%q1 =     0.0674;
%q2 =     0.1069;
q1 = 0.03;
q2 = 0.4;



% 2010 Initial Condition from TB Surveilance by Stat Can
TP0 = 6775765 - ReportedImmigration(1) ; % 6775765 2011 census  % https://en.wikipedia.org/wiki/Canada_immigration_statistics#2011_census
T0 = 1054; % reported in 2010
R0 = 0; 
% TP = total FB population

% initial % of populations from Guo Wu
E0 = 0.001735006* TP0;
L0 = 0.21218547 * TP0; 




%%  call and run solveGuoWu; unoptimized

bioParameters = [beta p w v a d n q1 q2];
initialConditions = [TP0, E0, L0, T0, R0 ];

PI = ReportedImmigration(1);

y = findSteadyState(bioParameters, initialConditions,pi);

yrescale = y/sum(y)*TP0;

incidence = getTBIncidenceRate(yrescale,p,w,v);
%%
