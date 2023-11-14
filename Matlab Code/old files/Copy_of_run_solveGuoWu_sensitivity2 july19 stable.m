%optimization toolbox%% Parameters Setup

%Parameters from York Paper
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada.
w = 0.4;               % (q1+q2) Pro portion of all LTBI immigrants pass through latent stage in the first 2.5 years
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
v = 0.0002;            % The rate of slow progression to active TB due to reactivation
a = 0.06;              % TB-caused death rate
d = 0.8;               % Constant rate of recovery by nature or treatment
n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039



%Actual Pi, annual new immigration to Canada every year 2010-2020
ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

%Actual Pi, annual new immigration to Canada every year 2001-2020
%Actual_Pi = [223840 199170 239083 244578 254374 238125 249622 245289 270581 259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];


% 2010 Initial Conditions from TB Surveilance by Stat Can
X0 = 5874633;
T0 = 1054;
R0 = 0;

%INIT_L0 = (INIT_ACTUAL_TB_RATE*(y0(1) + y0(2) + y0(4) + y0(5)) - 100000*p*w*y0(2))/(100000*v-INIT_ACTUAL_TB_RATE);

%y0 = [X0, INIT_E0, INIT_L0, T0, R0];

%INIT_ACTUAL_TB_RATE = ReportedTB(1);

ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020



%% call and run solveGuoWu

params1 = [beta p w v a d n];
params2 = [X0, T0, R0];

ELratio=10;

q1 = 0.000645;
q2 = ELratio*q1;
E0 = 35000;
%L0 = 7.310031644067797e+06;
L0 = ELratio*E0;
R0 = 0;
params3 = [q1, q2, E0, L0];

[XELTR, EstimatedTB] = solveGuoWu(params1, params2, params3, ReportedImmigration);

% this XELTR is ~ half hour of me fiddling with parameters

%% Now minimize using fminsearch

x0 = params3;

options = optimset('PlotFcns',@optimplotfval);
%xmin_FMS = fminsearch(@IncidenceError, x0) ; % about 10 seconds to run
%save('xmin_FMS.mat');

%https://www.mathworks.com/help/optim/ug/fmincon.html



load('xmin_FMS.mat')


[XELTR2, EstimatedTB2] = solveGuoWu(params1, params2, xmin_FMS, ReportedImmigration);

% XELTR2 is the "best"
%% Minimize using fmincon



%q1start = 0.0544;
%q2start = 0.3224;


%  $e_* = 22533$ and $l_* = 3016687$
%x0(3) = 22533;  %E0
%x0(4) = 3016687;




% intialize fmincon matrix constraints

% row1 is
% row1 is q1 + q2 <= 1
% row2 is E0 + L0 <= total population (we assume R0 is 0)

A  = [ 1 1 0 0];
b = [1  ];
% define equality constraints
Aeq = [];
beq = [];
% bounds
lb = zeros(size(x0));
ub = [Inf Inf Inf Inf];

%q1start =     0.0674;
%q2start =     0.1069;










% initialize time grid
Nt = length(ReportedTB);
Nt = Nt-1;
TimeGrid = 1:Nt;
TimeGrid0 = [0, TimeGrid];

% plot



q1list = [ 0.07 ];
q2list = [ 0.08  0.3 ];

%q1list = [ 0.005 ];

%q2list = [ 0.03 ];

numq1 = length(q1list);
numq2 = length(q2list);
% Calculate total population

% TP = total population
TP0 = 6775765 - ReportedImmigration(1) ; % 6775765 2011 census  % https://en.wikipedia.org/wiki/Canada_immigration_statistics#2011_census
% Nnl0 is the upper bound of initial $ of non-latent-TB People
Nnl0 = 6775765 - ReportedImmigration(1) - sum(params2) ; % upper bound initial # of non-E0 or L0
Pnl0 = Nnl0/TP0;

% initialize E0 and L0
%ELscale = 0.2; % applied to both E0 and L0
%ELratio = 0.1; % L0 uses (1-this)
%E0start = ELratio* Nnl0 * ELscale;
%L0start = (1-ELratio)* Nnl0 * ELscale;

E0list =[ 5e4 ] ;
L0list =[1e5 4e5 7e5 1e6] ;

numE0 = length(E0list);
numL0 = length(L0list);



% simulate the rest
xoptimal = zeros(numq1, numq2, numE0, numL0,4);

% close all;
%  set(0,'DefaultFigureWindowStyle','docked')
for ie0=1:length(E0list)
     E0start = E0list(ie0);
    for il0=1:length(L0list)
        L0start = L0list(il0);
        
%         figure
%         % plot reported incidence
%         plot(TimeGrid0,ReportedTB(1:end),'LineWidth',1, 'DisplayName','Reported TB Incidence');
%         %title(strcat('Canadian Foreign-born Population TB Incidence Rate, E0 L0 =', num2str(E0start), ',', num2str(L0start)))
%         title(strcat('TB Incidence Rate, E0 L0 =', num2str(E0start), ',', num2str(L0start)))
%         xlabel('Years After 2010')
%         ylabel('New Active TB Cases (per 100,000 population)')
%         
%         hold on
        
        for iq1=1:length(q1list)
            q1start = q1list(iq1);
            for iq2=1:length(q2list)
                
                q2start = q2list(iq2);
                
                
                x0(1) = q1start;
                x0(2) = q2start;
                x0(3) = E0start;
                x0(4) = L0start;
                %                x0(3) = E0start;
                %                x0(4) = L0start;
                
                
                xmin = fmincon(@IncidenceError, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run
                
                xoptimal(iq1,iq2,ie0,il0,:) = xmin;
                
%                 x0now = zeros(1,4);
%                 x0now(1) = xoptimal(iq1, iq2, ie0, il0, 1);
%                 x0now(2) = xoptimal(iq1, iq2, ie0, il0, 2);
%                 x0now(3) = xoptimal(iq1, iq2, ie0, il0, 3);
%                 x0now(4) = xoptimal(iq1, iq2, ie0, il0, 4);
%                 
%                 [XELTRnow, EstimatedTBnow] = solveGuoWu(params1, params2, x0now, ReportedImmigration);
%                 
%                 plot(TimeGrid0,EstimatedTBnow(1:end-1), 'DisplayName', strcat('FMC q1 q2 =', num2str(q1start), ',' , num2str(q2start) , 'x=', num2str(xmin,3) ))
%                 legend('show','Location','northwest');
%                 
                
            end
        end
        
       % hold off
    end
end

%% Plot optimal trajectories

close all;
 set(0,'DefaultFigureWindowStyle','docked')
for ie0=1:length(E0list)
     E0start = E0list(ie0);
    for il0=1:length(L0list)
        L0start = L0list(il0);
        
        figure
        % plot reported incidence
        plot(TimeGrid0,ReportedTB(1:end),'LineWidth',1, 'DisplayName','Reported TB Incidence');
        %title(strcat('Canadian Foreign-born Population TB Incidence Rate, E0 L0 =', num2str(E0start), ',', num2str(L0start)))
        title(strcat('TB Incidence Rate, E0 L0 =', num2str(E0start), ',', num2str(L0start)))
        xlabel('Years After 2010')
        ylabel('New Active TB Cases (per 100,000 population)')
        
        hold on
        
        for iq1=1:length(q1list)
            q1start = q1list(iq1);
            for iq2=1:length(q2list)
                
                q2start = q2list(iq2);
                
                
                x0(1) = q1start;
                x0(2) = q2start;
                x0(3) = E0start;
                x0(4) = L0start;
                %                x0(3) = E0start;
                %                x0(4) = L0start;
                

                
                x0now = zeros(1,4);
                x0now(1) = xoptimal(iq1, iq2, ie0, il0, 1);
                x0now(2) = xoptimal(iq1, iq2, ie0, il0, 2);
                x0now(3) = xoptimal(iq1, iq2, ie0, il0, 3);
                x0now(4) = xoptimal(iq1, iq2, ie0, il0, 4);
                
                [XELTRnow, EstimatedTBnow] = solveGuoWu(params1, params2, x0now, ReportedImmigration);
                
                plot(TimeGrid0,EstimatedTBnow(1:end-1), 'DisplayName', strcat('FMC q1 q2 =', num2str(q1start), ',' , num2str(q2start) , '; x=', num2str(x0now,3) ))
                legend('show','Location','northwest');
                
                
            end
        end
        
        hold off
    end
end


%% 

%xmin = fmincon(@IncidenceError, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run






%% Functions

%% Setup the ODE System
%P.4 of York Paper

% X <-> y(1)
% E <-> y(2)
% L <-> y(3)
% T <-> y(4)
% R <-> y(5)
function err = IncidenceError2(y,q1,E0)
    err = IncidenceError([q1,y(1),E0, y(2)]); 
end 

function err = IncidenceError(x)
% input: x is in R5, q1, q2, E0, L0, R0

q1 = x(1);
q2 = x(2);
E0 = x(3);
L0 = x(4);

% Parameters Setup

beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada.
w = 0.4;               % (q1+q2) Pro portion of all LTBI immigrants pass through latent stage in the first 2.5 years
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
v = 0.0002;            % The rate of slow progression to active TB due to reactivation
a = 0.06;              % TB-caused death rate
d = 0.8;               % Constant rate of recovery by nature or treatment
n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039
ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];
X0 = 5874633;
T0 = 1054;
TFP2011 = 6775765; % from wikipedia
R0 = TFP2011 - ReportedImmigration(1) - X0 - T0 - E0 - L0 ;
ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020
params1 = [beta p w v a d n];
params2 = [X0, T0, R0];

[~, EstimatedTB] = solveGuoWu(params1, params2, [q1, q2, E0, L0, R0], ReportedImmigration);

err = norm(EstimatedTB-ReportedTB);
end


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

%% To-dos
% Make pi into a smooth function
% vary q1 and q2
% E0 and L0
