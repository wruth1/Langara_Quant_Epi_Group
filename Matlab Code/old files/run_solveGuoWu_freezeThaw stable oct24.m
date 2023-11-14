%optimization toolbox

%% Parameters Setup.  Hard code parameter values.

% Observed Data

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
v = 0.0002;            % The rate of slow progression to active TB due to reactivation
% v = p*w/15;    %15x more likely to develop

a = 0.06;              % TB-caused death rate
d = 0.8;               % Constant rate of recovery by nature or treatment
n = 0.039;             %Natural removal rate, dX = dE = dL = dT = dR = 0.039
%q1 =     0.0674;
%q2 =     0.1069;
q1 = 0.03;
q2 = 0.37;



% Initial Population Conditions

%{
Total population is known.

T0 is reported and thus known.

E0, L0, R0 are fixed at the start, and then X0 is calculated so total pop
is correct.
%}

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
initialConditions = [TP0, E0, L0, T0, R0];


[XELTR, EstimatedIncidence] = solveGuoWu(bioParameters, initialConditions, ReportedImmigration);


% plot popluation vs time


t0 = 2010;
numYear = length(ReportedTB);
yearGrid = t0:t0+numYear;

figure 
plotPopulationVsTime(yearGrid,XELTR,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence, ReportedTB )
saveas(gcf,'PopVsTime_unoptimized.png')
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
[XELTR, EstimatedIncidence] = solveGuoWu(BP2, IC2, ReportedImmigration);

 

figure 
plotPopulationVsTime(yearGrid,XELTR,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence, ReportedTB )
% compareIncidence(xmin, bioParameters, initialConditions, ReportedImmigration, ReportedTB, yearGrid );
saveas(gcf,'PopVsTime_optimized.png')

% compare x

xcomparison = zeros(length(xmin5), 2);
xcomparison(:,1) = x0;
xcomparison(:,2) = xmin5;
display('q1, q2, E0, L0, R0:')
xcomparison 

%% Now thaw on L0

ReportedImmigration0 = ReportedImmigration(1);
ReportedTB0 = ReportedTB(1);

f1=@(x)IncidenceError1(x,BP2, IC2, ReportedImmigration, ReportedTB);

% minimizer starting value
x0 = xmin5(4); % load L0 from 5D solve



[xmin1,fval1,exitflag,output] = fminsearch(f1, x0 ) ; % about 10 seconds to run

BP3 = BP2;

IC3 = IC2;
IC3(3) = xmin1; % update L0

[XELTR, EstimatedIncidence] = solveGuoWu(BP3, IC3, ReportedImmigration);



figure 
plotPopulationVsTime(yearGrid,XELTR,ReportedImmigration);

subplot(2,3,6)
plotIncidence(yearGrid, EstimatedIncidence, ReportedTB )
% compareIncidence(xmin, bioParameters, initialConditions, ReportedImmigration, ReportedTB, yearGrid );
% saveas(gcf,'PopVsTime_optimized.png')

%% SENSITIVITY ANALYSIS
%% create lists

beta_list = [beta];
p_list = [p];
w_list = [w/2,w,2*w];

v_list = [v/4, v/2, v, 2*v];            % The rate of slow progression to active TB due to reactivation
a_list = [a];              % TB-caused death rate
d_list = [d];               % Constant rate of recovery by nature or treatment
n_list = [n]; 

q1_list = [q1];
q2_list = [q2];

TP0_list = [TP0];
E0_list = [E0];
L0_list = [L0/2,L0,L0*2];
T0_list = [T0];
R0_list = [0, 1E3, 1E5];

paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, q1_list, q2_list, TP0_list,E0_list, L0_list, T0_list,R0_list};


%% Big Data Structure XX
%{
XX stores all the data.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st column is vector |BP|+|IC|, starting parameters
- 2nd column is vector in R4, optimal [q1 q2 E0 L0]
- 3rd column is hessian (output of fmincon)
- 4th column is matrix, population XELTR vs time (using optimized params)
- 5th column is TB incidence (using optimized params)
- 6th column is error of TB incidence using optimized params
%}

% use Will's expand_grid.m to create every combination of parameters
paramgrid = expand_grid(paramcell{:}); 

% initialize XX
NumSims = size(paramgrid,1);
XX = cell(NumSims, 6);



% intialize fmincon matrix constraints
% row1 is q1 + q2 <= 1
A  = [1 1 0 0];
b = [1];
% define equality constraints
Aeq = [];
beq = [];
% bounds
lb = zeros(size(x0));
ub = [Inf Inf Inf Inf];


% fill in XX
for i = 1:NumSims

    % get parameters, save, load
    paramsi = paramgrid(i,:);
    XX{i,1} = paramsi;
    BPi = paramsi(1:9);
    ICi = paramsi(10:14);
    
    % setup optimizer inputs
    f5=@(x)IncidenceError(x,BPi, ICi, ReportedImmigration, ReportedTB);
    x0 = [BPi(8) BPi(9) ICi(2) ICi(3)]; %q1 q2 E0 L0

    % optimize, store results
    [xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f5, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run
    XX{i,2} = xmin5;
    XX{i,3} = hessian;
    
    % compute other XX output; pop vs time and incidence error

    % update parameters to optimal ones
    BPi(8) = xmin5(1); %q1
    BPi(9) = xmin5(2); %q2
    ICi(2) = xmin5(3); %E0
    ICi(3) = xmin5(4); %L0

    [XELTRi, TBIncidencei] = solveGuoWu(BPi, ICi, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = TBIncidencei;
    XX{i,6} = fval5;

end


save('XX.mat')
save('paramcell.mat')
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

%
function err = IncidenceError1(x, bioParameters, initialConditions, ImmigrationRate, ReportedTB)
% input: x is in R1, L0
% bioParameters and initialConditions don't have to be exactly correct;
% they will be updated by x

% load bioParameters
localBP = bioParameters;


%load initialConditions
localIC = initialConditions;
localIC(3)=x(1); % L0

[~, EstimatedTB] = solveGuoWu(localBP, localIC, ImmigrationRate);

   
err = norm((EstimatedTB-ReportedTB));
%err = norm((EstimatedTB-ReportedTB)./ReportedTB);
end

