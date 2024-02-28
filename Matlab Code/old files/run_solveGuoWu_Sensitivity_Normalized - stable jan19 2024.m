
%% SENSITIVITY ANALYSIS

addpath('functions\')
addpath('figure creation\')

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


% 2010 Initial Condition from TB Surveilance by Stat Can
TP0 = 6775765 - ReportedImmigration(1) ; % 6775765 2011 census  % https://en.wikipedia.org/wiki/Canada_immigration_statistics#2011_census
T0 = 1054; % reported in 2010
R0 = 0; 
% TP = total FB population

% initial % of populations from Guo Wu
E0 = 0.001735006* TP0;
L0 = 0.21218547 * TP0; 


% bioParameters = [beta p w v a d n q1 q2];
% initialConditions = [TP0, E0, L0, T0, R0];



%% create lists for sensitivity analysis

beta_list = [beta];
p_list = [p];

scale_param = 9/10; % assume <1
scale_list = scale_param.^[1,0,-1,-2];
w_list = w*scale_list;
v_list = v*scale_list;            % The rate of slow progression to active TB due to reactivation
% w_list = w;
% v_list = v;
a_list = [a];              % TB-caused death rate
d_list = [d];               % Constant rate of recovery by nature or treatment
n_list = [n]; 

q1_list = q1*scale_list;
q2_list = q2*scale_list;

TP0_list = [TP0];
E0_list = [E0];
L0_list = [L0];
T0_list = [T0];
R0_list = [1E6];

paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, q1_list, q2_list, TP0_list,E0_list, L0_list, T0_list,R0_list};

% save('Experiment_Parameters' )

load('normalization_minmax.mat')
normalization_minmax = normalize_minmax;

%%
%% Big Data Structure XX
%{
XX stores all the data.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st column is vector |BP|+|IC|, starting parameters
- 2nd column is vector in R4, optimal [q1 q2 E0 L0]
- 3rd column is output of fmincon, {exitflag, output, lambda, grad, hessian};
- 4th column is matrix, population XELTR vs time (using optimized params)
- 5th column is TB incidence (using optimized params)
- 6th column is error of TB incidence using optimized params
%}

% use Will's expand_grid.m to create every combination of parameters
paramgrid = expand_grid(paramcell{:}); 

% initialize XX
NumSims = size(paramgrid,1);
XX = cell(NumSims, 6);



% row1 is q1 + q2 <= 1
A  = [1 1 0 0 0];
b = [1];
% define equality constraints
Aeq = [];
beq = [];
% bounds
numx = 5;
lb = zeros(1, numx); % hard code 5D
ub = lb + Inf;

% load normlization_minmax.  this is used to rescale x, the variables we
% optimize

% fill in XX
for i = 1:NumSims

    % get parameters, save, load
    paramsi = paramgrid(i,:);
    XX{i,1} = paramsi;
    BPi = paramsi(1:9);
    ICi = paramsi(10:14);

% find initial conditions R0, E0, L0, T0, R0. use steady state
ysteady = findSteadyState(BPi, ICi, ReportedImmigration(1));
% rescale to correct total population.  localIC(1) has correct TP0
ysteady = ysteady*ICi(1)/sum(ysteady);
localIC = ysteady;
localIC(1) = ICi(1);
    
    % setup optimizer inputs
    % f5=@(x)IncidenceError5_SS(x,BPi, localIC, ReportedImmigration, ReportedTB);
    f5=@(X)IncidenceError5_SS(X,BPi, localIC, ReportedImmigration, ReportedTB, normalization_minmax);
    x0 = [BPi(8) BPi(9) localIC(2) localIC(3) localIC(5)]; %q1 q2 E0 L0 R0

    X0 = normalizex(x0, normalization_minmax );
    % optimize, store results
    % [xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f5, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run
    [Xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f5, X0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run

    xmin5 = unnormalizeX(Xmin5, normalization_minmax);

    XX{i,2} = xmin5;
    XX{i,3} = {exitflag, output, lambda, grad, hessian};
    
    % compute other XX output; pop vs time and incidence error

    % update parameters to optimal ones
    BPi(8) = xmin5(1); %q1
    BPi(9) = xmin5(2); %q2
    ICi(2) = xmin5(3); %E0
    ICi(3) = xmin5(4); %L0
    ICi(5) = xmin5(5); %R0

    [XELTRi, TBIncidencei] = solveGuoWu(BPi, ICi, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = TBIncidencei;
    XX{i,6} = norm(TBIncidencei-ReportedTB);

end


save("XX.mat", "XX")
save('paramcell.mat')
%% Functions



%% IncidenceError, function to optmmize put into fmincon
function err = IncidenceError5_SS(X, bioParameters, initialConditions, ImmigrationRate, ReportedTB, normalization_minmax)
% input: x is in R5, q1, q2, E0, L0, R0

x = unnormalizeX(X,normalization_minmax);

% load bioParameters
localBP = bioParameters;
localBP(8)=x(1); % q1
localBP(9)=x(2); % q2

%load initialConditions
localIC = initialConditions;
localIC(2)=x(3); % E0
localIC(3)=x(4); % L0
localIC(5)=x(5);

% ysteady = findSteadyState(localBP, localIC, ImmigrationRate(1));
% 
% % rescale to correct total population.  localIC(1) has correct TP0
% ysteady = ysteady*localIC(1)/sum(ysteady);
% 
% localIC = ysteady;
% localIC(1) = initialConditions(1);

[~, EstimatedTB] = solveGuoWu(localBP, localIC, ImmigrationRate);

   
err = norm((EstimatedTB-ReportedTB));
%err = norm((EstimatedTB-ReportedTB)./ReportedTB);
end

%% normalize and unnormalize

function X = normalizex(x, normalizing_minmax)
%% input x in R5, x=q1 q2 E0 L0 R0
% normalizing_minmax in R^2x5.  first row is min, last row is max. for each
% u in x, normalized U = (u-u_min)/(u_max-u_min)

X = (x - normalizing_minmax(1,:))./(normalizing_minmax(2,:)-normalizing_minmax(1,:));

% sometimes you get a small negative value for X.  replace with 0.

% X = X.*(X>=0);
end

function x = unnormalizeX(X, normalizing_minmax)
%% input x in R5, x=q1 q2 E0 L0 R0
% normalizing_minmax in R^2x5.  first row is min, last row is max. for each
% u in x, normalized U = (u-u_min)/(u_max-u_min)



x = (X ).*(normalizing_minmax(2,:)-normalizing_minmax(1,:))+ normalizing_minmax(1,:);

end

