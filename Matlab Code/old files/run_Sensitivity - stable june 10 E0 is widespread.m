
%% SENSITIVITY ANALYSIS

addpath('functions\')
addpath('figure creation\')

%% Parameters Setup.  Hard code parameter values.
load('Parameters20062020.mat')



%% create lists for sensitivity analysis

% scale_param = 9.5/10; % assume <1
% scale_list = scale_param.^[1,0,-1];


% beta_list = [beta];
% beta_list = beta_list * scale_list;
beta_list = [0.00000001 0.000007];




% w_list = w;
% w_list = w_list*scale_list;


w=0.5;
w_list = [0.5];
% p_list = [0.00849, 0.011] / w;
% p_list = [p_list, w*0.05

p_list = [4/100, 5/100, 10/100]/w; % p*w range is 2-15%

% w_list = w;
% v_list = v;
% v_list = v_list*scale_list;            % The rate of slow progression to active TB due to reactivation
% v_list = [0.00161, 0.00224, 0.00487];
v_list = [1e-4, 2e-4 5e-4 1.7e-3]; % Menzies mediansx2, Guo-Wu, Dowdy2013


a_list = [0.0522 0.06];              % TB-caused death rate
% a_list = a_list*scale_list;
% d_list = [d];               % Constant rate of recovery by nature or treatment
% d_list = d_list*scale_list;
d_list = [ 0.6, 0.7, 0.8];

% relapse rate
u_list = [3e-3];

% % uniform grid
% u_dowdy = 3e-3;
% u_aparico = 8e-5;
% numu = 6;
% stepu = (u_dowdy-u_aparico)/(numu-1);
% % u_list = [8e-5 3e-3]; %dowdy 
% u_list = u_aparico:stepu:u_dowdy;

n_list = [0.0071];  % natural removal rate, mode from Statcan
% n_list = n_list *scale_list;

q1_list = [q1];
q2_list = [q2];
q3_list = [0.1];


X0 = 2.435289381318261e+06;
E0 = 1.484646245171584e+04;

L0 = 3.096493365566524e+06;
T0 = 1.070931823865891e+03;
R0 =      6.392498588396340e+05;
TFP0 = X0 + E0+L0+T0+R0;

X0_list = [X0];
E0_list = [E0];
L0_list = [L0];
T0_list = [T0];
R0_list = [R0];

paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, u_list, q1_list, q2_list, q3_list, X0_list,E0_list, L0_list, T0_list,R0_list};


%% Big Data Structure XX
%{
XX stores all the data.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st column is vector |BP|+|IC|, starting parameters
    Parameter order:
     1 beta = bioParameters(1); %TB infectivity
     2 p  = bioParameters(2); % ~probability someone in E goes straight into E; pi in Guo-Wu
     3 w = bioParameters(3); % period of time new infectee considered E rather than L
     4 v = bioParameters(4); % rate people in L develop TB
     5 a = bioParameters(5); % TB death rate for people in T
     6 d = bioParameters(6); % non-TB death rate
     7 n = bioParameters(7); % recovery rate; delta in Guo-Wu
     8 u = bioParameters(8); % percentage people in R develop active TB 
     9 qE = q1 = bioParameters(9); % percent of immigrants into E
     10 qL = q2 = bioParameters(10); % percent of immigrants into L
     11 qR = q3 = bioParameters(11); % percent of immigrants into R
- 2nd column is vector in R6, optimal [q1 q2 X0 E0 L0 R0]
- 3rd column is output of fmincon, {exitflag, output, lambda, grad, hessian};
- 4th column is matrix, population XELTR vs time (using optimized params)
- 5th column is TB incidence (using optimized params)
- 6th column is TB prevalence
- 7th column is error of TB incidence using optimized params
%}

% use Will's expand_grid.m to create every combination of parameters
paramgrid = expand_grid(paramcell{:}); 

% initialize XX
NumSims = size(paramgrid,1);
XX = cell(NumSims, 9);


% optimize

% fill in XX
for i = 1:NumSims

    % get parameters, save, load
    paramsi = paramgrid(i,:);
    XX{i,1} = paramsi;
    BPi = paramsi(1:11);
    ICi = paramsi(12:end);

    TFP0 = sum(ICi);


    % update initial conditions

    % % find initial conditions R0, E0, L0, T0, R0. use steady state
    ysteady = findSteadyState2(BPi, ICi, ReportedImmigration(1));
    % rescale to correct total population.  
    ysteady = ysteady*TFP0/sum(ysteady);
    localIC = ysteady;

    

    % set up fmincon input

    % x = [q1 q2 q3 s]

    numx=9;
    A = zeros(1,numx);
    A(1,1:3) = [1 1 1]; % q1+q2+q3 <= 100%
    b(1) = 1 ;

    Aeq = [0 0 0 0 1 1 1 1 1]; % total foreign-born population
    
    beq = [TFP0];

    %T0 = 1081
    Aeq(2,1:numx) = [0 0 0 0 0 0 0 1 0];
    beq(2) = [1081];
 
    lb = zeros(1, numx); 
    ub = lb + Inf;

    % setup optimizer inputs
    

    f2=@(x)IncidenceError4(x,BPi, ReportedImmigration, ReportedIncidence);
    x0 = [BPi(9) BPi(10) BPi(11) BPi(8)]; %q1 q2 q3
    x0 = [x0, localIC]; 

    % x0 = [BPi(9) BPi(10) BPi(11) ]; %q1 q2 q3 sigma
    [xmin2,fval2,exitflag,output,lambda,grad,hessian] = fmincon(f2, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run



    XX{i,2} = xmin2;
    XX{i,3} = {exitflag, output, lambda, grad, hessian};
    
    % compute other XX output; pop vs time and incidence error

    % update parameters to optimal ones
    BPi(9) = xmin2(1); %q1
    BPi(10) = xmin2(2); %q2
    BPi(11) = xmin2(3); %q3
    BPi(8) = xmin2(4); % s

    newIC = xmin2(5:9);

    % find initial conditions R0, E0, L0, T0, R0. use steady state

    [XELTRi, EstimatedIncidence, EstimatedPrevalence] = solveGuoWu4(BPi, newIC, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = EstimatedIncidence;
    XX{i,6} = EstimatedPrevalence;
    XX{i,7} = fval2;
    XX{i,8} = localIC; % the initial conditions and sigma_new used to compute q_ELR
    % XX{i,9} = u_new;

% ysteady = findSteadyState2(localBP, initialConditions, ImmigrationRate(1));
% 


end

%% save data

path_name = '';
version_name = '28';
save("./data and results/XX.mat", "XX")
save("./data and results/paramcell.mat", "paramcell")
%% Functions






