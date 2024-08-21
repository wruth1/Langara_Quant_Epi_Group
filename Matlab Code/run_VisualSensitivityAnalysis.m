
%% SENSITIVITY ANALYSIS


addpath('functions\')
addpath('figure creation\')

%% Parameters Setup.  Hard code parameter values.
% load('Parameters20062020.mat')



%% create lists for sensitivity analysis
load('pi_extrapolated.mat')
load('qvec.mat')
load('GOODPARAM.mat')

bioParameters_start = GOODPARAM(1:8);
bioParameters_start(8) = GOODOUTPUT(4);
IC_start = GOODOUTPUT(5:9);

Params_start = [bioParameters_start];
scale_dP = 10/100;

num_params = length(Params_start);

paramcell = cell(num_params,1);

beta_list = [1, 0.5];
for k=1:num_params
    paramcell{k} = [(1-scale_dP)*Params_start(k) Params_start(k) (1+scale_dP)*Params_start(k)];
end

% paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, u_list, q1_list, q2_list, q3_list, X0_list,E0_list, L0_list, T0_list,R0_list, power_relapse_list};
%%
load('piW_base.mat')
qvec2 = qvecW_base;

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
     12 X0, 13 E0, 14 L0, 15 T0, 16 R0 % initial conditions
     17 power_relapse, the power to raise error_incidence proportion by
- 2nd column is vector in R9, optimal [qE qL qR sigma, X0 E0 L0 T0 R0]
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
    
    BPi = paramsi(1:8);
    ICi = IC_start;

        [XELTRi, EstimatedIncidence, EstimatedPrevalence] = solveGuoWu4_extrapolate(BPi, ICi, pi_extrapolated,qvec2);


    % [XELTRi, EstimatedIncidence, EstimatedPrevalence] = solveGuoWu4(BPi, newIC, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = EstimatedIncidence;
    XX{i,6} = EstimatedPrevalence;
    % XX{i,7} = fval2;
    % XX{i,8} = localIC; % the initial conditions and sigma_new used to compute q_ELR
    % XX{i,9} = u_new;

% ysteady = findSteadyState2(localBP, initialConditions, ImmigrationRate(1));
% 


end



%% save data
save("./data and results/XX_extra.mat", "XX")
save("./data and results/paramcell_extra.mat", "paramcell")
%% Functions






