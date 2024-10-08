
%% SENSITIVITY ANALYSIS

addpath('functions\')
addpath('figure creation\')

%% Parameters Setup.  Hard code parameter values.
load('Parameters20102020.mat')



%% create lists for sensitivity analysis

beta_list = [beta];
p_list = [p];

scale_param = 9/10; % assume <1
scale_list = scale_param.^[1,0,-1,-2];

w_list = w;
w_list = w_list*scale_list;

% w_list = w;
v_list = v;
v_list = v_list*scale_list;            % The rate of slow progression to active TB due to reactivation

a_list = [a];              % TB-caused death rate
d_list = [d];               % Constant rate of recovery by nature or treatment
n_list = [n]; 

q1_list = [q1];
q1_list = q1_list*scale_list;
q2_list = [q2];
q2_list = q2_list*scale_list;

X0_list = [X0];
E0_list = [E0];
L0_list = [L0];
T0_list = [T0];
R0_list = [R0];

paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, q1_list, q2_list, X0_list,E0_list, L0_list, T0_list,R0_list};

% save('Experiment_Parameters' )

% load('normalization_minmax6.mat')
% normalization_minmax = normalize_minmax;

%%
%% Big Data Structure XX
%{
XX stores all the data.  Each row corresponds to a combination of
parameters / initial conditions.

- 1st column is vector |BP|+|IC|, starting parameters
- 2nd column is vector in R6, optimal [q1 q2 X0 E0 L0 R0]
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






% load normlization_minmax.  this is used to rescale x, the variables we
% optimize
load 'normalize_minmax-x7-run18.mat';
normalization_minmax = normalize_minmax;

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
    

    % set up fmincon input

    % parse normalization factors
    Q1min = normalization_minmax(1,1);
    Q2min = normalization_minmax(1,2);

    Q1max = normalization_minmax(2,1);
    Q2max = normalization_minmax(2,2);

    % initial distribution
    X0min = normalization_minmax(1,3);
    E0min = normalization_minmax(1,4);
    L0min = normalization_minmax(1,5);
    T0min = normalization_minmax(1,6);
    R0min = normalization_minmax(1,7);

    X0max = normalization_minmax(2,3);
    E0max = normalization_minmax(2,4);
    L0max = normalization_minmax(2,5);
    T0max = normalization_minmax(2,6);
    R0max = normalization_minmax(2,7);
        

    numx = 7;
    % row1 is Q1 + Q2 <= 1
    A(1,1:numx)  = [ Q1max-Q1min, Q2max-Q2min, 0 0 0 0 0];
    b(1) = [1-Q1min-Q2min];


    % A(2,1:5) = [0, 0, E0max-E0min, L0max-L0min, R0max-R0min];
    % b(2) = [ICi(1) - ysteady(1) - ysteady(4) - E0min - L0min - R0min];% TP_0 - X0 - T0 - Emin - Lmin - Rmin

    % define equality constraints
    % X0 + E0 + L0 + R0 = TP0 - T0
    Aeq(1,1:numx) = [0, 0, X0max-X0min, E0max-E0min, L0max-L0min, T0max-T0min, R0max-R0min];
    TP0 = sum(ysteady);
    beq(1) = [TP0 - X0min - E0min - L0min - T0min - R0min];% TP_0 - T0 - Emin - Lmin - Rmin
    
    % pwE + vL = 14.1*(total population)

    % [beta p w v a d n q1 q2]
    % p_i = BPi(2);
    % w_i = BPi(3);
    % v_i = BPi(4);
    % 
    % ri0 = ReportedIncidence(1)/100000;
    % Aeq(2,1:numx) = [0, 0, ... 
    %     (X0max-X0min)*(-ri0), ...
    %     (E0max-E0min)*(p_i*w_i-ri0), ...
    %     (L0max-L0min)*(v_i-ri0), ...
    %     (T0max-T0min)*(-ri0), ...
    %     (R0max-R0min)*(-ri0)];
    % beq(2) = ri0*(X0min + E0min + L0min + T0min + R0min) - p_i*w_i*E0min - v_i*L0min;

    
    % bounds
    numx = 7;
    lb = zeros(1, numx); % hard code 5D
    ub = lb + Inf;

    % setup optimizer inputs
    f7=@(X)IncidenceError7_SS(X,BPi, localIC, ReportedImmigration, ReportedTB, normalization_minmax);
    x0 = [BPi(8) BPi(9) localIC(1) localIC(2) localIC(3) localIC(4) localIC(5)]; %q1 q2 X0 E0 L0 R0
    X0 = normalizex(x0, normalization_minmax );
    [Xmin7,fval7,exitflag,output,lambda,grad,hessian] = fmincon(f7, X0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run

    % X0 = normalizex(x0, normalization_minmax );
    % optimize, store results
    % [Xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f6, X0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run

    xmin7 = unnormalizeX(Xmin7, normalization_minmax);

    XX{i,2} = xmin7;
    XX{i,3} = {exitflag, output, lambda, grad, hessian};
    
    % compute other XX output; pop vs time and incidence error

    % update parameters to optimal ones
    BPi(8) = xmin7(1); %q1
    BPi(9) = xmin7(2); %q2
    ICi(1) = xmin7(3); %X0
    ICi(2) = xmin7(4); %E0
    ICi(3) = xmin7(5); %L0
    ICi(4) = xmin7(6);
    ICi(5) = xmin7(7); %R0

    [XELTRi, TBIncidencei] = solveGuoWu3(BPi, ICi, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = TBIncidencei;
    XX{i,6} = norm(TBIncidencei-ReportedTB);

end


save("XX.mat", "XX")
save('paramcell.mat')
%% Functions



%% IncidenceError, function to optmmize put into fmincon
function err = IncidenceError7_SS(X, bioParameters, initialConditions, ImmigrationRate, ReportedTB,normalization_minmax)
% input: x is in R6, q1, q2, X0, E0, L0, R0

x = unnormalizeX(X,normalization_minmax);

% load bioParameters
localBP = bioParameters;
localBP(8)=x(1); % q1
localBP(9)=x(2); % q2

%load initialConditions
localIC = initialConditions;
localIC(1)=x(3); %X0
localIC(2)=x(4); % E0
localIC(3)=x(5); % L0
localIC(4)=x(6); % T0
localIC(5)=x(7); %R0

% ysteady = findSteadyState(localBP, localIC, ImmigrationRate(1));
% 
% % rescale to correct total population.  localIC(1) has correct TP0
% ysteady = ysteady*localIC(1)/sum(ysteady);
% 
% localIC = ysteady;
% localIC(1) = initialConditions(1);

[~, EstimatedTB] = solveGuoWu3(localBP, localIC, ImmigrationRate);

   
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

