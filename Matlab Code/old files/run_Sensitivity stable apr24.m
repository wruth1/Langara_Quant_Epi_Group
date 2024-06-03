
%% SENSITIVITY ANALYSIS

addpath('functions\')
addpath('figure creation\')

%% Parameters Setup.  Hard code parameter values.
load('Parameters20102020.mat')
load('ReportedTB201020.mat')


%% create lists for sensitivity analysis

scale_param = 9.5/10; % assume <1
scale_list = scale_param.^[1,0,-1];


beta_list = [beta];
% beta_list = beta_list * scale_list;

p_list = [p];
p_list = p_list*scale_list;

w_list = w;
w_list = w_list*scale_list;

% w_list = w;
v_list = v;
v_list = v_list*scale_list;            % The rate of slow progression to active TB due to reactivation

a_list = [a];              % TB-caused death rate
a_list = a_list*scale_list;
d_list = [d];               % Constant rate of recovery by nature or treatment
d_list = d_list*scale_list;

n_list = [n];  % natural removal rate
% n_list = n_list *scale_list;

q1_list = [q1];
q2_list = [q2];

X0_list = [X0];
E0_list = [E0];
L0_list = [L0];
T0_list = [T0];
R0_list = [R0];

paramcell = {beta_list, p_list, w_list,  v_list, a_list, d_list, n_list, q1_list, q2_list, X0_list,E0_list, L0_list, T0_list,R0_list};

%% Steady state from Guo-Wu

PI = ReportedImmigration(1);

a1 = -beta*(n+w)*(n+v)*(n+a+d);
a2 = n*PI*(q1*w*(p*n+v)+q2*v*(n+w));
a3 = beta*(1-q2)*PI*(p*n+v)+beta*q2*PI*v*(n+w)-n*( (n+w)*(n+v)*(n+a+d) );

Tstar = (a2+sqrt(a2^2-4*a1*a3))/(-2*a1);

%% 
Tstar= (-(PI * beta * n * p * q2 * w) + (PI * beta * n * p * w) + (PI * beta * n * q2 * v) + (PI * beta * v * w) - (a * n ^ 3) - (a * n ^ 2 * v) - (a * n ^ 2 * w) - (a * n * v * w) - (d * n ^ 3) - (d * n ^ 2 * v) - (d * n ^ 2 * w) - (d * n * v * w) - (n ^ 4) - (n ^ 3 * v) - (n ^ 3 * w) - (n ^ 2 * v * w) + sqrt((PI ^ 2 * beta ^ 2 * n ^ 2 * p ^ 2 * q2 ^ 2 * w ^ 2 - 2 * PI ^ 2 * beta ^ 2 * n ^ 2 * p ^ 2 * q2 * w ^ 2 - 2 * PI ^ 2 * beta ^ 2 * n ^ 2 * p * q2 ^ 2 * v * w + PI ^ 2 * beta ^ 2 * n ^ 2 * p ^ 2 * w ^ 2 + 2 * PI ^ 2 * beta ^ 2 * n ^ 2 * p * q2 * v * w + PI ^ 2 * beta ^ 2 * n ^ 2 * q2 ^ 2 * v ^ 2 - 2 * PI ^ 2 * beta ^ 2 * n * p * q2 * v * w ^ 2 + 4 * PI * a * beta * n ^ 4 * p * q1 * w + 2 * PI * a * beta * n ^ 4 * p * q2 * w + 4 * PI * a * beta * n ^ 3 * p * q1 * v * w + 4 * PI * a * beta * n ^ 3 * p * q1 * w ^ 2 + 2 * PI * a * beta * n ^ 3 * p * q2 * v * w + 2 * PI * a * beta * n ^ 3 * p * q2 * w ^ 2 + 4 * PI * a * beta * n ^ 2 * p * q1 * v * w ^ 2 + 2 * PI * a * beta * n ^ 2 * p * q2 * v * w ^ 2 + 4 * PI * beta * d * n ^ 4 * p * q1 * w + 2 * PI * beta * d * n ^ 4 * p * q2 * w + 4 * PI * beta * d * n ^ 3 * p * q1 * v * w + 4 * PI * beta * d * n ^ 3 * p * q1 * w ^ 2 + 2 * PI * beta * d * n ^ 3 * p * q2 * v * w + 2 * PI * beta * d * n ^ 3 * p * q2 * w ^ 2 + 4 * PI * beta * d * n ^ 2 * p * q1 * v * w ^ 2 + 2 * PI * beta * d * n ^ 2 * p * q2 * v * w ^ 2 + 4 * PI * beta * n ^ 5 * p * q1 * w + 2 * PI * beta * n ^ 5 * p * q2 * w + 4 * PI * beta * n ^ 4 * p * q1 * v * w + 4 * PI * beta * n ^ 4 * p * q1 * w ^ 2 + 2 * PI * beta * n ^ 4 * p * q2 * v * w + 2 * PI * beta * n ^ 4 * p * q2 * w ^ 2 + 4 * PI * beta * n ^ 3 * p * q1 * v * w ^ 2 + 2 * PI * beta * n ^ 3 * p * q2 * v * w ^ 2 + 2 * PI ^ 2 * beta ^ 2 * n * p * v * w ^ 2 + 2 * PI ^ 2 * beta ^ 2 * n * q2 * v ^ 2 * w - 2 * PI * a * beta * n ^ 4 * p * w + 2 * PI * a * beta * n ^ 4 * q2 * v - 2 * PI * a * beta * n ^ 3 * p * v * w - 2 * PI * a * beta * n ^ 3 * p * w ^ 2 + 4 * PI * a * beta * n ^ 3 * q1 * v * w + 2 * PI * a * beta * n ^ 3 * q2 * v ^ 2 + 6 * PI * a * beta * n ^ 3 * q2 * v * w - 2 * PI * a * beta * n ^ 2 * p * v * w ^ 2 + 4 * PI * a * beta * n ^ 2 * q1 * v ^ 2 * w + 4 * PI * a * beta * n ^ 2 * q1 * v * w ^ 2 + 6 * PI * a * beta * n ^ 2 * q2 * v ^ 2 * w + 4 * PI * a * beta * n ^ 2 * q2 * v * w ^ 2 + 4 * PI * a * beta * n * q1 * v ^ 2 * w ^ 2 + 4 * PI * a * beta * n * q2 * v ^ 2 * w ^ 2 - 2 * PI * beta * d * n ^ 4 * p * w + 2 * PI * beta * d * n ^ 4 * q2 * v - 2 * PI * beta * d * n ^ 3 * p * v * w - 2 * PI * beta * d * n ^ 3 * p * w ^ 2 + 4 * PI * beta * d * n ^ 3 * q1 * v * w + 2 * PI * beta * d * n ^ 3 * q2 * v ^ 2 + 6 * PI * beta * d * n ^ 3 * q2 * v * w - 2 * PI * beta * d * n ^ 2 * p * v * w ^ 2 + 4 * PI * beta * d * n ^ 2 * q1 * v ^ 2 * w + 4 * PI * beta * d * n ^ 2 * q1 * v * w ^ 2 + 6 * PI * beta * d * n ^ 2 * q2 * v ^ 2 * w + 4 * PI * beta * d * n ^ 2 * q2 * v * w ^ 2 + 4 * PI * beta * d * n * q1 * v ^ 2 * w ^ 2 + 4 * PI * beta * d * n * q2 * v ^ 2 * w ^ 2 - 2 * PI * beta * n ^ 5 * p * w + 2 * PI * beta * n ^ 5 * q2 * v - 2 * PI * beta * n ^ 4 * p * v * w - 2 * PI * beta * n ^ 4 * p * w ^ 2 + 4 * PI * beta * n ^ 4 * q1 * v * w + 2 * PI * beta * n ^ 4 * q2 * v ^ 2 + 6 * PI * beta * n ^ 4 * q2 * v * w - 2 * PI * beta * n ^ 3 * p * v * w ^ 2 + 4 * PI * beta * n ^ 3 * q1 * v ^ 2 * w + 4 * PI * beta * n ^ 3 * q1 * v * w ^ 2 + 6 * PI * beta * n ^ 3 * q2 * v ^ 2 * w + 4 * PI * beta * n ^ 3 * q2 * v * w ^ 2 + 4 * PI * beta * n ^ 2 * q1 * v ^ 2 * w ^ 2 + 4 * PI * beta * n ^ 2 * q2 * v ^ 2 * w ^ 2 + PI ^ 2 * beta ^ 2 * v ^ 2 * w ^ 2 - 2 * PI * a * beta * n ^ 3 * v * w - 2 * PI * a * beta * n ^ 2 * v ^ 2 * w - 2 * PI * a * beta * n ^ 2 * v * w ^ 2 - 2 * PI * a * beta * n * v ^ 2 * w ^ 2 - 2 * PI * beta * d * n ^ 3 * v * w - 2 * PI * beta * d * n ^ 2 * v ^ 2 * w - 2 * PI * beta * d * n ^ 2 * v * w ^ 2 - 2 * PI * beta * d * n * v ^ 2 * w ^ 2 - 2 * PI * beta * n ^ 4 * v * w - 2 * PI * beta * n ^ 3 * v ^ 2 * w - 2 * PI * beta * n ^ 3 * v * w ^ 2 - 2 * PI * beta * n ^ 2 * v ^ 2 * w ^ 2 + a ^ 2 * n ^ 6 + 2 * a ^ 2 * n ^ 5 * v + 2 * a ^ 2 * n ^ 5 * w + a ^ 2 * n ^ 4 * v ^ 2 + 4 * a ^ 2 * n ^ 4 * v * w + a ^ 2 * n ^ 4 * w ^ 2 + 2 * a ^ 2 * n ^ 3 * v ^ 2 * w + 2 * a ^ 2 * n ^ 3 * v * w ^ 2 + a ^ 2 * n ^ 2 * v ^ 2 * w ^ 2 + 2 * a * d * n ^ 6 + 4 * a * d * n ^ 5 * v + 4 * a * d * n ^ 5 * w + 2 * a * d * n ^ 4 * v ^ 2 + 8 * a * d * n ^ 4 * v * w + 2 * a * d * n ^ 4 * w ^ 2 + 4 * a * d * n ^ 3 * v ^ 2 * w + 4 * a * d * n ^ 3 * v * w ^ 2 + 2 * a * d * n ^ 2 * v ^ 2 * w ^ 2 + 2 * a * n ^ 7 + 4 * a * n ^ 6 * v + 4 * a * n ^ 6 * w + 2 * a * n ^ 5 * v ^ 2 + 8 * a * n ^ 5 * v * w + 2 * a * n ^ 5 * w ^ 2 + 4 * a * n ^ 4 * v ^ 2 * w + 4 * a * n ^ 4 * v * w ^ 2 + 2 * a * n ^ 3 * v ^ 2 * w ^ 2 + d ^ 2 * n ^ 6 + 2 * d ^ 2 * n ^ 5 * v + 2 * d ^ 2 * n ^ 5 * w + d ^ 2 * n ^ 4 * v ^ 2 + 4 * d ^ 2 * n ^ 4 * v * w + d ^ 2 * n ^ 4 * w ^ 2 + 2 * d ^ 2 * n ^ 3 * v ^ 2 * w + 2 * d ^ 2 * n ^ 3 * v * w ^ 2 + d ^ 2 * n ^ 2 * v ^ 2 * w ^ 2 + 2 * d * n ^ 7 + 4 * d * n ^ 6 * v + 4 * d * n ^ 6 * w + 2 * d * n ^ 5 * v ^ 2 + 8 * d * n ^ 5 * v * w + 2 * d * n ^ 5 * w ^ 2 + 4 * d * n ^ 4 * v ^ 2 * w + 4 * d * n ^ 4 * v * w ^ 2 + 2 * d * n ^ 3 * v ^ 2 * w ^ 2 + n ^ 8 + 2 * n ^ 7 * v + 2 * n ^ 7 * w + n ^ 6 * v ^ 2 + 4 * n ^ 6 * v * w + n ^ 6 * w ^ 2 + 2 * n ^ 5 * v ^ 2 * w + 2 * n ^ 5 * v * w ^ 2 + n ^ 4 * v ^ 2 * w ^ 2))) / beta / (a * n ^ 2 + a * n * v + a * n * w + a * v * w + d * n ^ 2 + d * n * v + d * n * w + d * v * w + n ^ 3 + n ^ 2 * v + n ^ 2 * w + n * v * w) / 0.2e1;

% check 

% ysteady = findSteadyState(BPi, ICi, ReportedImmigration(1));


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
    

    % set up fmincon input

    numx=2;
    A(1,1:numx) = [1 1];
    b(1) = 1;

 
    numx = 2;
    lb = zeros(1, numx); 
    ub = lb + Inf;

    % setup optimizer inputs
    % f5=@(x)IncidenceError5_SS(x,BPi, localIC, ReportedImmigration, ReportedTB);
    f2=@(x)IncidenceError2_SS(x,BPi, localIC, ReportedImmigration, ReportedIncidence,ReportedPrevalence);
    x0 = [BPi(8) BPi(9)]; %q1 q2 X0 E0 L0 R0
    [xmin2,fval2,exitflag,output,lambda,grad,hessian] = fmincon(f2, x0 , A , b, [], [], lb, ub) ; % about 10 seconds to run

    % X0 = normalizex(x0, normalization_minmax );
    % optimize, store results
    % [Xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f6, X0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run

    % xmin5 = unnormalizeX(Xmin5, normalization_minmax);

    XX{i,2} = xmin2;
    XX{i,3} = {exitflag, output, lambda, grad, hessian};
    
    % compute other XX output; pop vs time and incidence error

    % update parameters to optimal ones
    BPi(8) = xmin2(1); %q1
    BPi(9) = xmin2(2); %q2
    % ICi(1) = xmin7(3); %X0
    % ICi(2) = xmin7(4); %E0
    % ICi(3) = xmin7(5); %L0
    % ICi(4) = xmin7(6); %T0
    % ICi(5) = xmin7(7); %R0
    ICi = ysteady;

    [XELTRi, TBIncidencei] = solveGuoWu3(BPi, ICi, ReportedImmigration);

    XX{i,4} = XELTRi;
    XX{i,5} = TBIncidencei;
    XX{i,6} = fval2;

end


save("XX.mat", "XX")
save('paramcell.mat')
%% Functions



%% IncidenceError, function to optmmize put into fmincon
function err = IncidenceError2_SS(x, bioParameters, initialConditions, ImmigrationRate, ReportedIncidence,ReportedPrevalence)
% input: x is in R2, q1, q2

localBP = bioParameters;
localBP(8)=x(1); % q1
localBP(9)=x(2); % q2

%load initialConditions
localIC = initialConditions;


[XELTR, EstimatedIncidence] = solveGuoWu3(localBP, localIC, ImmigrationRate);

   
err1 = norm((EstimatedIncidence-ReportedIncidence)./ReportedIncidence);

EstimatedPrevalence = XELTR(1:end-1,4)';

err2 = norm((EstimatedPrevalence-ReportedPrevalence)./ReportedPrevalence);

err = err1 + err2;

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

