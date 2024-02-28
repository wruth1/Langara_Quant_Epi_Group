
%% SENSITIVITY ANALYSIS


%% Parameters Setup.  Hard code parameter values.

ReportedImmigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309]; %immigration into canada 2010-2020
ReportedTB = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020



% bioParameters

%Parameters from York Paper
beta = 1*10^-8;        % Transmission rate within foreign-born population in Canada.
p = 0.05;              % P(progresses directly to active TB stage from early latent without treatment)
w = 0.4;               % (q1+q2) Pro portion of all LTBI immigrants pass through latent stage in the first 2.5 years
% v = 0.0002;            % The rate of slow progression to active TB due to reactivation
v = p*w/15;    %15x more likely to develop
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
R0 = 1e6; 
% TP = total FB population

% initial % of populations from Guo Wu
E0 = 0.001735006* TP0;
L0 = 0.21218547 * TP0; 


bioParameters = [beta p w v a d n q1 q2];
initialConditions = [TP0, E0, L0, T0, R0];

ysteady = findSteadyState(bioParameters, initialConditions, ReportedImmigration(1));


display(['The steady state total population is: ', num2str(sum(ysteady))]);

estimated_incidence = getTBIncidenceRate(ysteady, p, w, v);

display(['The steady state incidence is: ', num2str(estimated_incidence)]);
%%

% initial conditions are steady state; fudge X0 a bit
initialConditions = ysteady;
initialConditions(1) = TP0;

%%



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


% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)

BPi = bioParameters;
ICi = initialConditions;


    % setup optimizer inputs
    f5=@(x)IncidenceError5(x,BPi, ICi, ReportedImmigration, ReportedTB);
    x0 = [BPi(8) BPi(9) ICi(2) ICi(3) ICi(5)]; %q1 q2 E0 L0


problem = createOptimProblem('fmincon','objective',...
    f5,'x0',ICi,'Aineq',A,'bineq',b ,'Aeq',Aeq, 'beq', beq, 'lb',lb,'ub',ub);

ms = MultiStart;

[x,f] = run(ms,problem,20);

% fill in XX
% for i = 1:NumSims
% 
%     % get parameters, save, load
%     paramsi = paramgrid(i,:);
%     XX{i,1} = paramsi;
%     BPi = paramsi(1:9);
%     ICi = paramsi(10:14);
%     
%     % setup optimizer inputs
%     f5=@(x)IncidenceError5(x,BPi, ICi, ReportedImmigration, ReportedTB);
%     x0 = [BPi(8) BPi(9) ICi(2) ICi(3) ICi(5)]; %q1 q2 E0 L0
% 
%     % optimize, store results
%     [xmin5,fval5,exitflag,output,lambda,grad,hessian] = fmincon(f5, x0 , A , b, Aeq, beq, lb, ub) ; % about 10 seconds to run
%     XX{i,2} = xmin5;
%     XX{i,3} = {exitflag, output, lambda, grad, hessian};
%     
%     % compute other XX output; pop vs time and incidence error
% 
%     % update parameters to optimal ones
%     BPi(8) = xmin5(1); %q1
%     BPi(9) = xmin5(2); %q2
%     ICi(2) = xmin5(3); %E0
%     ICi(3) = xmin5(4); %L0
%     ICi(5) = xmin5(5); %R0
% 
%     [XELTRi, TBIncidencei] = solveGuoWu(BPi, ICi, ReportedImmigration);
% 
%     XX{i,4} = XELTRi;
%     XX{i,5} = TBIncidencei;
%     XX{i,6} = norm(TBIncidencei-ReportedTB);
% 
% end
% 
% 
% save('XX.mat')
% save('paramcell.mat')
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


