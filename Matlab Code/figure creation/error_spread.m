numSims = size(XX,1);

errors = zeros(numSims,1);
for k=1:numSims
    errors(k) = XX{k,6};
end

figure 
plot(errors)
title('Error vs experiment')

display(['The standard deviation of error is ',num2str(std(errors))])
display(['The min and max error are respectively ', num2str(min(errors)), ' ', num2str(max(errors)) ])

% , ' ' , num2str(max(errors))



%% analyze errors

[error_min,min_index] = min(errors);
[error_max,max_index] = max(errors);

% - 3rd column is output of fmincon, {exitflag, output, lambda, grad, hessian};"

fminconOutput_min =  XX{min_index, 3};

exitflag_min = fminconOutput_min{1};
% xoptimal_min = fminconOutput_min{2};
gradient_min = fminconOutput_min{4};
hessian_min = fminconOutput_min{5};

%
fminconOutput_max =  XX{max_index, 3};
exitflag_max = fminconOutput_max{1};
xoptimal_max = fminconOutput_max{2};
gradient_max = fminconOutput_max{4};
hessian_max = fminconOutput_max{5};

%% 

display('Error minimizing experiment: gradient, hessian, hessian eigenvalues, and exit info.')

display(gradient_min');
display(hessian_min);
display(eig(hessian_min));
display( fminconOutput_min{2});

%% 
display('Error maximizing experiment: gradient, hessian, and exit info.')


display(gradient_max');
display(hessian_max);
display(eig(hessian_max));
display( fminconOutput_max{2});



%% Compare optimized E0 L0 R0 to steady state


% load the starting parameters and optimized x
startingparams_min = XX{min_index,1};
xoptimal_min = XX{min_index,2};


q1_optimal = xoptimal_min(1);
q2_optimal = xoptimal_min(2);
E0_optimal = xoptimal_min(3);
L0_optimal = xoptimal_min(4);
R0_optimal = xoptimal_min(5);

% extract preoptimized parameters
bioParameters = startingparams_min(1:9);
TP0 = startingparams_min(10);
T0 = startingparams_min(13);

% update to optimal q1 and q2
bioParameters(8) = q1_optimal;
bioParameters(9) = q2_optimal;

% find steady state
IC = [TP0, E0_optimal, L0_optimal, T0, R0_optimal];

pi0 = ReportedImmigration(1);
% pi0 = 270635;

ICX = IC;
ICX(1) = IC(1)-sum(IC(2:end));

[y,k] = findSteadyState(bioParameters,IC,pi0);

y2 = y*sum(ICX)/sum(y);

display('Relative error of steady state compared to optimized')
display((y2-ICX)./ICX)