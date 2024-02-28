%% Analyze Sensitivity

%{
A few tools to analyze the results following sensitivity analysis.

1. Histogram of optimal values
2. Error vs iteration.
3. 
%}
%%
imageName = 'Run16_x6-v2';

NumSims=size(XX,1);
numSims = NumSims;


allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);
allX0 = zeros(NumSims,1);
allE0 = zeros(NumSims,1);
allL0 = zeros(NumSims,1);
allR0 = zeros(NumSims,1);

for i=1:NumSims
    currentx = XX{i,2};
    allq1(i)=currentx(1);
    allq2(i) = currentx(2);
    allX0(i) = currentx(3);
    allE0(i) = currentx(4);
    allL0(i) = currentx(5);
    allR0(i) = currentx(6);
end

numbins=10;
%% Plot Histograms
figure(1)

subplot(2,3,1)

% histfit(allq1,numbins)
histfit(allq1,numbins)
title('Distribution of optimal q1')

subplot(2,3,2)
histfit(allq2,numbins)
title('Distribution of optimal q2')


subplot(2,3,3)
histfit(allE0,numbins)
title('Distribution of optimal X0')

subplot(2,3,4)
histfit(allE0,numbins)
title('Distribution of optimal E0')

subplot(2,3,5)
histfit(allL0,numbins)
title('Distribution of optimal L0')

subplot(2,3,6)
histfit(allR0,numbins)
% histfit(allR0,numbins)
title('Distribution of optimal R0')
saveas(gcf,['Optimized parameter distribution across sensitivity analysis',imageName, '.png'])

%% Error spread

% numSims = size(XX,1);

errors = zeros(numSims,1);
for k=1:numSims
    errors(k) = XX{k,6};
end

figure(2)
plot(errors)
title('Error vs experiment')

display(['The standard deviation of error is ',num2str(std(errors))])
display(['The min and max error are respectively ', num2str(min(errors)), ' ', num2str(max(errors)) ])

%% Exit information for error-minimizing experiment

% find error minimiizing and maximizing experiments
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


display('Error minimizing experiment: gradient, hessian, hessian eigenvalues, and exit info.')

display(gradient_min');
display(hessian_min);
display(eig(hessian_min));
display( fminconOutput_min{2});


