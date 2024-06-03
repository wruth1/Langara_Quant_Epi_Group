%% Analyze Sensitivity

%{
A few tools to analyze the results following sensitivity analysis.

1. Histogram of optimal values
2. Error vs iteration.
3. 
%}

set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultContourLineWidth',1.5)
set(groot,'DefaultFunctionContourLineWidth',2)
set(groot,'defaultAxesFontSize',22)
set(groot,'defaultLegendFontSize',22)


addpath('data and results\')

%%

load('XX28.mat')
load('paramcell28.mat')


imageName = 'Run28';

NumSims=size(XX,1);
numSims = NumSims;





%% Houben comparison
% Compare q1 and q2, our estimated prevalence among new immigrants, against Houben's global prevalence of latent TB

allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);

for i=1:NumSims
    currentx = XX{i,2};
    allq1(i)=currentx(1);
    allq2(i) = currentx(2);
end

LTBI_newImmigrants = allq1+allq2;
% EoverL_newImmigrants = allq1./LTBI_newImmigrants    ;


numbins=10;
figure('units','normalized','outerposition',[0 0 1 1])


title('Comparing new immigrants to Houben estimates of global prevalence')
subplot(2,2,1)

% histfit(allq1,numbins)
histfit(LTBI_newImmigrants,numbins)
title('Prevalence among new immigrants (%)')
% average low high  21.8046423	17.28412035	28.60442942
Houben_mean = 21.80/100;
Houben_low = 17.28/100;
Houben_high = 28.60/100;




% Add vertical lines at the mean and +/- 1 standard deviation
xline(Houben_mean, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houben_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houben_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);

Houbenq1 = 0.725/100;

Houbenq1_low = 0.563/100;
Houbenq1_high = 0.920/100;


subplot(2,2,2)
histfit(allq1,numbins)
title('Infected within past 2 years (%)')
xline(Houbenq1, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);


subplot(2,2,3)
scatter(allq1,allq2,'DisplayName','q2 vs q1')
hold on

A1 = min(allq1);
% A1 = 0;
A2 =  max(allq1);
Houben_EoverL = 0.036036036;

% plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'DisplayName','Houben q1/(q1+q2)')
plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'LineStyle','--','DisplayName','Houben Mean')
legend('Location','NorthEast')
hold off


subplot(2,2,4)
scatter(allq1,allq2,'DisplayName','q2 vs q1')
hold on

A1 = 0;

plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'DisplayName','Houben Mean')
plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'LineStyle','--','DisplayName','Houben q1/(q1+q2)')

legend('Location','NorthWest')
hold off



saveas(gcf,['Histograms of q1 q2 compared to Houben ',imageName, '.png'])





%% Error spread

% numSims = size(XX,1);

errors = zeros(numSims,1);
for k=1:numSims
    errors(k) = XX{k,7};
end

figure('units','normalized','outerposition',[0 0 1 1])

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


display('Error minimizing experiment:')


display('Gradient: ')
display(gradient_min');

display('Hessian: ');
display(hessian_min);



display('Eigenvectors and eigenvalues:')
[V,D] = eig(hessian_min);
display(V)
display(diag(D))


display('Output info:')
display( fminconOutput_min{2});

%% Plot Incidence and Prevalence of error minimizing

% - 4th column is matrix, population XELTR vs time (using optimized params)
% - 5th column is TB incidence (using optimized params)

XELTR_min = XX{min_index, 4};
incidence_min = XX{min_index,5};
prevalence_min = XX{min_index,6};

load('.\data and results\ReportedTB201020.mat');
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,1,1)
plot(Years,incidence_min(1:end-1)','r--','DisplayName','Estimated');% prevalence
hold on
plot(Years,ReportedIncidence,'k','LineWidth',2, 'DisplayName','Reported')
title('Incidence vs Time')
legend('Location','East')
hold off

subplot(2,1,2)
plot(Years,prevalence_min(1:end-1),'r--','DisplayName','Estimated');% prevalence
hold on
plot(Years,ReportedPrevalence,'k','LineWidth',2, 'DisplayName','Reported')
title('Prevalence vs Time')
legend('Location','East')
hold off

saveas(gcf,['Prevalence Incidence',imageName, '.png'])

%% Plot test data - compare our model to

load('Error_aj.mat');

% compute error Paj
% total_tb = sum(XELTR_min(:,2:4));

Years2 = [Years, 2021];

estimate_Paj = zeros(length(Years2),1 ) ;

for k=1:length(Years2)
    estimate_Paj(k) = sum(XELTR_min(k,2:3))./sum(XELTR_min(k,:));
end

% estimated_Paj = total_tb./sum(XELTR_min(:,:));

figure('units','normalized','outerposition',[0 0 1 1])

plotPopulationVsTime(Years2,XELTR_min)

errorbar(years_aj, prevalence_aj, prevalence_aj-LL_aj, UL_aj-prevalence_aj,'DisplayName','Jordan Prevalence')
hold on
plot(Years2,estimate_Paj*100,'DisplayName','Estimated')
title('Comparison to Prevalence of Latent TB from Jordan')
legend('Location','SouthWest')
hold off
% plot(years_aj+2, prevalence_aj)
saveas(gcf,['Population vs Time of ErrorMinExp with Jordan Prevalence',imageName, '.png'])


%% 


