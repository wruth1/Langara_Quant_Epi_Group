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


addpath('functions\')
addpath('figure creation\')

addpath('data and results\')

%%

load('XX50.mat')
load('paramcell50.mat')
load('ReportedTB20062020.mat');

%%

imageName = 'Run50';

NumSims=size(XX,1);
numSims = NumSims;


close all;

%% Histograms of optimal

allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);
allq3 = zeros(NumSims,1);
alls = zeros(NumSims,1);
allPR = zeros(NumSims,1);



for k=1:NumSims
    currentx = XX{k,2};
    allq1(k)=currentx(1);
    allq2(k) = currentx(2);
    allq3(k) = currentx(3);
    % alls(k) = XX{i,9};

    % also extract allq3

    % currentparams = XX{i,1};
    % allq3(i) = currentparams(11);

    BPk = XX{k,1};
    XELTRk = XX{k, 4};

    p = BPk(2);
    w = BPk(3);
    v = BPk(4);
    u = BPk(8);

    % PR_timeseries = u*XELTRk(:,5)./(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5));
    PR=  sum(u*XELTRk(:,5))/sum(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) ;

    % PR_mean = mean(PR_timeseries);

   % proportionRelapse = mean( u*XELTRk(:,5)./(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) );


   allPR(k) = PR;

end

numbins=10;
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
histfit(allq1, numbins);
title('Distribution of optimal q_E')

subplot(2,2,2)
histfit(allq2, numbins);
title('Distribution of optimal q_L')

subplot(2,2,3)
histfit(allq3, numbins);
title('Distribution of optimized q_R')

subplot(2,2,4)
histfit(allPR, numbins);
title('Distribution of relapse proportion')

saveas(gcf,['Histograms of optimal parameters ',imageName, '.png'])


%% 

figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)
scatter(alls,allq1)
title('q_E vs relapse \sigma')
xlabel('relapse \sigma')
ylabel('q_E')


subplot(1,2,2)
scatter(alls,allq3)
title('q_R vs relapse \sigma')
xlabel('relapse \sigma')
ylabel('q_R')

saveas(gcf,['Scatter sigma vs qE qR ',imageName, '.png'])


%% Houben comparison
% Compare q1 and q2, our estimated prevalence among new immigrants, against Houben's global prevalence of latent TB


LTBI_newImmigrants = allq1+allq2 + allq3;
% EoverL_newImmigrants = allq1./LTBI_newImmigrants    ;


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
scatter(allq1,allq2+allq3,'DisplayName','q2 vs q1')
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
scatter(allq1,allq2+allq3,'DisplayName','(q2+q3) vs q1')
hold on

A1 = 0;

plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'DisplayName','Houben Mean')
plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'LineStyle','--','DisplayName','Houben q1/(q1+q2+q3)')

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

plot(sort(errors))
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
%% 

BPi_min = XX{min_index,1};
x_min = XX{min_index,2};
BPi_min(9) = x_min(1);
BPi_min(10) = x_min(2);
BPi_min(11) = x_min(3);
% BPi_min(8) = XX{min_index,9};

    

numt = length(ReportedImmigration);


INC = zeros(1,numt);
PREV = zeros(1,numt);
for k=1:numt
    Xi = XELTR_min(k,1);
    Ei = XELTR_min(k,2);
    Li = XELTR_min(k,3);
    Ti = XELTR_min(k,4);
    Ri = XELTR_min(k,5);

    

    INC(k) = BPi_min(2) * BPi_min(3) * Ei + BPi_min(4) * Li + BPi_min(8)*Ri;

    INC(k) = INC(k) * 100000/(Xi+Ei+Li+Ti+Ri);

    PREV(k) = Ti;

end

figure 
subplot(2,1,1)
plot(INC)
hold on
plot(incidence_min,'--')

subplot(2,1,2)
plot(PREV)
hold on
plot(prevalence_min,'--')



%% Plot test data - compare our model to

load('Error_aj.mat');

% compute error Paj
% total_tb = sum(XELTR_min(:,2:4));

Years2 = [Years, 2021];

estimate_Paj = zeros(length(Years2),1 ) ;

for k=1:length(Years2)
    estimate_Paj(k) = sum(XELTR_min(k,[2:3,5]))./sum(XELTR_min(k,:));
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


%% only examine plots that agree with Houben's prevalence estimate

lower_cutoff = Houben_low;
% upper_cutoff = (Houben_high+Houben_mean)/2;
upper_cutoff = Houben_high;

TOL = 1.0;
indices_Houben = find(LTBI_newImmigrants<=upper_cutoff & LTBI_newImmigrants>=lower_cutoff & errors<=TOL);

% indices_Houben(4) = [];

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(Years,ReportedIncidence, 'LineWidth',3)
title('Incidence vs time')
hold on
for k=1:length(indices_Houben)
    current_index = indices_Houben(k);
    plot(Years2, XX{current_index,5} )
end

subplot(2,2,2)
plot(Years,ReportedPrevalence, 'LineWidth',3)
title('Prevalence vs time')
hold on
for k=1:length(indices_Houben)
    current_index = indices_Houben(k);
    plot(Years2, XX{current_index,6} )
end

subplot(2,2,3)
errorbar(years_aj, prevalence_aj, prevalence_aj-LL_aj, UL_aj-prevalence_aj,'DisplayName','Jordan Prevalence')
title('Prevalence of latent TB, compared to Jordan')
hold on
for k=1:length(indices_Houben)
    current_index = indices_Houben(k);
    estimate_Paj = zeros(length(Years2),1 ) ;

    for j=1:length(Years2)
        XELTR = XX{current_index,4};
        estimate_Paj(j) = sum(XELTR(j,[2:3,5]))./sum(XELTR(j,:));
    end

    plot(Years2, estimate_Paj*100)
end

subplot(2,2,4)
histfit(alls(indices_Houben), numbins);
title('Distribution of s, reactivation of R')

saveas(gcf,['Incidence Prevalence Jordan vs time, Sigma dist',imageName, '.png'])

