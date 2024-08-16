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

imageName = 'Run63';
load('XX63.mat')
load('paramcell63.mat')
load('ReportedTB20062020.mat');

%% code initializations

ERR_TOL = 0.5;
numbins = 10;


close all;

%% extract data
NumSims=size(XX,1);
numSims = NumSims;


allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);
allq3 = zeros(NumSims,1);
alls = zeros(NumSims,1);

allX0 = zeros(NumSims,1);
allE0 = zeros(NumSims,1);
allL0 = zeros(NumSims,1);
allT0 = zeros(NumSims,1);
allR0 = zeros(NumSims,1);

% incidence proportions from Ng and Ricks

allPR = zeros(NumSims,1);
allPLR = zeros(NumSims,1); % all incidence proportion from latent reactivation (L or R)
allPR0 = zeros(NumSims,1); % above but only using initial populations

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
    % don't fit to u
    % u = BPk(8);

    % fit to u
    xmin = XX{k,2};
    u=xmin(4);


    alls(k) = u;
    % PR_timeseries = u*XELTRk(:,5)./(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5));
    PR=  sum(u*XELTRk(:,5))/sum(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) ;
    allPR(k) = PR;


% initial conditions incidence proportion relapse
    allX0(k)=currentx(5);
    allE0(k) = currentx(6);
    allL0(k) = currentx(7);
    allT0(k) = currentx(8);
    allR0(k) = currentx(9);
 
    PR0 = sum(v*XELTRk(1,3)+u*XELTRk(1,5))/sum(p*w*XELTRk(1,2)+v*XELTRk(1,3)+u*XELTRk(1,5)) ;
    allPR0(k) = PR0;

    % Rick proportion from latent reactivation

    PLR = sum(v*XELTRk(:,3)+u*XELTRk(:,5))/sum(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) ;
    allPLR(k) = PLR;



    
   % subplot(2,2,3)
   % plot(PR_timeseries)


end






for k=1:NumSims
    currentx = XX{k,2};



end



%% Get indices based on errors

% numSims = size(XX,1);

errors = zeros(numSims,1);
errors_incidence = zeros(numSims,1);
errors_PR = zeros(numSims,1);
for k=1:numSims
    errors(k) = XX{k,7};


    % incidence error
    EstimatedIncidence = XX{k,5};
    errors_incidence(k) = norm((EstimatedIncidence(1:end-1)'-ReportedIncidence)./ReportedIncidence);

    % PR error
    BP = XX{k,1};
    XELTR = XX{k,4};
    power_relapse = 4;

    xmin2 = XX{k,2};
    BP(9) = xmin2(1); %q1
    BP(10) = xmin2(2); %q2
    BP(11) = xmin2(3); %q3
    BP(8) = xmin2(4); % s
    errors_PR(k) = error_PR(XELTR, BP, power_relapse);
    % err = IncidenceError3(x, BP, initialConditions, ImmigrationRate, ReportedIncidence,ReportedPrevalence);
end





indices_smallerror = find(errors<ERR_TOL);

% try finding cases where sigma is appropriately large and error is small



figure('units','normalized','outerposition',[0 0 1 1])

subplot(1,2,1)
semilogy(sort(errors))
hold on
horizline = 0*sort(errors)+ERR_TOL;
% semilogy(1:length(errors), horizline )
yline(ERR_TOL, '--', ['Cutoff is ERR TOL =', num2str(ERR_TOL)], ...
    'FontSize', 14, ...
    'LabelHorizontalAlignment', 'right');
title('Error vs experiment, log')

subplot(1,2,2)
plot(sort(errors(indices_smallerror)))
title(['Error vs experiment, cutoff ', num2str(ERR_TOL)])

saveas(gcf,['Error vs experiment ',imageName, '.png'])


display(['Among ', num2str(numSims),' simulations, ', num2str(length(indices_smallerror)), ' had errors less than ', num2str(ERR_TOL)])
display(['The standard deviation of error is ',num2str(std(errors(indices_smallerror)))])
display(['The min and max error are respectively ', num2str(min(errors(indices_smallerror))), ' ', num2str(max(errors(indices_smallerror))) ])

%% Houben comparison
% Compare q1 and q2, our estimated prevalence among new immigrants, against Houben's global prevalence of latent TB


LTBI_newImmigrants = allq1+allq2 + allq3;
% EoverL_newImmigrants = allq1./LTBI_newImmigrants    ;


figure('units','normalized','outerposition',[0 0 1 1])


title('Comparing new immigrants to Houben estimates of global prevalence')
subplot(2,2,1)

% histfit(allq1,numbins)
histfit(LTBI_newImmigrants(indices_smallerror),numbins)
title('Prevalence among new immigrants (%)')
% average low high  21.8046423	17.28412035	28.60442942
Houben_mean = 21.80/100;
Houben_low = 17.28/100;
Houben_high = 28.60/100;

Houben_EoverL = 0.036036036;


% Add vertical lines at the mean and +/- 1 standard deviation
xline(Houben_mean, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houben_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houben_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);

Houbenq1 = 0.725/100;

Houbenq1_low = 0.563/100;
Houbenq1_high = 0.920/100;


subplot(2,2,2)
histfit(allq1(indices_smallerror),numbins)
title('Infected within past 2 years (%)')
xline(Houbenq1, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);


subplot(2,2,3)


% A1 = min(allq1(indices_smallerror));
A1 = 0;
% A1 = 0;
A2 =  max(allq1(indices_smallerror));



scatter(allq1(indices_smallerror),allq2(indices_smallerror)+allq3(indices_smallerror),'DisplayName','(qL+qR) vs qE')
hold on
title('Comparing estimates of prevalence among new immigrants to Houben estimates')

plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'color','r','LineStyle','--','LineWidth',3,'DisplayName','Houben total prevalence')
plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'LineStyle',':','LineWidth',3,'DisplayName','Houben proportion recently infected')
legend('Location','NorthEast')

plot( [A1,A2], [Houben_low-A1,Houben_low-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')
plot( [A1,A2], [Houben_high-A1,Houben_high-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')


hold off

subplot(2,2,4)

histfit(alls(indices_smallerror), numbins);
title('Distribution of reactivation rate \sigma')




saveas(gcf,['Houben comparison to our model ',imageName, '.png'])




%% get indices based on Houben's estimates

lower_cutoff = Houben_low;
% upper_cutoff = (Houben_high+Houben_mean)/2;
upper_cutoff = Houben_high;

TOL = 1.0;
ERR_TOL2 = 2;
indices_Houben = find(LTBI_newImmigrants<=upper_cutoff & LTBI_newImmigrants>=lower_cutoff & errors<ERR_TOL2);

% indices_Houben=good_indices;

% indices_Houben(4) = [];

%% plot comparisons to Houben

load('Error_aj.mat');
Years2 = [Years,2021];

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
    % current_index = indices_Houben(k);
    % plot(Years2, XX{current_index,6} )
    current_index = indices_Houben(k);
 
    XELTRk = XX{current_index,4};

    plot(Years2, XELTRk(:,4))
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
% errorbar(years_aj, prevalence_aj, prevalence_aj-LL_aj, UL_aj-prevalence_aj,'DisplayName','Jordan Prevalence')
title('Early latent vs time')
hold on
for k=1:length(indices_Houben)
    current_index = indices_Houben(k);
 
    XELTRk = XX{current_index,4};

    plot(Years2, XELTRk(:,2))
end


saveas(gcf,['Incidence Prevalence Jordan vs time, Sigma dist',imageName, '.png'])




%% Histograms of optimal qELR, sigma; and the proportion of relapse

good_indices = indices_smallerror;
% good_indices = indices_Houben;

numbins=10;
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
histfit(allq1(good_indices), numbins);
title('Distribution of optimal q_E')

subplot(2,2,2)
histfit(allq2(good_indices), numbins);
title('Distribution of optimal q_L')

subplot(2,2,3)
histfit(allq3(good_indices), numbins);
title('Distribution of optimized q_R')

subplot(2,2,4)
histfit(alls(good_indices), numbins);
title('Distribution of s, relapse rate')

saveas(gcf,['Histograms of optimal parameters ',imageName, '.png'])

% numbins=10;



%% Histograms of optimal Y



numbins=10;
figure('units','normalized','outerposition',[0 0 1 1])
set(gca,'FontSize',5)

subplot(2,3,1)
histogram(allX0(good_indices), numbins);
title('Distribution of optimal X_0')

subplot(2,3,2)
histogram(allE0(good_indices), numbins);
title('Distribution of optimal E_0')

subplot(2,3,3)
histogram(allL0(good_indices), numbins);
title('Distribution of optimized L_0')

subplot(2,3,4)
histogram(allT0(good_indices), numbins);
title('Distribution of optimized T_0')

subplot(2,3,5)
histogram(allR0(good_indices), numbins);
title('Distribution of optimized R_0')



saveas(gcf,['histogram of optimal initial conditions ',imageName, '.png'])

numbins=10;


%% plot histogram of incidence proportions
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
histogram(allPR(indices_smallerror), numbins);
title('Proportion of incidence from relapse (%)')
xline(NgRelapseFraction, '--', {'Relapse Proportion'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);

subplot(2,2,2)
histogram(allPR0(good_indices), numbins);
title('Proportion of incidence from relapse at year 0 (%)')

% Rick data
subplot(2,2,3)
histogram(allPLR(good_indices), numbins);
title('Distribution of incidence from latent reactivation (%)')

saveas(gcf,['histogram of incidence proportions',imageName, '.png'])

% %% Proportion of incidence histograms
% figure('units','normalized','outerposition',[0 0 1 1])
% 
% histfit(allPR(indices_smallerror), numbins);
% title('Proportion of incidence from relapse (%)')
%     NgRelapseFraction = 110/1120;
% xline(NgRelapseFraction, '--', {'Relapse Proportion'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);
% 
% saveas(gcf,['Histograms of relapse proportion ',imageName, '.png'])
% 
% 

%% 

figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
scatter(alls(good_indices),allq1(good_indices))
title('q_E vs relapse \sigma')
xlabel('relapse \sigma')
ylabel('q_E')

subplot(2,2,2)
scatter(alls(good_indices),allq2(good_indices))
title('late latent among new immigrants q_L vs relapse \sigma')
xlabel('relapse \sigma')
ylabel('q_L')

subplot(2,2,3)
scatter(alls(good_indices),allq3(good_indices))
title('recovered among new immigrants q_R vs relapse \sigma')
xlabel('relapse \sigma')
ylabel('q_R')

saveas(gcf,['Scatter sigma vs qELR ',imageName, '.png'])





%% Error minimizing experiment's exit info, gradient, hessians

% find error minimiizing and maximizing experiments
[error_min,min_index] = min(errors_incidence);
[error_max,max_index] = max(errors_incidence);

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
BPi_min(8) = x_min(4);

    

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



