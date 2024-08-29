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
set(groot,'defaultAxesFontSize',32)
set(groot,'defaultLegendFontSize',32)


addpath('functions\')
addpath('figure creation\')

addpath('data and results\')

%%

runNum = '93';
imageName = runNum;
xxname = ['XX',runNum];
load([xxname,'.mat'])
paramcellname = ['paramcell',runNum];
load([paramcellname,'.mat'])
load('ReportedTB20062020.mat');

%% code initializations

ERR_TOL = 1;
numbins = 10;

% E over total TB in 2021 should be close to Jordan
JordanEoverTB_mean = 1/488;
JordanEoverTB_high =1/185;
JordanEoverTB_low = 1/1039;

Houben_mean = 21.80/100;
Houben_low = 17.28/100;
Houben_high = 28.60/100;

Houben_EoverL = 0.036036036;

Houbenq1 = 0.725/100;
Houbenq1_low = 0.563/100;
Houbenq1_high = 0.920/100;

RickReactivation = 83.7/100;

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

% error
errors = zeros(numSims,1);
errors_incidence = zeros(numSims,1);
errors_PR = zeros(numSims,1);
errors_E2021 = zeros(numSims,1);

% exit flags
allexitflag = zeros(NumSims,1);

% comparisons to Jordan
numYears = length(Years)+1;
allPrevalenceLTBI = zeros(NumSims,numYears);
EtoAllTB2021 = zeros(NumSims,1);
EtoAllTB0 = zeros(NumSims,1);


% incidence proportions from Ng 110/1120 

allPR = zeros(NumSims,1);
allPR0 = zeros(NumSims,1); % above but only using initial populations

% incidence proportion among new immigrants
all_qPR =  zeros(NumSims,1);

% Ricks 83.7%
% incidence proportion from latent reactivation (L or R)
allPLR = zeros(NumSims,1); % all incidence proportion from latent reactivation (L or R)
allPLR0 = zeros(NumSims,1);
allJE0 = zeros(NumSims,1);

% E2021 = E0 + integral, E2021/(E+L+T+R2021) = 1/488
allE2021ratio = zeros(NumSims,1);

for k=1:NumSims
    currentparams = XX{k,1};
    currentx = XX{k,2};
    XELTRk = XX{k, 4};

    % currentx is a y; unscale
    % currentxm = [BPk(9) BPk(10) BPk(11) BPk(8), BPk(12:end)]
    % currentx = currentx.*currentxm;


    allq1(k)=currentx(1);
    allq2(k) = currentx(2);
    allq3(k) = currentx(3);
    alls(k) = currentx(4);

    % initial conditions incidence proportion relapse
    % allX0(k)=currentx(5);
    % allE0(k) = currentx(6);
    % allL0(k) = currentx(7);
    % allT0(k) = currentx(8);
    % allR0(k) = currentx(9);
 
    allX0(k) = XELTRk(1,1);
    allE0(k) = XELTRk(1,2);
    allL0(k) = XELTRk(1,3);
    allT0(k) = XELTRk(1,4);
    allR0(k) = XELTRk(1,5);


    % compute incidence proportions
    BPk = XX{k,1};
    BPk(8) = currentx(4); % update sigma


    p = BPk(2);
    w = BPk(3);
    v = BPk(4);
    % don't fit to u
    % u = BPk(8);

    
    u=currentx(4);

    % PR_timeseries = u*XELTRk(:,5)./(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5));
    % PR=  sum(u*XELTRk(:,5))/sum(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) ;
    PR = IncidenceRelapseProportion(XELTRk,BPk);
    allPR(k) = PR;


    % PR0 = sum(v*XELTRk(1,3)+u*XELTRk(1,5))/sum(p*w*XELTRk(1,2)+v*XELTRk(1,3)+u*XELTRk(1,5)) ;
    PR0 = IncidenceRelapseProportion(XELTRk(1,:),BPk);
    allPR0(k) = PR0;


    XELTRk_q = XELTRk;
    XELTRk_q(2) = allq1(k);
    XELTRk_q(3) = allq2(k);
    XELTRk_q(5) = allq3(k);
    % q_PR = IncidenceRelapseProportion(XELTRk_q,BPk);
    % q_PR = IncidenceRelapseProportion([0 allq1(k) allq2(k) 0 allq3(k)],BPk);
    q_PR = u*allq3(k)/(p*w*allq1(k)+v*allq2(k) + u*allq3(k));
    all_qPR(k) = q_PR;

    % Rick proportion from latent reactivation

    PLR = sum(v*XELTRk(:,3)+u*XELTRk(:,5))/sum(p*w*XELTRk(:,2)+v*XELTRk(:,3)+u*XELTRk(:,5)) ;
    allPLR(k) = PLR;

    PLR0 = sum(v*XELTRk(1,3)+u*XELTRk(1,5))/sum(p*w*XELTRk(1,2)+v*XELTRk(1,3)+u*XELTRk(1,5)) ;
    allPLR0(k) = PLR0;

    % Jordan2021 E2021/(E+L+T+R) = 1/488

    % Jordan E2021 ratio, but applied for E0
    JE0 = XELTRk(1,2)/sum(XELTRk(1,[2:5]));
    allJE0(k) = JE0;


    % E2021_num = sum(XELTRk(:,2));
    % E2021_denom = sum(XELTRk(:,2))+sum(XELTRk(:,3))+sum(XELTRk(:,4))+sum(XELTRk(:,5));
    % allE2021ratio(k) = E2021_num/E2021_denom;

    
   % subplot(2,2,3)
   % plot(PR_timeseries)

   XX3 = XX{k,3};

   allexitflag(k) = XX3{1};

   % prevalence of LTBI
   allPrevalenceLTBI(k,:) = (sum(XELTRk(:,[2:5])')'./ (sum(XELTRk(:,:)')') )';
    %make into percentage
    % totalPops = (sum(XELTRk(:,:)')');
% allPrevalenceLTBI(k) = allPrevalenceLTBI(k)./totalPops;
    EtoAllTB2021(k) = XELTRk(end,2)/sum(XELTRk(end,[2:5]));
    EtoAllTB0(k) = XELTRk(1,2)/sum(XELTRk(1,[2:5]));
    
    power_relapsek = currentparams(17);
    % EtoAllTB2021(k) = error_E2021(XELTR,  power_relapsek);

   % errors
    errors(k) = XX{k,7};


    % incidence error
    EstimatedIncidence = XX{k,5};
    errors_incidence(k) = norm((EstimatedIncidence(1:end-1)'-ReportedIncidence)./ReportedIncidence);

    % PR error
    BP = XX{k,1};
    XELTR = XX{k,4};
    

    xmin2 = XX{k,2};
    BP(9) = xmin2(1); %q1
    BP(10) = xmin2(2); %q2
    BP(11) = xmin2(3); %q3
    BP(8) = xmin2(4); % s
    errors_PR(k) = error_PR(XELTR, BP, power_relapsek);
    % err = IncidenceError3(x, BP, initialConditions, ImmigrationRate, ReportedIncidence,ReportedPrevalence);

    % errors Jordan E2021
    errors_E2021(k) = error_E2021(XELTR,power_relapsek);
end





%% Get indices based on errors and other checks

% numSims = size(XX,1);


indices_smallerror = find(errors<ERR_TOL);

% try finding cases where sigma is appropriately large and error is small



% exit flag is one
indices_exitflag1 = find(allexitflag==1);
% good_indices = exitflag1;

% Houben indices

LTBI_newImmigrants = allq1+allq2 + allq3;

Houben_mean = 21.80/100;
Houben_low = 17.28/100;
Houben_high = 28.60/100;

Houben_EoverL = 0.036036036;


TOL = 1.0;
% ERR_TOL2 = 2;
% indices_Houben = find(LTBI_newImmigrants<=upper_cutoff & LTBI_newImmigrants>=lower_cutoff &  errors<ERR_TOL);



% begin with all indices are good. turn off ones that fail

good_indices = ones(numSims,1);

good_indices(find(errors>=ERR_TOL)) = 0;
% exitflag correct
% good_indices(find(allexitflag==0)) = 0;
% good_indices(find(allexitflag==2)) = 0;

% allPLR >90% eliminate
% good_indices(find(allPLR>0.90)) = 0;


%allPR0 < 8% eliminate
% good_indices(find(allPR0<8/100)) = 0;
% > 12% eliminate
% good_indices(find(allPR0>12/100)) = 0;


% if annual LTBI prevalence ever < Jordan min or > Jordan max, eliminate
% kcount = 0;
% Jordan_min = min(LL_aj);
% Jordan_max  = max(UL_aj);
% for k=1:numSims
%     % only eliminate if less than 
% 
%     prevalenceLTBI = allPrevalenceLTBI(k,:);
%     if min(prevalenceLTBI) < Jordan_min/100
%       good_indices(k) = 0;
%       kcount=kcount+1;
%     elseif max(prevalenceLTBI) > Jordan_max/100
%       good_indices(k) = 0;
%       kcount=kcount+1;
%     end
% end


% E over total TB in 2021 should be close to Jordan 2021
good_indices(find(EtoAllTB2021<JordanEoverTB_low)) = 0;
good_indices(find(EtoAllTB2021>JordanEoverTB_high)) = 0;

% good_indices(find(EtoAllTB0(k) <JordanEoverTB_low/2)) = 0;
% good_indices(find(EtoAllTB0>JordanEoverTB_high*2)) = 0;


% E over total TB in 2001 should be close to Jordan 2021
% E2021_upper = 1/185;
% E2021_lower = 1/1039;
%  good_indices(find(allJE0>E2021_upper))=0;
%  good_indices(find(allJE0<E2021_lower))=0;

% proportion incidence from latent should be close to Ricks US count
% RickDouble = 1-(1-RickReactivation)*2;
 % good_indices(find(allPLR0<RickDouble))=0;
 % good_indices(find(allPLR0<0.6))=0;


 % Houben qE+qL+qR
upper_cutoff = Houben_high;
lower_cutoff = Houben_low;
% new immigrants LTBI comparable to Houben
good_indices(find(LTBI_newImmigrants>upper_cutoff ))=0;
good_indices(find(LTBI_newImmigrants<lower_cutoff ))=0;
% good_indices(find(LTBI_newImmigrants>0.5))=0;
% 
% good_indices(find(LTBI_newImmigrants<0.05))=0;

% Houben qE 
good_indices( find(allq1>Houbenq1_high) )=0;
good_indices( find(allq1<Houbenq1_low) )=0;


good_indices = find(good_indices);


%% export data from good indices

numgood = length(good_indices);

numparams = length(XX{1,1});
numoutput = length(XX{1,2});

GOODPARAMS = zeros(numgood,numparams);
GOODOUTPUT = zeros(numgood,numoutput);
for k=1:numgood
    GOODPARAMS(k,:) = XX{good_indices(k),1};
    GOODOUTPUT(k,:) = XX{good_indices(k),2};
end
% save('GOODPARAMS.mat','GOODPARAMS','GOODOUTPUT')

% compute prevalence amplification from 2011 onward

prev_amplification = zeros(numgood,1);

warmup_index = 5; % first 5 years are warmup

for k=1:numgood 
    prevalence_estimate = XX{good_indices(k),6};
    


    prev_amplification(k) = mean( prevalence_estimate(warmup_index+1:end-1)./ReportedPrevalence(warmup_index+1:end)' );
end

display(['with recovery rate d=',num2str(GOODPARAMS(1,6)), ' the amplification is ', mean(num2str(prev_amplification(1:2)))])
display(['with recovery rate d=',num2str(GOODPARAMS(3,6)), ' the amplification is ', mean(num2str(prev_amplification(3:4)))])


%% 

numparam = length(XX{1,1});
validparams = zeros(length(good_indices),numparam);

numx=length(XX{1,2});
validoutputs = zeros(length(good_indices),numx);


for k=1:length(good_indices)
    currentparam = XX{good_indices(k),1};
    validparams(k,:) = currentparam;

    currentoutput = XX{good_indices(k),2};
    validoutputs(k,:) = currentoutput;
end
colLabels = {'beta', 'p','w','v','a','d','n','sigma','qE','qL','qR','X0','E0','L0','T0','R0','prevamp','power relapse'};
my_matrix2latex(validparams, 'validparams.tex',  'columnLabels',colLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
colLabels = {'qE','qL','qR', '\sigma','X0','E0','L0','T0','R0'};
my_matrix2latex(validoutputs, 'validoutputs.tex',  'columnLabels',colLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

%% plot error vs experiment


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
histfit(LTBI_newImmigrants(good_indices),numbins)
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
histfit(allq1(good_indices),numbins)
title('Infected within past 2 years (%)')
xline(Houbenq1, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(Houbenq1_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);


subplot(2,2,3)


% A1 = min(allq1(indices_smallerror));
A1 = 0;
% A1 = 0;
A2 =  max(allq1(good_indices));



scatter(allq1(good_indices),allq2(good_indices)+allq3(good_indices),'DisplayName','(qL+qR) vs qE')
hold on
title('Comparing estimates of prevalence among new immigrants to Houben estimates')

plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'color','r','LineStyle','--','LineWidth',3,'DisplayName','Houben total prevalence')
plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'LineStyle',':','LineWidth',3,'DisplayName','Houben proportion recently infected')
legend('Location','NorthEast')

plot( [A1,A2], [Houben_low-A1,Houben_low-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')
plot( [A1,A2], [Houben_high-A1,Houben_high-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')


hold off

subplot(2,2,4)

histfit(alls(good_indices), numbins);
title('Distribution of reactivation rate \sigma')




% Comparisons to Houben
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,1,1)
title('Comparing new immigrants to Houben estimates of global prevalence')
hold on

% histfit(allq1,numbins)
histogram(LTBI_newImmigrants(good_indices),'BinWidth',0.02,'DisplayName','LTBI prevalence')
title('LTBI prevalence among new immigrants (%)')
xlabel('%')


% Add vertical lines at the mean and +/- 1 standard deviation
xline(Houben_mean, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houben_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houben_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')


% also plot proportion of recently infected
histogram(allq1(good_indices),'BinWidth',0.002,'DisplayName','LTBI % recently infected')

xline(Houbenq1, '-','LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houbenq1_low, '--', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houbenq1_high, '--',  'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')

legend('Location','NorthWest')



subplot(2,1,2)


% A1 = min(allq1(indices_smallerror));
A1 = 0;
% A1 = 0;
A2 =  max(allq1(good_indices));



scatter(allq1(good_indices),allq2(good_indices)+allq3(good_indices),'DisplayName','(qL+qR) vs qE')
hold on
title('LTBI prevalence among new immigrants')
xlabel('LTBI recently infected (%)')
ylabel('LTBI not recently infected (%)')
plot( [A1,A2], [Houben_mean-A1,Houben_mean-A2], 'color','r','LineStyle','--','LineWidth',3,'DisplayName','Houben total prevalence')
plot([A1,A2], (1-Houben_EoverL)*[A1,A2]/Houben_EoverL,'LineStyle',':','LineWidth',3,'DisplayName','Houben proportion recently infected')
legend('Location','SouthEast')

plot( [A1,A2], [Houben_low-A1,Houben_low-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')
plot( [A1,A2], [Houben_high-A1,Houben_high-A2], 'color','r','LineStyle','--','LineWidth',1,'HandleVisibility','off')

saveas(gcf,['Houben comparison to our model ',imageName, '.png'])

%% calibration cube: compare to Houben LTBI and Jordan's E2021 ratio


% Comparisons to Houben
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,1,1)
title('Comparing new immigrants to Houben estimates of global prevalence')
hold on

% histfit(allq1,numbins)
histogram(LTBI_newImmigrants(good_indices),'BinWidth',0.02,'DisplayName','LTBI prevalence')
title('LTBI prevalence among new immigrants (%)')
xlabel('%')


% Add vertical lines at the mean and +/- 1 standard deviation
xline(Houben_mean, '-', {'Houben Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houben_low, '--', {'Houben Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houben_high, '--', {'Houben High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')


% also plot proportion of recently infected
histogram(allq1(good_indices),'BinWidth',0.002,'DisplayName','LTBI % recently infected')

xline(Houbenq1, '-','LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houbenq1_low, '--', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')
xline(Houbenq1_high, '--',  'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15,'HandleVisibility','off')

legend();
% legend('Location','NorthWest')

% Comparisons to Jordan
subplot(2,1,2)


histogram(EtoAllTB2021(good_indices)*100,numbins)
title('Proportion of TB in 2021 that were recently infected (E over total TB)')
xlabel('%')

% Add vertical lines at the mean and +/- 1 standard deviation
xline(JordanEoverTB_mean*100, '-', {'Jordan Mean'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(JordanEoverTB_low*100, '--', {'Jordan Low'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);
xline(JordanEoverTB_high*100, '--', {'Jordan High'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',15);


saveas(gcf,['Calibration ranges Houben Jordan',imageName, '.png'])


%% many plots of incidence and prevalence

load('Error_aj.mat');
Years2 = [Years,2021];

% good_indices = d
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
title('Incidence vs time')
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    plot(Years2, XX{current_index,5} )
end
plot(Years,ReportedIncidence,'k--', 'LineWidth',5)


subplot(2,1,2)
title('Prevalence vs time')
hold on
for k=1:length(good_indices)
    % current_index = indices_Houben(k);
    % plot(Years2, XX{current_index,6} )
    current_index = good_indices(k);
 
    XELTRk = XX{current_index,4};

    plot(Years2, XELTRk(:,4))
end
plot(Years,ReportedPrevalence,'k--', 'LineWidth',5)

% ylim([min(ReportedPrevalence)/3,max(ReportedPrevalence)*5])

saveas(gcf,['Incidence Prevalence vs time, many plots',imageName, '.png'])

%% XELR vs time

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
title('susceptible X vs time')
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    XELTRk = XX{current_index,4};
    plot(Years2, XELTRk(:,1) )
end

subplot(2,2,2)
title('early latent E vs time')
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    XELTRk = XX{current_index,4};
    plot(Years2, XELTRk(:,2) )
end

subplot(2,2,3)
title('late latent L vs time')
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    XELTRk = XX{current_index,4};
    plot(Years2, XELTRk(:,3) )
end

subplot(2,2,4)
title('recovered R vs time')
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    XELTRk = XX{current_index,4};
    plot(Years2, XELTRk(:,5) )
end

saveas(gcf,['XELR vs time, many plots',imageName, '.png'])


%% how much more is prevalence than reported prevalence?

% prev_amp=zeros(length(good_indices),1);
% for k=1:length(good_indices)
%     current_index = good_indices(k);
% 
%     XELTRk = XX{current_index,4};
% 
%     prev_amp(k) = mean( XELTRk(1:end-1,4)./ReportedPrevalence');
% end


%% comparison to Jordan 
figure('units','normalized','outerposition',[0 0 1 1])

% subplot(2,1,1)
title(['LTBI prevalence, compared to Jordan''s estimate'])
hold on
for k=1:length(good_indices)
    current_index = good_indices(k);
    estimate_Paj = zeros(length(Years2),1 ) ;

    for j=1:length(Years2)
        XELTR = XX{current_index,4};
        estimate_Paj(j) = sum(XELTR(j,[2:3,5]))./sum(XELTR(j,:));
    end

    plot(Years2, estimate_Paj*100)
end
errorbar(years_aj, prevalence_aj, prevalence_aj-LL_aj, UL_aj-prevalence_aj,'DisplayName','Jordan Prevalence')

saveas(gcf,['Jordan LTBI ',imageName, '.png'])





%% Histograms of optimal qELR, sigma; and the proportion of relapse



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
title('Distribution of optimal q_R')

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
xlim([0,max(allE0(good_indices))])

subplot(2,3,3)
histogram(allL0(good_indices), numbins);
title('Distribution of optimized L_0')

subplot(2,3,4)
histogram(allT0(good_indices), numbins);
title('Distribution of optimized T_0')

subplot(2,3,5)
histogram(allR0(good_indices), numbins);
title('Distribution of optimized R_0')
xlim([0,max(allR0(good_indices))])


saveas(gcf,['histogram of optimal initial conditions ',imageName, '.png'])

numbins=10;


%% plot histogram of incidence proportions
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
histogram(allPR(good_indices), numbins);
title('Proportion of incidence from relapse (%)')
NgRelapseFraction = 110/1120;
xline(NgRelapseFraction, '--', {'Relapse %'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);

subplot(2,2,2)
histogram(allPR0(good_indices), numbins);
title('at year 0 (%)')
xline(NgRelapseFraction, '--', {'Relapse %'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);

% Rick data = 83.7
RickReactivation = 83.7/100;

subplot(2,2,3)
histogram(allPLR(good_indices), numbins);
title('Proportion of incidence from latent reactivation (%)')
xline(RickReactivation, '--', {'Reactivation %'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);


subplot(2,2,4)
histogram(allPLR0(good_indices), numbins);
title('at year 0 (%)')
xline(RickReactivation, '--', {'Reactivation %'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);

saveas(gcf,['Incidence Proportion Histograms',imageName, '.png'])

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

%% plot histogram of incidence proportions using q
figure('units','normalized','outerposition',[0 0 1 1])

histogram(all_qPR(good_indices), numbins);
title('Relapse proportion among new immigrants (%)')
NgRelapseFraction = 110/1120;
xline(NgRelapseFraction, '--', {'Relapse %'}, 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'FontSize',24);

saveas(gcf,['Incidence Proportion New Immigrants ',imageName, '.png'])






%% Error minimizing experiment's exit info, gradient, hessians
% 
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




