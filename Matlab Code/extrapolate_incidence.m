
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultContourLineWidth',1.5)
set(groot,'DefaultFunctionContourLineWidth',2)
set(groot,'defaultAxesFontSize',32)
set(groot,'defaultLegendFontSize',32)


addpath('functions\')
addpath('figure creation\')

addpath('data and results\')

%% extrapolate immigration



load('MoreImmigration.mat')
p = polyfit(Years20062022,ReportedImmigration20062022,1);


Years_long = 2006:2035;
numt_long = length(Years_long);
pi_extrapolated = zeros(1,numt_long);

for k=1:numt_long
    pi_extrapolated(k) = p(1)*Years_long(k)+p(2);
end

figure('units','normalized','outerposition',[0 0 1 1])

plot(Years20062022,ReportedImmigration20062022)
hold on
plot(Years_long, pi_extrapolated)
saveas(gcf,['Extrapolated immigration rate.png'])



save('pi_extrapolated.mat','pi_extrapolated');


load('qvec.mat')



%% extrapolate incidence



load('GOODPARAM.mat')


bioParameters = GOODPARAM(1:8);
bioParameters(8) = GOODOUTPUT(4);
IC = GOODOUTPUT(5:9);





%% Estimate incidence in 2030
[XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2);

% load incidence data to plot
load('ReportedTB20062020.mat');


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
title('Incidence vs time')
plot(Years_long, TBIncidence_long(1:end-1)')
hold on

plot(Years,ReportedIncidence,'k--', 'LineWidth',5)

subplot(2,1,2)
title('Prevalence vs time')
plot(Years_long, TBPrevalence_long(1:end-1)')

hold on

plot(Years,ReportedPrevalence,'k--', 'LineWidth',5)

saveas(gcf,['Incidence Prevalence extrapolated.png'])

%% compare Houben's mean burden among new immigrants vs fitted to incidence

qvec2 = zeros(numt_long,4);

qvec2(1:15,:) = qvec_t(6:end,:);
for k=16:numt_long
    qvec2(k,:) = qvecbar;
end


figure('units','normalized','outerposition',[0 0 1 1])

plot(Years_long, qvec2)
hold on
qE = GOODOUTPUT(1);
qL = GOODOUTPUT(2);
qR = GOODOUTPUT(3);
qX = 1-qE-qL-qR;
plot(Years_long, ones(size(Years_long))*qX,'--')
plot(Years_long, ones(size(Years_long))*qE,'--')
plot(Years_long, ones(size(Years_long))*qL,'--')
plot(Years_long, ones(size(Years_long))*qR,'--')


title('q_X, q_E, q_L, q_R vs time')
legend('q_X','q_E','q_L','q_R', 'q_X optimized','q_E optimized', 'q_L optimized', 'q_R optimized')

saveas(gcf,['qvec vs time.png'])