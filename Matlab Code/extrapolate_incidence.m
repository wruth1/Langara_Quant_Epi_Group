
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultContourLineWidth',1.5)
set(groot,'DefaultFunctionContourLineWidth',2)
set(groot,'defaultAxesFontSize',32)
set(groot,'defaultLegendFontSize',32)


addpath('functions\')
% addpath('figure creation\')

addpath('data and results\')
addpath('Global Sensitivity Analysis\')


%% extrapolate immigration



load('MoreImmigration.mat')
regression_coeff = polyfit(Years20062022(1:end-3),ReportedImmigration20062022(1:end-3),1); %truncate 2020,2021,2022


Years_long = 2006:2035;
numt_long = length(Years_long);
pi_linear = zeros(1,numt_long);

for k=1:numt_long
    pi_linear(k) = regression_coeff(1)*Years_long(k)+regression_coeff(2);
end

pi_extrapolated = zeros(numt_long,1);
indexcounter = 0;
for t=Years_long
    indexcounter = indexcounter+1;
    pi_extrapolated(indexcounter) = pi_extra(t);
end

%% plot 
figure('units','normalized','outerposition',[0 0 1 1])

title('Immigration vs time')
hold on

plot(Years_long, pi_linear,'ro--', 'DisplayName','linear')
plot(Years_long,pi_extrapolated,'bSquare-','DisplayName','extrapolated')
plot(Years20062022,ReportedImmigration20062022,'k--','DisplayName','reported','LineWidth',5)
legend()
saveas(gcf,['Extrapolated immigration rate.png'])

hold off
%%
% pi_extrapolated = zeros(length(Years_long),1);
% numdata = length(ReportedImmigration20062022);
% pi_extrapolated(1:numdata) = ReportedImmigration20062022;
% pi_extrapolated(numdata+1:end) = pi_linear(numdata+1:end);

pi_extrapolated = zeros(numt_long,1);
indexcounter = 0;
for t=Years_long
    indexcounter = indexcounter+1;
    pi_extrapolated(indexcounter) = pi_extra(t);
end


% pi_extrapolated = pi_extra()
plot(Years_long,pi_extrapolated')

save('pi_extrapolated.mat','pi_linear');


load('qvec.mat')



%% extrapolate incidence



load('GOODPARAM.mat')


bioParameters = GOODPARAM(1:8);
bioParameters(8) = GOODOUTPUT(4);
IC = GOODOUTPUT(5:9);

Params = [bioParameters, IC];
scale_dP = 10/100;

% ParamMesh = createUniformParameters(Params,scale_dP);



%% Estimate incidence in 2030
[XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_linear,qvec2);

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