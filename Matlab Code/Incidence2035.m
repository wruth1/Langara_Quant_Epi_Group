

% TBI_c = zeros(numc,1);
% for c=1:numc
%     TBI_c(c) = WHOincidence(WHOregion_c{c});
% end


addpath('functions\')
addpath('figure creation\')

addpath('data and results\')


load('Reported_piCW_20012020.mat')


%% create qvec 2006-2020



piW_20062020 = piW_20012020(6:end,:);
Years_20062020 = (2006:2020)';

numt_short = length(Years_20062020);

% qvec20062020 = build_qvecW_t(piW_20062020,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w, INF_w,questdlg)
% qvecW_SEA2 = build_qvecW_t(piW_SEA,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );

%% make qvec2_t20062035 by stapling

numq = 4; % number of compartments for q: qX qE qL qR

Years_long = 2006:2035;
numt_long = length(Years_long);
numt_short = length(Years_20062020);
numt_extra = numt_long-numt_short;

%% assume SEA and AMR only

% assume only immigrants from SEA region, highest prevalence of LTBI
qvec_SEA = zeros(numt_long,numq);
% assume only immigrants from AMR region, highest prevalence of LTBI
qvec_AMR = zeros(numt_long,numq);



% load('pi_extrapolated.mat');
% piW_SEA
numw=6;

Years_extra = Years_long(end-numt_extra+1:end);
piW_extra_SEA = zeros(numt_extra,numw);
piW_extra_AMR = zeros(numt_extra,numw);

% SEA is index 5
piW_extra_SEA(:,5) = pi_extra(Years_extra);
piW_extra_AMR(:,2) = pi_extra(Years_extra);

piW_SEA = [piW_20062020; piW_extra_SEA];
piW_AMR = [piW_20062020; piW_extra_AMR];

% create qvec
load('Houben2014_w.mat')
Houben_weight = 0;
qLqR_ratio = 0.46;

% Houben_weight= Houben_weight_experimental;

% qvecW_SEA = build_qvecW_t(piW_SEA,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );
% qvecW_AMR = build_qvecW_t(piW_AMR,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );

DecayRate = 1/100;
qvecW_SEA_decay = build_qvecW_t_decay(Years_long,piW_SEA,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio ,DecayRate );
qvecW_AMR_decay = build_qvecW_t_decay(Years_long,piW_AMR,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio,DecayRate  );



%%  build baseline

qvec_baseline = zeros(numt_long,numq);

PI_W = sum(piW_20062020);
PI_total = sum(PI_W);

% this is vector in R^6. sums to 1. proportion of immigrants from each Who geoegraphic region
vecpi_baseline = PI_W/PI_total; 

piW_extra_base = zeros(numt_extra,numw);

% total number of immigrants is whatever is extrapolated
pi_extrapolated_extra = pi_extra(Years_extra);
% then distribution is always vecpi_baseline
for i=1:numt_extra
    piW_extra_base(i,:) = vecpi_baseline;
end

piW_base = [piW_20062020; piW_extra_base];
% qvecW_base = build_qvecW_t(piW_base,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );
qvecW_base_decay = build_qvecW_t_decay(Years_long,piW_base,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio ,DecayRate );


% save('piW_base.mat','piW_base','qvecW_base','pi_extrapolated');




%% compute incidence in 2035



load('GOODPARAM.mat')


bioParameters = GOODPARAM(1:8);
bioParameters(8) = GOODOUTPUT(4);
IC = GOODOUTPUT(5:9);


% find Houben weight
qLTBI_computed = sum(GOODOUTPUT(1:3));

% observe qLTBI > mean
Houben_mean = 21.80/100;
Houben_low = 17.28/100;
Houben_high = 28.60/100;

Houben_weight_experimental = (qLTBI_computed-Houben_mean)/(Houben_high-Houben_mean);





%%
pi_extrapolated = pi_extra(Years_long);
[XELTR_SEA, TBIncidence_SEA, TBPrevalence_SEA] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvecW_SEA_decay);
% [XELTR_SEA2, TBIncidence_SEA2, TBPrevalence_SEA_decay] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvecW_SEA_decay);

[XELTR_AMR, TBIncidence_AMR, TBPrevalence_AMR] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvecW_AMR_decay);
[XELTR_base, TBIncidence_base, TBPrevalence_base] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvecW_base_decay);
%% plot incidence
% load incidence data to plot
load('ReportedTB20062020.mat');


figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,2,1)
title('Incidence vs time')

hold on
plot(Years_long, TBIncidence_SEA(1:end-1),'b:','DisplayName','SEA')
% plot(Years_long, TBIncidence_SEA2(1:end-1),'k','DisplayName','SEA2')
% plot(Years_long, TBIncidence_SEA3(1:end-1),'ro-')


plot(Years_long, TBIncidence_AMR(1:end-1)','ro--','DisplayName','AMR')
plot(Years_long, TBIncidence_base(1:end-1)','mx--','DisplayName','base')
plot(Years,ReportedIncidence,'k--', 'LineWidth',5,'DisplayName','reported')
incidence_2015 = ReportedIncidence(find(Years==2015));

%plot(Years_long,0.1*incidence_2015*ones(size(Years_long)))

legend('Location','northwest')

saveas(gcf,['Extrapolated incidence.png'])

display(['The incidence in 2015 is ', num2str(incidence_2015)])
display(['In 2035 assuming baseline is ', num2str(TBIncidence_base(end-1)),', a ', num2str(100*(incidence_2015/TBIncidence_base(end-1)-1)), '% reduction'])
display(['In 2035 assuming SEA is ', num2str(TBIncidence_SEA(end-1)),', a ', num2str(100*(incidence_2015/TBIncidence_SEA(end-1)-1)), '% reduction'])
display(['In 2035, assuming AMR is ', num2str(TBIncidence_AMR(end-1)),', a ', num2str(100*(incidence_2015/TBIncidence_AMR(end-1)-1)), '% reduction'])
%% plot prevalence
figure('units','normalized','outerposition',[0 0 1 1])

title('Prevalence vs time')

hold on
plot(Years_long, TBPrevalence_SEA(1:end-1)','b:','DisplayName','SEA')
% plot(Years_long, TBPrevalence_SEA_decay(1:end-1)','k','DisplayName','SEA2')
plot(Years_long, TBPrevalence_AMR(1:end-1)','ro-','DisplayName','AMR')
plot(Years_long, TBPrevalence_base(1:end-1)','mx--','DisplayName','base')

plot(Years,ReportedPrevalence,'k--', 'LineWidth',5,'DisplayName','reported')
legend('Location','northwest')

saveas(gcf,['Extrapolated prevalence.png'])

%%
figure('units','normalized','outerposition',[0 0 1 1])
title(['LTBI prevalence, compared to Jordan''s estimate'])
% subplot(2,2,3)
load('Error_aj.mat');

hold on 


LTBI_SEA = compute_LTBI(XELTR_SEA)*100;
LTBI_AMR = compute_LTBI(XELTR_AMR)*100;
LTBI_base = compute_LTBI(XELTR_base)*100;
plot(Years_long, LTBI_SEA(1:end-1),'b:','DisplayName','SEA')

plot(Years_long, LTBI_AMR(1:end-1)','ro-','DisplayName','AMR')
plot(Years_long, LTBI_base(1:end-1)','mx--','DisplayName','base')
plot(years_aj,prevalence_aj,'k--', 'DisplayName','Jordan estimate','LineWidth',5)
% errorbar(years_aj, prevalence_aj, prevalence_aj-LL_aj, UL_aj-prevalence_aj,'DisplayName','95% UI')

legend('Location','northwest')
saveas(gcf,['Extrapolated LTBI prevalence.png'])

