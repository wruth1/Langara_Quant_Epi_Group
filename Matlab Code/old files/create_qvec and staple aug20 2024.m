

%% import NationalHobenPlus
addpath('data and results')
M = readtable("NationalHoubenPlus");
% delete Canada
M(25,:) = [];
numc = 167; %number of countries


Names_c= M{:,1};

WHOregion_c = table2array(M(:,"WHOREGION"));

Population_c = table2array(M(:,"POPULATION"));
LTBI_c = table2array(M(:,"ALLLTBI___"))/100;
LTBIhigh_c = M(:,"ALLLTBI_High____");
LTBIlow_c = M(:,"ALLLTBI_Low____");

INF1_c = table2array(M(:,"x1stINF_2Years___"))/100;
INFany_c = table2array(M(:,"AnyINF_2Years___"))/100;

% import incidence

IncidenceTable = readtable("TBIncidenceNational.csv");
IncidenceTable(25,:) = []; % delete Canada
%
TBI_c = table2array(IncidenceTable(:,2));

%% import annual immigration from each country

Annual_Immigration_table = readtable("AnnualImmigration.csv");

Immigration_years = table2array(Annual_Immigration_table(1,2:end));
numt = length(Immigration_years);

pi_c = table2array(Annual_Immigration_table(2:end,2:end));



% TBI_c = zeros(numc,1);
% for c=1:numc
%     TBI_c(c) = WHOincidence(WHOregion_c{c});
% end



%% compute XELR among new immigrants, year 2014, using 46% method


%index for 2014
i14 = find(Immigration_years==2014);
RP = 10/100; % relapse proportion

XELRN_c = zeros(numc,4); %columns respectively mean X, E, L, R

plocal = 2e-1;
wlocal = 0.5;
vlocal = 1.7e-3;
sigmalocal =      7.844347348586576e-05; % average of Run 93's 4 good indices
NgRelapseFraction = 110/1120;

qLqR_ratio = 0.46;


qvec_t = zeros(numt,4); 
for t = 1:numt
    qvec = zeros(1,4);
    for c=1:numc
        
        % solve Ax = b
        A = zeros(4);
        b = zeros(4,1);
    
        % E + L + R = 22.4% * pi_c(2014)
        A(1,:) = [0 1 1 1];
        b(1) = LTBI_c(c)*pi_c(c,t);
    
        % E/(X+E+L+R) = 1.31%
        A(2,:) = [0 1 0 0 ] - [1 1 1 1]*INFany_c(c);
        b(2) = 0;
    
        % X/(X+k*(L+R)) = 1.24%/1.31%
        % K = 0.50; % reinfection protection
        % A(3,:) = [1 0 0 0] - [1 0 K K]*INF1_c(c)/INFany_c(c);
        % b(3) = 0;
    
        %(1+ RP/1-RP)pwE + (1+RP/1-RP)vL = pi_c / 100k * TB Incidence
        % A(3,2) = plocal*wlocal*(1+ RP/(1-RP));
        % A(3,3) = vlocal*(1+ RP/(1-RP));
        % % TBI_local = WHOincidence(WHOregion_c{c});
        % b(3) =  TBI_c(c)*pi_c(c,i_t)/100000;
     
        % L/R = 46%
        A(3,:) = [0 0 1 -qLqR_ratio];
        b(3)=0;
    
        % sigmaR / ( pwE + vL + sigma R ) = 110/1120;
        % A(3,:) = [0 0 0 sigmalocal] - [0 plocal*wlocal vlocal sigmalocal]*NgRelapseFraction;
        % b(3) = 0;
    
        % X+E+L+R = pi_c
        A(4,:) = [1 1 1 1];
        b(4) = pi_c(c,t);
    
        XELR_local = A\b;
    
        XELRN_c(c,:) = XELR_local';
    end
    qvec(1) = sum(XELRN_c(:,1))/sum(pi_c(:,t));
    qvec(2) = sum(XELRN_c(:,2))/sum(pi_c(:,t));
    qvec(3) = sum(XELRN_c(:,3))/sum(pi_c(:,t));
    qvec(4) = sum(XELRN_c(:,4))/sum(pi_c(:,t));

    qvec_t(t,:) = qvec;
end

% in particular, the average proportions are
qvecbar = mean(qvec_t);

save('qvec.mat','qvec_t','qvecbar')

%% 

figure('units','normalized','outerposition',[0 0 1 1])

plot(Immigration_years, qvec_t)
title('q_X, q_E, q_L, q_R vs time')
legend('q_X','q_E','q_L','q_R')

%% 

%% load WHO region data 
M2 = readtable("RegionHouben.csv");

% pull prevalence of LTBI mean, high, low
LTBI_w = table2array(M2(:,"LTBIPrevalenceMean"))/100;
LTBIhigh_w = table2array(M2(:,"LTBIPrevalenceHigh"))/100;
LTBIlow_w = table2array(M2(:,"LTBIPrevalenceLow"))/100;

INF_w = table2array(M2(:,"RecentMean"))/100;
INFhigh_w = table2array(M2(:,"RecentHigh"))/100;
INFlow_w  = table2array(M2(:,"RecentLow"))/100;


WHOregion_c_numeric = grp2idx(categorical(WHOregion_c));


numw = 6;


%% find q  for WHO regions

% ImmigrantsW_Time();

load('pi_w20012020.mat')

qvecW_t = zeros(numt,4); 
for t = 1:numt
    qvec2 = zeros(1,4);

    

    for w=1:numw
        
        % solve Ax = b
        A = zeros(4);
        b = zeros(4,1);
    
        % E + L + R = 22.4% * pi_w(2014)
        A(1,:) = [0 1 1 1];
        b(1) = LTBI_w(w)*pi_w(w,t);
    
        % E/(X+E+L+R) = 1.31%
        A(2,:) = [0 1 0 0 ] - [1 1 1 1]*INF_w(w);
        b(2) = 0;
    
        % X/(X+k*(L+R)) = 1.24%/1.31%
        % K = 0.50; % reinfection protection
        % A(3,:) = [1 0 0 0] - [1 0 K K]*INF1_c(c)/INFany_c(c);
        % b(3) = 0;
    
        %(1+ RP/1-RP)pwE + (1+RP/1-RP)vL = pi_c / 100k * TB Incidence
        % A(3,2) = plocal*wlocal*(1+ RP/(1-RP));
        % A(3,3) = vlocal*(1+ RP/(1-RP));
        % % TBI_local = WHOincidence(WHOregion_c{c});
        % b(3) =  TBI_c(c)*pi_c(c,i_t)/100000;
     
        % L/R = 46%
        A(3,:) = [0 0 1 -qLqR_ratio];
        b(3)=0;
    
        % sigmaR / ( pwE + vL + sigma R ) = 110/1120;
        % A(3,:) = [0 0 0 sigmalocal] - [0 plocal*wlocal vlocal sigmalocal]*NgRelapseFraction;
        % b(3) = 0;
    
        % X+E+L+R = pi_c
        A(4,:) = [1 1 1 1];
        b(4) = pi_w(w,t);
    
        XELR_local = A\b;
    
        XELRN_w(w,:) = XELR_local';
    end
    qvec2(1) = sum(XELRN_w(:,1))/sum(pi_w(:,t));
    qvec2(2) = sum(XELRN_w(:,2))/sum(pi_w(:,t));
    qvec2(3) = sum(XELRN_w(:,3))/sum(pi_w(:,t));
    qvec2(4) = sum(XELRN_w(:,4))/sum(pi_w(:,t));

    qvecW_t(t,:) = qvec2;
end

% in particular, the average proportions are
% qvec2bar = mean(qvec2_t);


% Years20062020 = 2006:2020;
qvec_20062020 = qvecW_t(6:end,:);
% siez

%% make qvec2_t20062035 by stapling

numC = 4; % number of compartments. qX qE qL qR

Years_long = 2006:2035;
numt_long = length(Years_long);
numt_short = size(qvec_20062020,1);



qvec_baseline = zeros(numt_long-numt_short,numC);
% assume only immigrants from SEA region, highest prevalence of LTBI
qvec_SEA = zeros(numt_long-numt_short,numC);
% assume only immigrants from AMR region, highest prevalence of LTBI
qvec_AMR = zeros(numt_long-numt_short,numC);

numt_extra = numt_long-numt_short;


% load('pi_extrapolated.mat');
% piW_SEA





Years_extra = Years_long(end-numt_extra+1:end);
piW_SEA = zeros(numt_extra,numw);
piW_AMR = zeros(numt_extra,numw);

% SEA is index 5
piW_SEA(:,5) = pi_LinearRegression(Years_extra);
piW_AMR(:,2) = pi_LinearRegression(Years_extra);



%qvecW_extra = zeros(numt_extra,4); 

qvecW_SEA = zeros(numt_extra,4); 

for t = 1:numt_extra
    qvec_local = zeros(1,4);

    

    for w=1:numw
        
        % solve Ax = b
        A = zeros(4);
        b = zeros(4,1);
    
        % E + L + R = 22.4% * pi_w(2014)
        A(1,:) = [0 1 1 1];
        b(1) = LTBI_w(w)*piW_SEA(t,w);
    
        % E/(X+E+L+R) = 1.31%
        A(2,:) = [0 1 0 0 ] - [1 1 1 1]*INF_w(w);
        b(2) = 0;
  
        % L/R = 46%
        A(3,:) = [0 0 1 -qLqR_ratio];
        b(3)=0;
    
        % sigmaR / ( pwE + vL + sigma R ) = 110/1120;
        % A(3,:) = [0 0 0 sigmalocal] - [0 plocal*wlocal vlocal sigmalocal]*NgRelapseFraction;
        % b(3) = 0;
    
        % X+E+L+R = pi_c
        A(4,:) = [1 1 1 1];
        b(4) = piW_SEA(t,w);
    
        XELR_local = A\b;
    
        XELRN_w(w,:) = XELR_local';
    end
    qvec_local(1) = sum(XELRN_w(:,1))/sum(piW_SEA(t,:));
    qvec_local(2) = sum(XELRN_w(:,2))/sum(piW_SEA(t,:));
    qvec_local(3) = sum(XELRN_w(:,3))/sum(piW_SEA(t,:));
    qvec_local(4) = sum(XELRN_w(:,4))/sum(piW_SEA(t,:));

    qvecW_SEA(t,:) = qvec_local;
end

%% 

%%
qvecW_SEA2 = build_qvecW_t(piW_SEA,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );
qvecW_SEA3 = build_qvecW_t(piW_SEA,LTBI_w,INF_w,  qLqR_ratio  );
qvecW_AMR = build_qvecW_t(piW_AMR,(1-Houben_weight)*LTBI_w+Houben_weight*LTBIhigh_w,INF_w,  qLqR_ratio  );

%
qvec_SEA20062035 = zeros(numt_long, numC);
qvec_SEA20062035(1:numt_short,:) = qvec_20062020;
qvec_SEA20062035(numt_short+1:end,:) = qvecW_SEA2;
qvec3_SEA20062035 = qvec_SEA20062035;
qvec3_SEA20062035(numt_short+1:end,:) = qvecW_SEA3;

qvec_AMR20062035 = zeros(numt_long, numC);
qvec_AMR20062035(1:numt_short,:) = qvec_20062020;
qvec_AMR20062035(numt_short+1:end,:) = qvecW_AMR;


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

Houben_weight = (qLTBI_computed-Houben_mean)/(Houben_high-Houben_mean);

%% 


[XELTR, TBIncidence_SEA, TBPrevalence_SEA] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec_SEA20062035);
[XELTR, TBIncidence_SEA3, TBPrevalence_SEA3] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec3_SEA20062035);
[XELTR, TBIncidence_AMR, TBPrevalence_AMR] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec_AMR20062035);

% load incidence data to plot
load('ReportedTB20062020.mat');


figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
title('Incidence vs time')

hold on
plot(Years_long, TBIncidence_SEA(1:end-1),'b:')
plot(Years_long, TBIncidence_SEA3(1:end-1),'ro-')

% plot(Years_long, TBIncidence_AMR(1:end-1)','ro-')
plot(Years,ReportedIncidence,'k--', 'LineWidth',5)
incidence_2015 = ReportedIncidence(find(Years==2015));
%plot(Years_long,0.1*incidence_2015*ones(size(Years_long)))

legend('SEA','AMR','Reported','Location','NorthWest')
subplot(2,1,2)
title('Prevalence vs time')

hold on
plot(Years_long, TBPrevalence_SEA(1:end-1)','b:')
plot(Years_long, TBPrevalence_AMR(1:end-1)','ro-')
plot(Years,ReportedPrevalence,'k--', 'LineWidth',5)

saveas(gcf,['Incidence Prevalence extrapolated.png'])
