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