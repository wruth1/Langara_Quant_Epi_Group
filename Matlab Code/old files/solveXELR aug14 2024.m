%% import NationalHobenPlus

M = readtable("NationalHoubenPlus");
% delete Canada
M(25,:) = [];
numc = 167; %number of countries


Names_c= M{:,1};

WHOregion_c = table2array(M(:,"WHOREGION"));

Population_c = table2array(M(:,"POPULATION"));
LTBI_c = table2array(M(:,"ALLLTBI___"))/100;
% LTBI_high = M(:,"ALLLTBI_High____");
% LTBI_low = M(:,"ALLLTBI_Low____");
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

%% compute XELR among new immigrants, year 2014, using method 2


%index for 2014
i14 = find(Immigration_years==2014);
% RP = 10/100; % relapse proportion

XELplusR = zeros(numc,3); %columns respectively mean X, E, L + R

plocal = 2e-1;
wlocal = 0.5;
vlocal = 1.7e-3;
sigmalocal =      7.844347348586576e-05; % average of Run 93's 4 good indices
NgRelapseFraction = 110/1120;


for c=1:numc
    
    XELplusR(c,1) = ((1-LTBI_c(c)))*pi_c(c,i14);
    
    XELplusR(c,2) = INFany_c(c)*pi_c(c,i14);

    XELplusR(c,3) = (LTBI_c(c)-INFany_c(c))*pi_c(c,i14);
end

%
pi2014 = sum(pi_c(:,i14));

qvec2 = zeros(1,4);

qvec2(1) = sum(XELplusR(:,1))/pi2014;
qvec2(2) = sum(XELplusR(:,2))/pi2014;

% 2 by 2 system
M = zeros(2);
b = zeros(2,1);

% qL + qR is known
M(1,:) = [1 1];
b(1) = sum(XELplusR(:,3))/pi2014;

M(2,:) = [-NgRelapseFraction*vlocal (1-NgRelapseFraction)*sigmalocal];
b(2) = NgRelapseFraction*plocal*wlocal*qvec2(2);

sollocal = M\b;
qvec2(3:4) = sollocal;