% call this to create pi_c, annual immigration by country; and pi_w, by WHO
% geographic region


%% import NationalHobenPlus
addpath('data and results')
M = readtable("NationalHoubenPlus");
% delete Canada
M(25,:) = [];
numc = 167; %number of countries


Names_c= M{:,1};

WHOregion_c = table2array(M(:,"WHOREGION"));
WHOregion_c_numeric = grp2idx(categorical(WHOregion_c)); % hash

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

WHOregion_c_numeric = grp2idx(categorical(WHOregion_c));

%% collate immigration by country into immigration by WHO region
numw = 6;
pi_w = zeros(numw, numt);


region1_indices = find(WHOregion_c_numeric==1);
region2_indices = find(WHOregion_c_numeric==2);
region3_indices = find(WHOregion_c_numeric==3);
region4_indices = find(WHOregion_c_numeric==4);
region5_indices = find(WHOregion_c_numeric==5);
region6_indices = find(WHOregion_c_numeric==6);

for t=1:numt
    pi_w(1,t) = sum(pi_c(region1_indices,t));
    pi_w(2,t) = sum(pi_c(region2_indices,t));
    pi_w(3,t) = sum(pi_c(region3_indices,t));
    pi_w(4,t) = sum(pi_c(region4_indices,t));
    pi_w(5,t) = sum(pi_c(region5_indices,t));
    pi_w(6,t) = sum(pi_c(region6_indices,t));
end
% TBI_c = zeros(numc,1);
% for c=1:numc
%     TBI_c(c) = WHOincidence(WHOregion_c{c});
% end

%%  save data

Years_20012020 = Immigration_years;
piC_20012020 = pi_c;
piW_20012020 = pi_w;

save('Reported_piCW_20012020', 'Years_20012020','piC_20012020','piW_20012020')
