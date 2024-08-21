M = readtable("NationalHoubenPlus");
% delete Canada
M(25,:) = [];
numc = 167; %number of countries


Names_c= M{:,1};

WHOregion_c = table2array(M(:,"WHOREGION"));

Population_c = table2array(M(:,"POPULATION"));
% LTBI_c = table2array(M(:,"ALLLTBI___"))/100;
% LTBIhigh_c = M(:,"ALLLTBI_High____");
% LTBIlow_c = M(:,"ALLLTBI_Low____");
% 
% INF1_c = table2array(M(:,"x1stINF_2Years___"))/100;
% INFany_c = table2array(M(:,"AnyINF_2Years___"))/100;
% 


WHOregion_c_numeric = grp2idx(categorical(WHOregion_c));
%% import annual immigration from each country

Annual_Immigration_table = readtable("AnnualImmigration.csv");

Immigration_years = table2array(Annual_Immigration_table(1,2:end));
numt = length(Immigration_years);

pi_c = table2array(Annual_Immigration_table(2:end,2:end));



numw = 6;
%% collate immigration by country into immigration by WHO region
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

save('pi_w20012020.mat','pi_w')

%% load WHO region data 
M2 = readtable("RegionHouben.csv");

% pull prevalence of LTBI mean, high, low
LTBI_w = table2array(M2(:,"LTBIPrevalenceMean"))/100;
LTBIhigh_w = table2array(M2(:,"LTBIPrevalenceHigh"))/100;
LTBIlow_w = table2array(M2(:,"LTBIPrevalenceLow"))/100;

INF_w = table2array(M2(:,"RecentMean"))/100;
INFhigh_w = table2array(M2(:,"RecentHigh"))/100;
INFlow_w  = table2array(M2(:,"RecentLow"))/100;


% WHOregion_c_numeric = grp2idx(categorical(WHOregion_c));



% plot immigration from each WHO geographic region vs time

% build string_w

numw=6;


string_w = {'AFR','AMR', 'EMR','EUR','SEAR','WPR'}; 


% Example data
string_w = {'AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR'};
numbers=LTBI_w';

% Combine string_w with the numbers using cellfun
legend_labels = cellfun(@(s, n) [s ' ' num2str(n) '%'], string_w, num2cell(numbers*100), 'UniformOutput', false);


% legend_labels = strcat(string_w, ' (', string(num2str(numbers(:))), ')');

%% trace KEZIA clean up this figure too
figure('units','normalized','outerposition',[0 0 1 1])
plot(Immigration_years,pi_w)
title('Number of migrants by WHO geographic region vs year')
legend(legend_labels,'Location','NorthWest')

saveas(gcf,['migrants by WHO geographic region vs year.png'])
