function pi_t = pi_extra(tvec)
%{
INPUT: t, R^numt, vector of integers, each t >= 2006
OUTPUT: pi_t, the number of immigrants.  t=<=2022 uses reported data, while
t>2022 extrapolates by fitting to 2006 to 2019
%}

% load('Reported_piCW_20012020');
% Years_20012022 = 2001:2022;
numt = length(tvec);
load('MoreImmigration.mat')
% regression_coeff = polyfit(Years20062022(1:end-3),ReportedImmigration20062022(1:end-3),1); %truncate 2020,2021,2022
regression_coeff = [6.034505494504717e+03  -1.187219945054788e+07];
for i=1:numt
    t=tvec(i);
    if t<2006
        error('t is less than 2006')
    elseif t<= 2019 %t needs to be an integer % ignore 2020 onward
        pi_t(i) = ReportedImmigration20062022(find(Years20062022==t));
    else
        pi_t(i)= regression_coeff(1)*t+regression_coeff(2);
    end

end
