function NumImmigrants = pi_LinearRegression(t)


% load('MoreImmigration.mat')
% coeffs = polyfit(Years20062022,ReportedImmigration20062022,1);

% coeffs = [8.623187533440683e+03 -1.708071594970642e+07]; % fit 2006-2022
coeffs = [6.034505494504717e+03 -1.187219945054788e+07]; %fit 2006-2019

    NumImmigrants = coeffs(1)*t + coeffs(2);

end