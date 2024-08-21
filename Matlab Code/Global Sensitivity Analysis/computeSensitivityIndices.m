%% computeSensitivityIndices
function [S , ST,  f0]=computeSensitivityIndices(YA,YB,YC)

% extract input
numSims = size(YA,1);
numOutput = size(YA,2);
numInput = size(YC,3);
% experiment k
% output j
% input i

% initialize outputs
S = zeros(numInput, numOutput);
ST = S;

% YC = zeros(numExp, numOutput, numInput);


f0 = zeros(1,numOutput);
for j=1:numOutput
    f0(j) = (1/numSims)*sum( YA(:,j) );

    for i=1:numInput
        YACi = YA(:,j).*YC(:,j,i);
    
        % first order sensitivity indices
        numerator = sum(YACi)/numSims - f0(j)^2;
        denominator = sum( YA(:,j).*YA(:,j) )/numSims - f0(j)^2;
        S(i,j) = numerator/denominator;

        % total-effect indices
        YBCi = YB(:,j).*YC(:,j,i);
        numerator2 = sum(YBCi)/numSims - f0(j)^2;
        ST(i,j) = 1-numerator2/denominator;
    end
end

