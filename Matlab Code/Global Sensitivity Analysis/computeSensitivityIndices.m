%% computeSensitivityIndices
function [S , ST,  f0]=computeSensitivityIndices(YA,YB,YC_ab,YC_ba)

% extract input
numSims = size(YA,1);
numOutput = size(YA,2);
numInput = size(YC_ab,3);
% experiment k
% output j
% input i

% initialize outputs
S = zeros(numInput, numOutput);
ST = S;

% YC = zeros(numExp, numOutput, numInput);


f0 = zeros(1,numOutput);
for o=1:numOutput
    f0(o) = (1/numSims)*sum( YA(:,o) );
    eqn21 = (1/numSims)*sum( YA(:,o).*YB(:,o) );

    for i=1:numInput
        YACi = YA(:,o).*YC_ab(:,o,i);
    
        % first order sensitivity indices
        numerator = sum(YACi)/(numSims-1) - f0(o)^2;
        denominator = sum( YA(:,o).*YA(:,o) )/numSims - f0(o)^2;
        
        % S(i,o) = numerator/denominator;

        denominator = var(YA(:,o));



   


        % eqn (12)
        Uj = 1/(numSims-1)*sum(YA(:,o).*YC_ba(:,o,i));


        % eqn (18b)
        Unj = 1/(numSims-1)*sum( YA(:,o).*YC_ab(:,o,i) );
        % ST(i,j) = 1-(Unj-f0(j)^2)/denominator;

        % eqn (21) and (22), for E^2(Y)
        E2_ab = (1/numSims)*(sum(YA(:,o).*YB(:,o))); %21
        E2_aa = (1/numSims)*(sum(YA(:,o).*YA(:,o))); %22
        % 
        E2_f0 = f0(o)^2;
        

        
     allY = [YA(:,o); YB(:,o)];
        E2_montecarlo = mean(allY)^2;
        var_montecarlo = mean(allY.^2) - E2_montecarlo;

        denominator = var_montecarlo;
        S(i,o) = (Uj-eqn21)/denominator;

        % total-effect indices
        YBCi = YB(:,o).*YC_ab(:,o,i);
        numerator2 = sum(YBCi)/(numSims) - f0(o)^2;
        

        
        S(i,o) = (Uj-E2_ab)/var_montecarlo;
        ST(i,o) = 1-( Unj - E2_aa )/var_montecarlo;
        % ST(i,o) = 1-numerator2/denominator;

    end
end

