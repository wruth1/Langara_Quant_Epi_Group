function RandomParams = createRandomParameters(Params,scale_dP)
%{
    input: 
    Params, in R^Nx1

    output: R^(Nx2) [Params *(1-scale_dP), Params*(1+scale_dP)]
%}
    
    numParams = length(Params);
    
    Range = zeros(numParams,2);
    
    Range(:,1) = Params *(1-scale_dP);
    Range(:,2) = Params*(1+scale_dP);

    %  % 1 beta = bioParameters(1); %TB infectivity
    % RANGES =  [0.00000001 * 0.9, 0.00000001 * 1.1; % beta
    % 4/100/0.5, 10/100/0.5; % p
    % 0.5, 0.5; %w
    % 0.00224 0.00487; %v
    % 0.0522 0.06; %a tb mortality
    % 0.6 0.8; % d active TB recovery
    % 0.0071-0.001 0.0071+0.001; % natural death rate
    % 8e-5 3e-3; % sigma relapse rate
    % 0.563/100, 0.920/100; %qE initial condition
    % 5/100 25/100; % qL
    % 5/100 25/100; % qR
    % 2.435279313142126e+06 2.435279313142126e+06; %X0    
    % 1.484646245171584e+04 1.484646245171584e+04; %E0
    % 3.096493365566524e+06 3.096493365566524e+06; % L0
    % 1081 1081; %T0
    % 6.392498588396340e+05 6.392498588396340e+05; %R0
    % 1.1, 1+1/3; % prevalence amplification initial
    % 8 8]; % power relapse


    numparams = size(Range,1);

    randvec = rand(numparams,1);


    RandomParams = Range(:,1)+randvec.*(Range(:,2)-Range(:,1));

    RandomParams = RandomParams';
end