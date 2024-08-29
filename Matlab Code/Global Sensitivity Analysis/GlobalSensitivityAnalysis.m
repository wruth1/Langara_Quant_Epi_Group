%% 

addpath('Global Sensitivity Analysis\')
addpath('functions')

%% 
load('pi_extrapolated.mat')
load('qvec.mat')
load('GOODPARAM.mat')


bioParameters_start = GOODPARAM(1:8);
bioParameters_start(8) = GOODOUTPUT(4);
IC_start = GOODOUTPUT(5:9);

Params_start = [bioParameters_start];
scale_dP = 10/100;

%% load qvec

load('piW_base.mat')
qvec2_t = qvecW_base;

%% Initialize data

numExp = 500000;

numInput = length(Params_start);
numOutput = 2;


% A and B are matrices of data
A = zeros(numExp, numInput);
B = A;

for k=1:numExp
    A(k,:) = createRandomParameters(Params_start,scale_dP);
    B(k,:) = createRandomParameters(Params_start,scale_dP);
end

%% load initial data
% TFP0 = 6186950;
% 
% X0 =    2.435279313142126e+06;
% E0 = 1.484646245171584e+04;
% 
% L0 = 3.096493365566524e+06;
% % T0 = 1.070931823865891e+03;
% T0 = 1081;
% R0 =      6.392498588396340e+05;
% TFP0 = X0 + E0+L0+T0+R0;
% 
% IC = [X0, E0, L0, T0, R0];

%% 
% save outputs
YA = zeros(numExp, numOutput);
YB = YA;
% YA_other = cell(numExp,1);
% YB_other = cell(numExp,1);

for n=1:numExp
    params_n = A(n,:);
    bioParameters = params_n(1:8);
    % IC = params_n(9:13);
    IC = IC_start;
    [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2_t);
   
    
    YA(n,:) = [TBIncidence_long(end-1),TBPrevalence_long(end-1)];

end

% YB
for n=1:numExp
    params_n = B(n,:);

    bioParameters = params_n(1:8);
    % IC = params_n(9:13);
    IC = IC_start;

    [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2_t);
   
    
    YB(n,:) = [TBIncidence_long(end-1),TBPrevalence_long(end-1)];

end

%% now compute YC
YC_ab = zeros(numExp, numOutput, numInput);
YC_ba = zeros(numExp, numOutput, numInput);
% YC_other = cell(numExp,numInput);

for i=1:numInput
    % C is B, but with ith column replaced by A
    C_ba = B; % mostly b
    C_ba(:,i) = A(:,i);

    % instead, C is A, but with ith column replaced by B
    C_ab = A; % mostly A
    C_ab(:,i) = B(:,i);

    for n=1:numExp
        params_n = C_ab(n,:);

        bioParameters = params_n(1:8);
        IC = IC_start;

        [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2_t); 
        YC_ab(n,:,i) = [TBIncidence_long(end-1),TBPrevalence_long(end-1)];
    
        % now for 
        params_n = C_ba(n,:);

        bioParameters = params_n(1:8);
        IC = IC_start;
        [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2_t); 
        YC_ba(n,:,i) = [TBIncidence_long(end-1),TBPrevalence_long(end-1)];
            

    end
end


%% 
% save('SensitivityIndices.mat')

save SensitivityIndices3.mat A B C_ba C_ab YA YB YC_ba YC_ab;
% 
% YA1 = YA;
% YB1 = YB;
% YC_ba1 = YC_ba;
% YC_ab1 = YC_ab;
%% 
[S,ST,f0] = computeSensitivityIndices(YA,YB,YC_ab,YC_ba);


  rowLabels = {'beta', 'p','w','v','a','d','n','sigma','X0','E0','L0','T0','R0'};
  columnLabels = {'Incidence','Prevalence'};

my_matrix2latex(S, 'outS.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

my_matrix2latex(ST, 'outST.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.5f', 'size', 'tiny'); 


%% see if S ST is stabilizing as numSims goes to infinity

Sblock = zeros(numInput, numExp);
STblock = Sblock;

for s = 1:numExp % s for each sim
    [Slocal, STlocal, f0local] = computeSensitivityIndices(YA(end-s+1:end,:),YB(end-s+1:end,:),YC_ab(end-s+1:end,:),YC_ba(end-s+1:end,:));
    Sblock(:,s) = Slocal(:,1);
    STblock(:,s) = STlocal(:,1);
end
%%
figure
start = 2500;
plot(Sblock(4,start:end))