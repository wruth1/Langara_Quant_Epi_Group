%% 
load('pi_extrapolated.mat')
load('qvec.mat')
load('GOODPARAM.mat')


bioParameters_start = GOODPARAM(1:8);
bioParameters_start(8) = GOODOUTPUT(4);
IC_start = GOODOUTPUT(5:9);

Params_start = [bioParameters_start];
scale_dP = 1/100;

%% load qvec

load('piW_base.mat')
qvec2 = qvecW_base;

%% Initialize data

numExp = 100;

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
   
    
    YA(n,:) = [TBIncidence_long(end),TBPrevalence_long(end)];

end

% YB
for n=1:numExp
    params_n = B(n,:);

    bioParameters = params_n(1:8);
    % IC = params_n(9:13);
    IC = IC_start;

    [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2);
   
    
    YB(n,:) = [TBIncidence_long(end),TBPrevalence_long(end)];

end

%% now compute YC
YC = zeros(numExp, numOutput, numInput);
% YC_other = cell(numExp,numInput);

for i=1:numInput
    % C is B, but with jth column replaced by A
    C = B;
    C(:,i) = A(:,i);

    for n=1:numExp
        params_n = C(n,:);

        bioParameters = params_n(1:8);
        % IC = params_n(9:13);
            IC = IC_start;

        [XELTR, TBIncidence_long, TBPrevalence_long] = solveGuoWu4_extrapolate(bioParameters, IC, pi_extrapolated,qvec2);
       
        
        YC(n,:,i) = [TBIncidence_long(end),TBPrevalence_long(end)];


    end
end


%% 
save('SensitivityIndices.mat')

save SensitivityIndices2.mat A B C YA YB YC;
%% 
[S,ST,f0] = computeSensitivityIndices(YA,YB,YC);


  rowLabels = {'beta', 'p','w','v','a','d','n','sigma','X0','E0','L0','T0','R0'};
  columnLabels = {'Incidence','Prevalence'};

my_matrix2latex(S, 'outS.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

my_matrix2latex(ST, 'outST.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

