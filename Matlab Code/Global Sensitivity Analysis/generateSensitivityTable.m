%generateSensitivityTable

% load("SensitivityIndices2.mat");

[S,ST,f0] = computeSensitivityIndices(YA,YB,YC);


%% df 
% YC_other = cell(numExp,numInput);
% 
% termination_flag = zeros(numExp,numInput);
% 
% for j=1:numInput
% 
% 
%     for i=1:numExp
%         YCotherij = YC_other{i,j};
%         termination_flag(i,j)=YCotherij{2};
%     end
% end
% 
% 
% % number of experiments
% length(termination_flag(:))
% 
% % number that failed to terminate
% length(find(termination_flag==0))



%%

% matrix2latex(S,'Sname')

  rowLabels = {'beta', 'p','w','v','a','d','n','sigma','X0','E0','L0','T0','R0'};
  columnLabels = {'Incidence','Prevalence'};

my_matrix2latex(S, 'outS.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

my_matrix2latex(ST, 'outST.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 

%%

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
    % 
    % my_matrix2latex(RANGES, 'out-ranges.tex', 'rowLabels', rowLabels, 'columnLabels', {'LB','UB'}, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny'); 
