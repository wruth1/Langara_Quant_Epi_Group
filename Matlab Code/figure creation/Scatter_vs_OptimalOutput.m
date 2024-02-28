%% load data

% addpath('data and results')
% load('Experiment_Output.mat','Experiment_Output')
% load('Experiment_Parameters.mat');
Experiment_Output = XX;

%%
NumSims=size(Experiment_Output,1);

% ScatterData = zeros(NumSims,2);


% j = 14;


% Error = Experiment_Output{1:NumSims,6};

% Scatter Data in R^{NumSims x 2 x numOptimal}.
% 2nd coordinate 2 is for x and y values of scatter plot



numOptimal = length(Experiment_Output{1,2});
ScatterData = zeros(NumSims, 2, numOptimal);

for j=1:numOptimal
    for k = 1:NumSims
        % OptimalOutput = Error(k);
        
        % scatter Error vs Optimal Output
        localY =  Experiment_Output{k,2};
        ScatterData(k,1,j) = localY(j);
        % y is error
        ScatterData(k,2,j) = Experiment_Output{k,6};
    end
end

OutputNames = {'q1','q2','E0','L0','R0' };

for j=1:numOptimal
    subplot(2,3,j)
    scatter(ScatterData(:,1,j),ScatterData(:,2,j));
    title('Error vs Optimized ', OutputNames{j})
end
% title = ('')

saveas(gcf,'ScatterPlotErrorVsOptimalOutput.png')

%% prune outliers

% pull Errors from ScatterData
Error = ScatterData(:,2,1);
% after staring at scatter plot, we prune error above 6
ErrorTOL = 9;

prune_indices = find(Error>ErrorTOL);

prunedScatterData = ScatterData;
prunedScatterData(prune_indices,:,:) = [];

figure(2)
for j=1:numOptimal
    subplot(3,2,j)
    scatter(prunedScatterData(:,1,j),(prunedScatterData(:,2,j)));
    title('Error vs Optimized ', OutputNames{j})
end
saveas(gcf,'ScatterPlotErrorVsOptimalOutput_Cut2.png')

% pull Errors from ScatterData
Error = ScatterData(:,2,1);
% after staring at scatter plot, we prune error above 6
ErrorTOL = 6;

prune_indices = find(Error>ErrorTOL);

prunedScatterData = ScatterData;
prunedScatterData(prune_indices,:,:) = [];

figure(2)
for j=1:numOptimal
    subplot(3,2,j)
    scatter(prunedScatterData(:,1,j),(prunedScatterData(:,2,j)));
    title('Error vs Optimized ', OutputNames{j})
end
saveas(gcf,'ScatterPlotErrorVsOptimalOutput_Cut3.png')
