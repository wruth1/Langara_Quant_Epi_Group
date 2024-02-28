%% load data

% addpath('data and results')
% load('Experiment_Output.mat','Experiment_Output')
% load('Experiment_Parameters.mat');
%%
NumSims=size(Experiment_Output,1);

ScatterData = zeros(NumSims,2);


j = 14;

for k = 1:NumSims
    localparams = Experiment_Output{k,1};

    % extract which index
    ScatterData(k,:) = [localparams(j), Experiment_Output{k,6}];

end

scatter(ScatterData(:,1),ScatterData(:,2));
% title = ('')

saveas(gcf,'ScatterPlotErrorVsParameter.png')