%% load data

% addpath('data and results')
load('Experiment_Output.mat','Experiment_Output')
load('Experiment_Parameters.mat');
%%

% extract error

NumSims=size(Experiment_Output,1);

errorvec = zeros(NumSims,1);

for k=1:NumSims
    localerror = Experiment_Output{k,6};
    errorvec(k) = localerror;
end

errorvecsorted=sort(errorvec);
plot(errorvecsorted(1:NumSims-50))

%% 


