%% load data

% load('Experiment_Output.mat')
% load('Experiment_Parameters.mat');
% load('')

NAME = 'Run 2';
%%


NumSims=size(Experiment_Output,1);



allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);
allE0 = zeros(NumSims,1);
allL0 = zeros(NumSims,1);
allR0 = zeros(NumSims,1);

for i=1:NumSims
    currentx = Experiment_Output{i,2};
    allq1(i)=currentx(1);
    allq2(i) = currentx(2);
    allE0(i) = currentx(3);
    allL0(i) = currentx(4);
    allR0(i) = currentx(5);
end

numbins=20;
%%
% figure
subplot(2,3,1)

histogram(allq1,numbins)
title('Distribution of optimal q1')

subplot(2,3,2)
histogram(allq2,numbins)
title('Distribution of optimal q2')

subplot(2,3,3)
histogram(allE0,numbins)
title('Distribution of optimal E0')

subplot(2,3,4)
histogram(allL0,numbins)
title('Distribution of optimal L0')

subplot(2,3,5)
histogram(allR0,numbins)
% histogram(allR0,numbins)
title('Distribution of optimal R0')
saveas(gcf,['Optimized parameter distribution across sensitivity analysis ', NAME, '.png'])

%% 
