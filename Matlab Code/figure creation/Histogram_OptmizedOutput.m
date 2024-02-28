%% load data

% load('XX.mat')
% load('Experiment_Parameters.mat');
% load('')
%%
imageName = 'Run15_newbounds';

NumSims=size(XX,1);



allq1 = zeros(NumSims,1);
allq2 = zeros(NumSims,1);
allX0 = zeros(NumSims,1);
allE0 = zeros(NumSims,1);
allL0 = zeros(NumSims,1);
allR0 = zeros(NumSims,1);

for i=1:NumSims
    currentx = XX{i,2};
    allq1(i)=currentx(1);
    allq2(i) = currentx(2);
    allX0(i) = currentx(3);
    allE0(i) = currentx(4);
    allL0(i) = currentx(5);
    allR0(i) = currentx(6);
end

numbins=20;
%% Plot Histograms
figure(1)

subplot(2,3,1)

% histfit(allq1,numbins)
histfit(allq1,numbins)
title('Distribution of optimal q1')

subplot(2,3,2)
histfit(allq2,numbins)
title('Distribution of optimal q2')


subplot(2,3,3)
histfit(allE0,numbins)
title('Distribution of optimal X0')

subplot(2,3,4)
histfit(allE0,numbins)
title('Distribution of optimal E0')

subplot(2,3,5)
histfit(allL0,numbins)
title('Distribution of optimal L0')

subplot(2,3,6)
histfit(allR0,numbins)
% histfit(allR0,numbins)
title('Distribution of optimal R0')
saveas(gcf,['Optimized parameter distribution across sensitivity analysis',imageName, '.png'])



%% plot sorted 

OutputNames = {'q1','q2','E0','L0','R0' };

numSims = length(allq1);

figure(2)

% figure
subplot(2,3,1)

semilogy(sort(allq1),'DisplayName','Optimized')
title('Sorted log optimal q1')
hold on
for i=1:length(q1_list)
    semilogy(q1_list(i).*ones(numSims,1),'--','DisplayName',['I.C. ', num2str(i)])
end
% legend('Location','SouthEast')

subplot(2,3,2)
semilogy(sort(allq2),'DisplayName','Optimized')
hold on
for i=1:length(q2_list)
    semilogy(q2_list(i).*ones(numSims,1),'--','DisplayName',['I.C. ', num2str(i)])
end
title('Sorted log optimal q2')
% legend('Location','SouthEast')


subplot(2,3,3)
semilogy(sort(allE0),'DisplayName','Optimized')
hold on
for i=1:length(E0_list)
    semilogy(E0_list(i).*ones(numSims,1),'--','DisplayName',['I.C. ', num2str(i)])
end
title('Sorted log optimal E0')
% legend('Location','SouthEast')

subplot(2,3,4)
semilogy(sort(allL0),'DisplayName','Optimized')
hold on
for i=1:length(L0_list)
    semilogy(L0_list(i).*ones(numSims,1),'--','DisplayName',['I.C. ', num2str(i)])
end
title('Sorted log optimal L0')
% legend('Location','NorthWest')

subplot(2,3,5)
semilogy(sort(allR0),'DisplayName','Optimized')
hold on
for i=1:length(R0_list)
    semilogy(R0_list(i).*ones(numSims,1),'--','DisplayName',['I.C. ', num2str(i)])
end
title('Sorted log optimal R0')
legend('Location','SouthEast')

saveas(gcf,['Optimized parameter log sorted with starting values', imageName, '.png'])


