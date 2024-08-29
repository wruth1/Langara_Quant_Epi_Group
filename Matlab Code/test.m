
numInput = 8;
numExp=size(YA,1);

Sblock = zeros(numInput, numExp);
STblock = Sblock;

for s = 1:numExp % s for each sim
    [Slocal, STlocal, f0local] = computeSensitivityIndices(YA(end-s+1:end,:),YB(end-s+1:end,:),YC_ab(end-s+1:end,:,:),YC_ba(end-s+1:end,:,:));
    Sblock(:,s) = Slocal(:,1);
    STblock(:,s) = STlocal(:,1);
end


%%


set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultContourLineWidth',1.5)
set(groot,'DefaultFunctionContourLineWidth',2)
set(groot,'defaultAxesFontSize',32)
set(groot,'defaultLegendFontSize',32)

figure('units','normalized','outerposition',[0 0 1 1])

start = 500;
plot(start:numExp,Sblock(4,start:end), 'DisplayName','v')
hold on
plot(start:numExp,Sblock(2,start:end), '--', 'DisplayName','p')
legend
title('First-order Sensitivity index vs number of Monte Carlo experiments')
xlabel('Number of Monte Carlo experiments')
ylabel('Sensitivity index')

saveas(gcf,['Sensitivity index vs number MC .png'])
