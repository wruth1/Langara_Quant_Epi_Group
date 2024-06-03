%%

set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultContourLineWidth',1.5)
set(groot,'DefaultFunctionContourLineWidth',2)
set(groot,'defaultAxesFontSize',22)
set(groot,'defaultLegendFontSize',22)


%% plot incidence

% load data and variable names
% load('XX3.mat')
load('paramcell47.mat')
load('XX47.mat')
load('paramnames.mat')

imageName = 'Run47';

[beta_list, p_list, w_list, v_list, a_list, d_list, n_list, s_list, q1_list, q2_list, q3_list,X0_list,E0_list, L0_list, T0_list,R0_list] = paramcell{1:16};
paramgrid = expand_grid(paramcell{:});

numSims = size(paramgrid,1);

% t0 = Years(1);
% yearGrid = t0:1:t0+length(ReportedIncidence);
yearGrid = [Years, 2021];

%% find out which parameters we did sensitivity analysis on
close all;

SA_indices = []; % sensitivity analysis indices

for k=1:length(paramcell)
    if length(paramcell{k})>1
        SA_indices = [SA_indices, k];
    end
end

numSA = length(SA_indices);

cell_SA = cell(numSA,1);



%% Multi Trajectory Plot - Incidence



% figure('units','normalized','outerposition',[0 0 1 1])

plot_counter = 0;

for i = 1:numSA % number of parameters

    paramIndex = SA_indices(i);

    param_list = paramcell{paramIndex};
    numparams = length(param_list);




    figure('units','normalized','outerposition',[0 0 1 1])

    % initialize avgIncidences
    avgIncidences = zeros(numparams,length(XX{1,5}));


    % subplot(2,2,plot_counter)


    for j=1:numparams

        paramnamenow = paramnames{paramIndex};

        numIndices = numSims/length(param_list);


        % numparams = length(param_list{j});


        paramnow = param_list(j);
        indices= row_index(paramgrid(:,paramIndex), paramnow);




        if j<=3
            subplot(2,2,j)
            plot(Years,ReportedIncidence, 'LineWidth',3)
            hold on
            title(['Incidence for fixed ',paramnames{paramIndex},'=',num2str(paramnow)])
        end
        % also store the average
        avgIncidence = zeros(size(XX{1,5})); % store average

        for k=1:numIndices
            plot( yearGrid, XX{indices(k),5} ,'LineWidth',0.25)
            avgIncidence = avgIncidence + XX{indices(k),5};
        end

        avgIncidence = avgIncidence./numIndices;

        avgIncidences(j,:) = avgIncidence;

    end
    hold off

    % save average incidences
    cell_SA{i} = avgIncidences;

    %figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,4)

    plot(Years, ReportedIncidence, 'LineWidth',3, 'DisplayName','Reported')
    hold on
    title('averages')
    for k=1:numparams
        plot(yearGrid,avgIncidences(k,:), 'DisplayName',[paramnamenow,'=',num2str(param_list(k))])
    end
    legend('Location','NorthWest')
    hold off

    saveas(gcf,['Incidence_', paramnames{paramIndex}, imageName, '.png'])


end

% Multiple parameters


figure('units','normalized','outerposition',[0 0 1 1])
for i =1:numSA

    paramIndex = SA_indices(i);
    paramnamenow = paramnames{paramIndex};
    
    param_list = paramcell{paramIndex};
    numparams = length(param_list);

    subplot(3,2,i)

    plot(Years, ReportedIncidence, 'LineWidth',3, 'DisplayName','Reported')
    hold on
    title(['Incidence for varying ', paramnamenow])

    avgIncidences = cell_SA{i};
    for k=1:numparams
        plot(yearGrid,avgIncidences(k,:), 'DisplayName',[paramnamenow,'=',num2str(param_list(k))])
    end
    legend('Location','SouthEast')


end
    saveas(gcf,['AvgIncidence_', imageName, '.png'])

%% Multi for prevalence



% plot_counter = 0;

for i = 1:numSA % number of parameters

    paramIndex = SA_indices(i);

    param_list = paramcell{paramIndex};
    numparams = length(param_list);




    % figure('units','normalized','outerposition',[0 0 1 1])

    % initialize avgPrevalences
    avgPrevalences = zeros(numparams,length(XX{1,5}));


    % subplot(2,2,plot_counter)

figure('units','normalized','outerposition',[0 0 1 1])

    for j=1:numparams

        paramnamenow = paramnames{paramIndex};

        numIndices = numSims/length(param_list);


        % numparams = length(param_list{j});


        paramnow = param_list(j);
        indices= row_index(paramgrid(:,paramIndex), paramnow);




        if j<=3
            subplot(2,2,j)
            plot(Years,ReportedPrevalence, 'LineWidth',3)
            hold on
            title(['Prevalence for fixed ',paramnames{paramIndex},'=',num2str(paramnow)])
        end
        % also store the average



        avgPrevalence = zeros(size(XX{1,5}))'; % store average

        for k=1:numIndices
            XELTRi = XX{indices(k),4};
            Ti = XELTRi(:,4)';

            plot( yearGrid, Ti ,'LineWidth',0.25)
            avgPrevalence = avgPrevalence + Ti;
        end

        avgPrevalence = avgPrevalence/numIndices;

        avgPrevalences(j,:) = avgPrevalence;

    end
    hold off

    % save average incidences
    cell_SA{i} = avgPrevalences;

    %figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,4)

    plot(Years, ReportedPrevalence, 'LineWidth',3, 'DisplayName','Reported')
    hold on
    title('averages')
    for k=1:numparams
        plot(yearGrid,avgPrevalences(k,:), 'DisplayName',[paramnamenow,'=',num2str(param_list(k))])
    end
    legend('Location','NorthWest')
    hold off

    saveas(gcf,['Prevalence_', paramnames{paramIndex}, imageName, '.png'])


end

% Multiple parameters


figure('units','normalized','outerposition',[0 0 1 1])
for i =1:numSA

    paramIndex = SA_indices(i);
    paramnamenow = paramnames{paramIndex};
    
    param_list = paramcell{paramIndex};
    numparams = length(param_list);

    subplot(3,2,i)

    plot(Years, ReportedPrevalence, 'LineWidth',3, 'DisplayName','Reported')
    hold on
    title(['Prevalence for varying ', paramnamenow])

    avgPrevalences = cell_SA{i};
    for k=1:numparams
        plot(yearGrid,avgPrevalences(k,:), 'DisplayName',[paramnamenow,'=',num2str(param_list(k))])
    end
    legend('Location','SouthEast')


end
    saveas(gcf,['AvgPrevalence_', imageName, '.png'])


