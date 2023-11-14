%% plot incidence

% load data and variable names
% load('XX3.mat')
% load('paramcell.mat')
% load('paramnames.mat')

[beta_list, p_list, w_list, v_list, a_list, d_list, n_list, q1_list, q2_list, TP0_list,E0_list, L0_list, T0_list,R0_list] = paramcell{1:14};
paramgrid = expand_grid(paramcell{:}); 

numSims = size(paramgrid,1);

yearGrid = t0:t0+length(ReportedTB);
%% hard code across w

% Sensitivity analysis across w
% 

% find indices corresponding to w

numw = length(w_list);

% initialize avgIncidences

avgIncidences = zeros(numw,length(XX{1,5}));

parameterIndex = 3; % w

figure


for j=1:numw
    numIndices = numSims/length(w_list);

    

    wnow = w_list(j);
    indices= row_index(paramgrid(:,parameterIndex), wnow);





    subplot(2,2,j)
    plot(yearGrid(1:end-1),ReportedTB, 'LineWidth',3)
    hold on
    title(['w=',num2str(wnow)])
   
    % also store the average    
    avgIncidence = zeros(size(XX{1,5})); % store average

    for k=1:numIndices
        plot( yearGrid, XX{indices(k),5} ,'LineWidth',0.25)
        avgIncidence = avgIncidence + XX{indices(k),5};
    end

    avgIncidence = avgIncidence./numIndices;
    
    avgIncidences(j,:) = avgIncidence;
end
subplot(2,2,4)


    plot(yearGrid(1:end-1), ReportedTB, 'LineWidth',3)
    hold on 
    title('averages')
    plot(yearGrid,avgIncidences)
    legend('Reported Incidence')

saveas(gcf,'Sensitivtyw.png')

%% choose a parameter to be fixed.

% paramIndex = 14; % L0

for paramIndex = 1:14 % number of parameters
    
    param_list = paramcell{paramIndex};
    
    % find indices corresponding to w
    numparams = length(param_list);
    
    if numparams >1

    % initialize avgIncidences
    avgIncidences = zeros(numparams,length(XX{1,5}));
    
    
    figure
    
    
    for j=1:numparams
        paramnamenow = paramnames{paramIndex};

        numIndices = numSims/length(param_list);
    
        
    
            paramnow = param_list(j);
            indices= row_index(paramgrid(:,paramIndex), paramnow);
        
        
        
        
        
            subplot(2,2,j)
            plot(yearGrid(1:end-1),ReportedTB, 'LineWidth',3)
            hold on
            title([paramnames{paramIndex},'=',num2str(paramnow)])
           
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
        
        %figure 
        subplot(2,2,4) %TRACE HARD CODED
        
            plot(yearGrid(1:end-1), ReportedTB, 'LineWidth',3, 'DisplayName','Reported Incidence')
            hold on 
            title('averages')
            for k=1:numparams
                plot(yearGrid,avgIncidences(k,:), 'DisplayName',[paramnamenow,'=',num2str(param_list(k))])
            end
            legend('Location','NorthWest')
        
        saveas(gcf,['Sensitivty_', paramnames{paramIndex}, '.png'])
    end
end