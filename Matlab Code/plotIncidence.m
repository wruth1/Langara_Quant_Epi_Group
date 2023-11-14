%% compareIncidence 

function plotIncidence(tgrid, EstimatedIncidence, ReportedIncidence )
  
    
    plot(tgrid, EstimatedIncidence, 'DisplayName','Estimated Incidence')
    hold on
    plot(tgrid(1:end-1), ReportedIncidence, 'r--', 'DisplayName','Reported Incidence')
    grid on;
    title('Incidence vs Time');
    xlabel('time'); ylabel('recovered');
    %legend('show','Location','south');
    legend('EstimatedIncidence','ReportedIncidence', 'Location','south')
end 
