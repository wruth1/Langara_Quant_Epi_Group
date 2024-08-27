%before running, import data.csv, import everything except the first row (exclude column titles)  
    
    userInput =input('Enter a year from 2001 to 2020:\n');
    yearInterval = 6+(floor((userInput-2000)/5));
    year = userInput-1991;
    totalImmigrantsThatYear = data{1, year};
    worldmap('World');
    load coastlines;
    geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black');
    plabel('off');
    mlabel('off');
    prevalencePercentage=(data{:, 4})./(data{:,5});
    for i=1:size(data,1)
        if (data{i,yearInterval}>0)
             propByCountry=(data{i, yearInterval}) / sum(data{:, yearInterval});
             immigrantsByCountry=propByCountry*totalImmigrantsThatYear;
             prevalence_normalized = (prevalencePercentage(i)) / (max(prevalencePercentage));
             rgbValues = [1, 1-prevalence_normalized, 1-prevalence_normalized];
             plotm(data{i, 30}, data{i,31}, 'Marker', 'o', 'MarkerFaceColor', rgbValues, 'MarkerEdgeColor', 'black', 'MarkerSize',(0.2*sqrt(immigrantsByCountry)));
              title(['Tuberculosis Prevalence in Origin Countries of Arriving Canadian Immigrants: ', sprintf('%d', userInput)],'FontSize', 14);
     
         
        end
    
    end

    %adding country borders
    %save the shp file to your device and edit for your local path
    borders = shaperead('C:\Users\timl9\OneDrive\Documents\MATLAB\110m_cultural\ne_110m_admin_0_countries.shp', 'UseGeoCoords', true);
    geoshow(borders, 'DisplayType', 'polygon', 'FaceColor', 'none', 'EdgeColor', 'black');


%Immigration size Legend
textm (15, -175, 'Migrant Number Size Reference', 'FontWeight', 'Bold', 'FontSize',  9 );
plotm(0, -148, 'Marker', 'o', 'MarkerEdgeColor','black', 'MarkerSize', (0.2*sqrt(40000)));
textm(0,-130, "40,000");
plotm(-17, -150, 'Marker', 'o', 'MarkerEdgeColor','black', 'MarkerSize', (0.2*sqrt(20000)));
textm(-17, -132, "20,000");
plotm(-30, -154, 'Marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerSize', (0.2*sqrt(10000)));
textm(-30, -135, "10,000");




%Prevalence Legend
custom_colormap = [linspace(1, 1, 256); linspace(1, 0, 256); linspace(1, 0, 256)]';
% Create colorbar
colorbar;

% Apply custom colormap
colormap(custom_colormap);



% Set tick marks and labels on the colorbar
colorbar_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1];
colorbar_tick_labels = {'0', '12', '24', '36', '48', '60'};
colorbar('YTick', colorbar_ticks, 'YTickLabel', colorbar_tick_labels);

annotation('textbox', [0.78, 0.1, 0.2, 0.05], 'String', 'LTBI Prevalence (%)', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');



