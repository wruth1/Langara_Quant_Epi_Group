%before running:
%import data.csv to your matlab environment, import everything except the first row (exclude column titles)
%save the countries.shp file to your device, edit the following line to the filepath for the file on your device
    borders = shaperead('C:\Users\timl9\OneDrive\Documents\MATLAB\110m_cultural\ne_110m_admin_0_countries.shp', 'UseGeoCoords', true);

    prevalence2014=(data{:, 4})./(data{:,5});
    year =input('Enter a year from 2001 to 2020:\n');
    yearInterval = 6+(floor((year-2001)/5));
    yearColumn = userInput-1991;
    totalImmigrantsThatYear = data{1, yearColumn};
    prevalence=prevalenceThatYear(year, prevalence2014);

    worldmap('World');
    load coastlines;
    geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black');
    plabel('off');
    mlabel('off');
    for i=1:size(data,1)
        if (data{i,yearInterval}>0)
             propByCountry=(data{i, yearInterval}) / sum(data{:, yearInterval});
             immigrantsByCountry=propByCountry*totalImmigrantsThatYear;
             % normalize with the the highest prevalence in 2000
             prevalence_normalized = (prevalence(i)) / (max(prevalenceThatYear(2001, prevalence2014)));
             rgbValues = [1, 1-prevalence_normalized, 1-prevalence_normalized];
             plotm(data{i, 30}, data{i,31}, 'Marker', 'o', 'MarkerFaceColor', rgbValues, 'MarkerEdgeColor', 'black', 'MarkerSize',(0.2*sqrt(immigrantsByCountry)));
              title(['Tuberculosis Prevalence in Origin Countries of Arriving Canadian Immigrants: ', sprintf('%d', year)],'FontSize', 14);
              
         
        end
    
    end

    %adding country borders
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
colorbar_tick_labels = {'0', '14', '28', '42', '56', '70'};
colorbar('YTick', colorbar_ticks, 'YTickLabel', colorbar_tick_labels);
annotation('textbox', [0.78, 0.1, 0.2, 0.05], 'String', 'LTBI Prevalence (%)', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');

% returns the list of estimated country prevalences in a given year, accounting for decay
function prevalence=prevalenceThatYear(year, prevalence2014)
    prevalence = prevalence2014*(0.99)^(year-2014);
end
