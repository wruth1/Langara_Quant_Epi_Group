% import data.csv, everything except first row (exclude column headers)
% Define the filename for the GIF
filename = 'tb_prevalence_and_immigration.gif';

% Load coastlines and borders once to reuse
% Save shp file and edit for local file path in shaperead function
load coastlines;
borders = shaperead('C:\Users\timl9\OneDrive\Documents\MATLAB\110m_cultural\ne_110m_admin_0_countries.shp', 'UseGeoCoords', true);

% Set the figure size (width x height in pixels)
figWidth = 1200;
figHeight = 900;

% Loop over each year from 2001 to 2020
for year = 2001:2020
    % Calculate year intervals
    yearInterval = 6+(floor((userInput-2000)/5));
    year = userInput-1991;
    totalImmigrantsThatYear = data{1, year};
    
    % Initialize the world map
    figure('Visible', 'off', 'Position', [100, 100, figWidth, figHeight]); % Set the figure size
    worldmap('World');
    geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'black');
    plabel('off');
    mlabel('off');
    
    % Plot data for each country
    prevalence = lookuptable.TOTALPREV;
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
    
    % Add country borders
    geoshow(borders, 'DisplayType', 'polygon', 'FaceColor', 'none', 'EdgeColor', 'black');
    
    % Immigration size legend
    textm (15, -175, 'Migrant Number Size Reference', 'FontWeight', 'Bold', 'FontSize',  9 );
    plotm(0, -148, 'Marker', 'o', 'MarkerEdgeColor','black', 'MarkerSize', (0.2*sqrt(40000)));
    textm(0,-130, "40,000");
    plotm(-17, -150, 'Marker', 'o', 'MarkerEdgeColor','black', 'MarkerSize', (0.2*sqrt(20000)));
    textm(-17, -132, "20,000");
    plotm(-30, -154, 'Marker', 'o', 'MarkerEdgeColor', 'black', 'MarkerSize', (0.2*sqrt(10000)));
    textm(-30, -135, "10,000");


    % Prevalence legend
    custom_colormap = [linspace(1, 1, 256); linspace(1, 0, 256); linspace(1, 0, 256)]';
    colormap(custom_colormap);
    colorbar('YTick', [0, 0.2, 0.4, 0.6, 0.8, 1], 'YTickLabel', {'0', '12', '24', '36', '48', '60'});
    annotation('textbox', [0.78, 0.1, 0.2, 0.05], 'String', 'LTBI Prevalence (%)', ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    
    % Capture the plot as an image
    frame = getframe(gcf);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if year == 2001
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.6);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.6);
    end
    
    close(gcf); % Close the figure
end

disp('GIF created successfully!');
