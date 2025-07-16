kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient41'};
allMRFocs = {};
%%
cd('/Users/zjohnson/Documents/MATLAB/KranionExports/');
for i = 1:length(kranionExports)
    cd(kranionExports{i})
    xmlFile = 'model.xml'; % Specify the path to your XML file
    xmlDoc = xmlread(xmlFile);
    
    % Find all <Attribute> elements
    attributes = xmlDoc.getElementsByTagName('ThermometryPhase');
    
    % Initialize an array to hold the extracted values
    attributeData = [];
    
    % Loop through each <Attribute> element
    for j = 0:attributes.getLength-1
        % Get the current <Attribute> element
        attributeNode = attributes.item(j);
        
        % Get the key and type attributes (optional, for verification)
        key = char(attributeNode.getAttribute('key'));
        type = char(attributeNode.getAttribute('type'));
        
        % Find the child nodes for <x>, <y>, <z>
        xValue = str2double(attributeNode.getElementsByTagName('float').item(0).getTextContent);
        yValue = str2double(attributeNode.getElementsByTagName('float').item(1).getTextContent);
        zValue = str2double(attributeNode.getElementsByTagName('float').item(2).getTextContent);
        
        % Store the values in the array
        attributeData = [attributeData; xValue, yValue, zValue];
    end
    
    
    clear mrFocsInHAS;
    KRXpath = '/Users/zjohnson/Documents/MATLAB/KranionExports/';
    [mrFocsInHAS] = mrToHAS([kranionExports{i},'.krx'], KRXpath, attributeData);
    allMRFocs{i} = mrFocsInHAS;
    save('mrFocsInHAS.mat','mrFocsInHAS');
    cd('/Users/zjohnson/Documents/MATLAB/KranionExports/');
end 

%%
fig_seven = figure;
hold on;
axis equal;


text(0, 3.2, 'A', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold'); % Anterior
text(0, -3.4, 'P', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold'); % Posterior
text(-3.5, 0, 'R', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold'); % Left
text(3.4, 0, 'L', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold'); % Right


% Blue unit circle
fadedColor = [0.6 0.6 0.6];

% Cross axes
plot([-3 3], [0 0], 'Color', fadedColor);
plot([0 0], [-3 3], 'Color', fadedColor);

% Concentric circles
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'Color', fadedColor);        % Unit circle
plot(2*cos(theta), 2*sin(theta), 'Color', fadedColor);    % Radius 2 circle
plot(3*cos(theta), 3*sin(theta), 'Color', fadedColor); 

xlabel('$\Delta x$ (mm)');
ylabel('$\Delta y$ (mm)');
axis([-4 4 -4 4]);
xticks(-4:1:4);
yticks(-4:1:4);
grid on;

% Change directory to data folder
baseDir = '/Users/zjohnson/Documents/MATLAB/KranionExports/';
cd(baseDir);

% Store all points from all folders
allPlottedPoints = [];

% Loop through each folder
for i = 1:length(kranionExports)
    disp(kranionExports{i});
    cd(kranionExports{i});

    for j = 1:length(allMRFocs{i})
        clear FocCTijk
        load(['SonicationParameters', num2str(j), '.mat'], 'FocCTijk');
        load(['SonicationDataJul8', num2str(j), '.mat'], 'temps');
        load('Modl.mat');

        % Find max temperature index
        [~, indexOfMax] = max(temps(:));
        [x, y, z, t] = ind2sub(size(temps), indexOfMax);
        [xdim, ydim, ~, ~] = size(temps);

        % Compute relative coordinates for plotting
        xPlot = ceil(xdim/2) - (FocCTijk(1) - allMRFocs{i}{j}(1)) - x;
        yPlot = ceil(ydim/2) - (FocCTijk(2) - allMRFocs{i}{j}(2)) - y;

        % Store the point
        allPlottedPoints = [allPlottedPoints; xPlot, yPlot];
    end

    cd(baseDir);
end

% Round and count repeated points across all folders
roundedPoints = round(allPlottedPoints * 1000) / 1000;
[uniquePts, ~, ic] = unique(roundedPoints, 'rows');
counts = accumarray(ic, 1);

% Avoid overlapping labels
labelPositions = [];

% Plot each unique point and its count label
for k = 1:size(uniquePts, 1)
    xPt = uniquePts(k, 1);
    yPt = uniquePts(k, 2);

    % Plot dot ONCE
    plot(xPt, yPt, 'o', 'MarkerSize', 14, 'LineWidth', 1.2, 'Color', 'b');

    % Try multiple offsets to avoid overlap
    offsets = [ 0.2  0.2;
               -0.2  0.2;
                0.2 -0.2;
               -0.2 -0.2;
                0.3  0;
                0    0.3;
               -0.3  0;
                0   -0.3];

    for o = 1:size(offsets, 1)
        xText = xPt + offsets(o, 1);
        yText = yPt + offsets(o, 2);

        % Check if offset already used
        if isempty(labelPositions) || ~any(all(abs(labelPositions - [xText yText]) < 0.1, 2))
            break;
        end
    end

    % Record label position
    labelPositions = [labelPositions; xText yText];

    % Add count label
    text(xText, yText, sprintf('%d', counts(k)), ...
        'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end

% Finalize figure settings
setFigureProperties(fig_seven, 'spatialError');
xticks(-4:1:4);
yticks(-4:1:4);
grid on;

