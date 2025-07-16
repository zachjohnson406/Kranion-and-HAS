%% Figure 1 (Modeling flow diagram - just the model)

SonNum = 8;
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
cd('Patient17');
fullfilepath = dir;
load(['SonicationDataMay23SwappedDiploeNewModl',num2str(SonNum)]);
load(['SonicationParameters',num2str(SonNum)]);
load('Modl.mat');
pout = sum(pressar,4);
Vhas = prepHAS_Model(Vhas,0,0,FocCTijk,0);
Foc_ijk = FocCTijk+1;
[maxforplot, i] = max(temps(:));
[xsim ysim zsim timeForFig] = ind2sub(size(temps), i);
sliceTemps = (temps(:,:,:,timeForFig));
Vhas(Vhas == 3) = 2;Vhas(Vhas == 5) = 3;Vhas(Vhas == 6) = 4;Vhas(Vhas == 7) = 5;

fig_one = figure;
colormap(parula(5));

% Get image dimensions
[ny, nx, nz] = size(Vhas);
slice_data = fliplr(flipud(Vhas(:,:,92)));

% Define spatial resolution (adjust these values based on your actual voxel size)
% Assuming 1mm per pixel - modify these if your resolution is different
dx = 1; % mm per pixel in x direction
dy = 1; % mm per pixel in y direction

% Create coordinate vectors centered at origin
x_coords = (1:nx) * dx - (nx/2 + 0.5) * dx;  % Center at 0
y_coords = (1:ny) * dy - (ny/2 + 0.5) * dy;  % Center at 0

% Display image with proper coordinates
imagesc(x_coords, y_coords, slice_data);
ax = gca;
axis image;

% Set colorbar
cb = colorbar;
cb.Ticks = [1, 2, 3, 4, 5];
cb.TickLabels = {'Water', 'Fat', 'Brain', 'Bone', 'Diploe'};

% Labels
ylabel('P/A (mm)');
xlabel('R/L (mm)');

% Set symmetric ticks around center
max_extent_x = max(abs(x_coords));
max_extent_y = max(abs(y_coords));

% Create symmetric tick marks (adjust spacing as needed)
tick_spacing = 50; % mm
x_ticks = [-tick_spacing*2, -tick_spacing, 0, tick_spacing, tick_spacing*2];
y_ticks = [-tick_spacing*2, -tick_spacing, 0, tick_spacing, tick_spacing*2];

% Filter ticks to only show those within the image bounds
x_ticks = x_ticks(x_ticks >= min(x_coords) & x_ticks <= max(x_coords));
y_ticks = y_ticks(y_ticks >= min(y_coords) & y_ticks <= max(y_coords));

ax.XTick = x_ticks;
ax.YTick = y_ticks;


setFigureProperties(fig_one,'modlFigure');
%% Model with pressure
fig_two = figure; % Added semicolon
ax1 = axes(fig_two); % Create axes handle, not figure handle
imgWithPressure(fliplr(flipud((Vhas(:,:,92)))), fliplr(((abs(pout(:,:,Foc_ijk(3)))))), '', ax1);
setFigureProperties(fig_two,'modlWithPressureFigure');
%% Model with tempature
fig_three = figure;
imgWithTemp(fliplr(flipud(Vhas(:,:,92))),fliplr(flipud(sliceTemps(:,:,Foc_ijk(3)))));
setFigureProperties(fig_three,'modlWithTempFigure');
%% MR thermometry
[maxT,tempsAll,t,x,y] = getThermometryInfo('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000017/20190306092122MR/Analysis/4208-2019-03-06-09-25-52',17);
%%
fig_four = figure;
thermometry = permute(tempsAll{SonNum}(:,:,t{SonNum}),[2,1,3]); %so
ax1 = axes;
thermometry_resized = imresize(thermometry, [size(Vhas, 1), size(Vhas, 2)],'bilinear');
imagesc(ax1, (thermometry-37) ,[0 30]); axis image, colormap hot;
set(ax1, 'Position',[.07 .11 .685 .815]);
cb = colorbar('Position',[.75 .11 .0675 .815]);
clim([0, 30]);
ylabel(cb, '$^\circ$C', 'Interpreter', 'latex');
xlabel(ax1, 'L (mm)');
ylabel(ax1, 'P (mm)');
ax1.XTick = [50 100 150 200];
ax1.YTick = [50 100 150 200];
title("MR thermometry")
fname = 'modlWithPressureFig';
setFigureProperties(fig_four,'MRThermometry');
%% Tempature over time
kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient39','Patient41'};
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
cd(kranionExports{1});
load('MRTIdata.mat')
fig_five = figure;
SonNum=7;
load(['SonicationDataJul8',num2str(SonNum)]);
plotTempOverTime(temps,tempsAll,timevec,SonNum,x,y);
setFigureProperties(fig_five,'tempOverTime');


%% Measured vs. simulated temps
fig_six = figure;
maxTempFile = 'maxTempsSimulationJul8.mat';
tempOverTimeFile = 'SonicationDataJul8';
plotMeasuredVsSimulated(maxTempFile,tempOverTimeFile,3,3);
%%
picturewidth = 27; % set this parameter and keep it forever
hw_ratio = 1.2; % feel free to play with this ratio
%set(findall(fig,'-property','FontSize'),'FontSize',17); % adjust fontsize to your document
set(findall(fig_six,'-property','Box'),'Box','off'); % optional
set(findall(fig_six,'-property','Interpreter'),'Interpreter','latex');
set(findall(fig_six,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_six,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
pos = get(fig_six,'Position');
set(fig_six,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
%%
print(fig_six,'measuredVsSimulatedTempatureRise','-dpdf','-painters');
%%

%%
fig_three = figure;
tile = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% Define voxel size (adjust as needed)
voxel_size = [1, 1]; % [dx, dy] in mm

% First tile
ax1 = nexttile(tile,1);
subPlotWithTemp(fliplr(flipud(Vhas(:,:,92))), fliplr(flipud(sliceTemps(:,:,Foc_ijk(3)))), '', ax1, voxel_size);

% Second tile
ax2 = nexttile(tile,2);
[maxT, tempsAll, t, x, y] = getThermometryInfo('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000017/20190306092122MR/Analysis/4208-2019-03-06-09-25-52',17);
thermometry = permute(tempsAll{SonNum}(:,:,t{SonNum}), [2,1,3]);
thermometry_resized = imresize(thermometry, [size(Vhas,1), size(Vhas,2)], 'bilinear');

% Create coordinate vectors for second tile (same as first tile)
[ny, nx] = size(thermometry_resized);
dx = voxel_size(1); % mm per pixel in x direction
dy = voxel_size(2); % mm per pixel in y direction
x_coords = (1:nx) * dx - (nx/2 + 0.5) * dx;  % Center at 0
y_coords = (1:ny) * dy - (ny/2 + 0.5) * dy;  % Center at 0

% Plot with proper coordinates
imagesc(ax2, x_coords, y_coords, (thermometry_resized - 37), [0 30]);
axis(ax2, 'image'); % Changed from 'square' to 'image' for consistency
colormap(ax2, 'hot');

% Add colorbar
cb = colorbar(ax2);
cb.Label.String = '$^\circ$C';
cb.Label.Interpreter = 'latex';
cb.Label.VerticalAlignment = 'middle';
cb.Label.Position = [0.5 30+0.75];
cb.Label.Rotation = 0;
cb.FontSize = 13;

% Set axis labels
xlabel(ax2, 'R/L (mm)','FontSize',13);
ylabel(ax2, 'P/A (mm)','FontSize',13);
title(ax2, 'MR thermometry','FontSize',17);

% Create symmetric tick marks (same as first tile)
tick_spacing = 50; % mm
x_ticks = [-tick_spacing*2, -tick_spacing, 0, tick_spacing, tick_spacing*2];
y_ticks = [-tick_spacing*2, -tick_spacing, 0, tick_spacing, tick_spacing*2];

% Filter ticks to only show those within the image bounds
x_ticks = x_ticks(x_ticks >= min(x_coords) & x_ticks <= max(x_coords));
y_ticks = y_ticks(y_ticks >= min(y_coords) & y_ticks <= max(y_coords));

ax2.XTick = x_ticks;
ax2.YTick = y_ticks;

% Third tile
ax3 = nexttile(tile,3);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
cd('Patient17')
load('MRTIdata.mat')
load(['SonicationDataJul8',num2str(SonNum)]);
plotTempOverTime(temps,tempsAll,timevec,SonNum,x,y);


%%
picturewidth = 27; % set this parameter and keep it forever
hw_ratio = 0.4; % feel free to play with this ratio
%set(findall(fig,'-property','FontSize'),'FontSize',17); % adjust fontsize to your document
set(findall(fig_three,'-property','Box'),'Box','off'); % optional
set(findall(fig_three,'-property','Interpreter'),'Interpreter','latex');
set(findall(fig_three,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_three,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
pos = get(fig_three,'Position');
set(fig_three,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
%%
cd('/Users/zjohnson/Documents/MATLAB/KranionFigures');
fname = 'oneSonicationMRTIvsSim';
print(fig_three,fname,'-dpdf','-painters','-fillpage');
%%
fig_eight = figure;
tile = tiledlayout(2, 2);%, 'Padding', 'compact', 'TileSpacing', 'compact');

[ny, nx, nz] = size(Vhas);
dx = 1; % mm per pixel in x direction
dy = 1; % mm per pixel in y direction
x_coords = (1:nx) * dx - (nx/2 + 0.5) * dx; % Center at 0
y_coords = (1:ny) * dy - (ny/2 + 0.5) * dy; % Center at 0

ax1 = nexttile(tile, 1);
% Read PNG image
png_image = imread('kranionScreenOne.png');
% Display PNG with same coordinate system
image(ax1, png_image);
axis off 
axis image


title(ax1, 'a');

ax2 = nexttile(tile, 2);
slice_data = fliplr(flipud(Vhas(:,:,92)));
% Display image with proper coordinates
imagesc(ax2, x_coords, y_coords, slice_data);
colormap(ax2, parula(5));
axis(ax2, 'image');
% Set colorbar
cb = colorbar(ax2);
cb.Ticks = [1, 2, 3, 4, 5];
cb.TickLabels = {'Water', 'Fat', 'Brain', 'Bone', 'Diploe'};
% Labels and ticks
ylabel(ax2, 'P/A (mm)');
xlabel(ax2, 'R/L (mm)');
ax2.XTick = x_ticks;
ax2.YTick = y_ticks;
title(ax2, 'b');

ax3 = nexttile(tile, 3);
imgWithPressure(fliplr(flipud((Vhas(:,:,92)))), fliplr(((abs(pout(:,:,Foc_ijk(3)))))), '', ax3);
title(ax3, 'c');

ax4 = nexttile(tile, 4);
% Read PNG image
png_image2 = imread('kranionScreenTwp.png');

imagesc(ax4,png_image2);
axis image
axis off 
title(ax4, 'd');

%%
picturewidth = 34; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(fig_eight,'-property','FontSize'),'FontSize',17); % adjust fontsize to your document
set(findall(fig_eight,'-property','Box'),'Box','off'); % optional
set(findall(fig_eight,'-property','Interpreter'),'Interpreter','latex');
set(findall(fig_eight,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
set(fig_eight,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
pos = get(fig_eight,'Position');
set(fig_eight,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
%% transducer fig 
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
cd('Patient17');
load('ERFA8.mat');
%%
fig_nine = figure;
hold on;
for j = 1:7
    scatter3(Seg(j).ElemLoc_Cart(:,1), Seg(j).ElemLoc_Cart(:,2), Seg(j).Center_Cart(:,3),100,'black');
    % Calculate the center of the current segment's dots
    center_x = mean(Seg(j).ElemLoc_Cart(:,1));
    center_y = mean(Seg(j).ElemLoc_Cart(:,2));
    center_z = mean(Seg(j).Center_Cart(:,3));
    % Add text label at the center with background
    text(center_x, center_y, center_z, num2str(j), 'BackgroundColor', [0.9 0.9 0.9], 'FontSize',25,'Color', 'black', 'FontWeight', 'bold', 'EdgeColor', 'black');
end
axis image;
grid on;
ylabel('P/A','FontSize',17);
xlabel('R/L','FontSize',17);
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
%%
setFigureProperties(fig_nine,'transducer')