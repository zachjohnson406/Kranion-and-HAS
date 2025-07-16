%%
kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient39','Patient41'};
cValues = [2660	3440 4200];
%%
for kk = 1:length(kranionExports)
    KRXfile = kranionExports{kk};
    KRXpath = '/Users/zjohnson/Documents/MATLAB/KranionExports/';
    FUSF_HAS_Prep_loop(KRXfile,KRXpath);
end
% for i = 1:length(kraniponExports)
%     KRXfile = kranionExports{i};
%     KRXpath = '/Users/zjohnson/Documents/MATLAB/KranionExports/';
%     [~, ~, ~] = FUSF_HAS_Prep_interate_to_converge(KRXfile,KRXpath);
% end 
%% plot models
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports) 
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('Modl.mat');
    load('SonicationParameters1.mat');
    load('maxTempsSimulationMay29badModl.mat','pHAS');
    new_Modl = cleanHASModel(Vhas,FocCTijk(3));
    %new_Modl = prepHAS_Model(new_Modl,0,0,FocCTijk,Sonication.XTilt);
    Modl = new_Modl;
    [xdim, ydim, zdim] = size(Modl);
    centerX = ceil(xdim/2);
    centerY = ceil(ydim/2);
    centerZ = FocCTijk(3);
    
    figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.6]);
    
    % XY slice (Z center)
    subplot(2, 3, 1);
    imagesc(Modl(:, :, centerZ));
    axis image;
    axis off;
    hold on;
    plot(centerY, centerX, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    title('XY Slice');
    
    % XZ slice (Y center)
    subplot(2, 3, 2);
    imagesc(squeeze(Modl(:, centerY, :)));
    axis image;
    axis off;
    hold on;
    plot(centerZ, centerX, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    title('XZ Slice');
    
    % YZ slice (X center)
    subplot(2, 3, 3);
    imagesc(squeeze(Modl(centerX, :, :)));
    axis image;
    axis off;
    hold on;
    plot(centerZ, centerY, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    title('YZ Slice');
    
    
    Modl = Vhas;%prepHAS_Model(Vhas,0,0,FocCTijk,Sonication.XTilt);
    subplot(2, 3, 4);
    imagesc(Modl(:, :, centerZ));
    axis image;
    axis off;
    hold on;
    plot(centerY, centerX, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    
    % XZ slice (Y center)
    subplot(2, 3, 5);
    imagesc(squeeze(Modl(:, centerY, :)));
    axis image;
    axis off;
    hold on;
    plot(centerZ, centerX, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    
    % YZ slice (X center)
    subplot(2, 3, 6);
    imagesc(squeeze(Modl(centerX, :, :)));
    axis image;
    axis off;
    hold on;
    plot(centerZ, centerY, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    cd('..');
end


%%
figure;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports) 
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('maxTempsSimulationJul8.mat');
    load('MRTIdata.mat');
    
    for j = 1:length(maxTempsSimulation)
        if isempty(maxTempsSimulation{j})
            maxTempsSimulation{j} = single(0);
        end
    end

    
    maxSimTemps = cell2mat(maxTempsSimulation);
    maxTa = cell2mat(maxT);
    if strcmp(kranionExports{i}, 'Patient22')
        maxTa = maxTa(1:15);
        maxTa(4) = 44.3712;
    end
    if strcmp(kranionExports{i}, 'Patient17')
        maxTa(6) = 54.5122;
        maxTa(4) = 46.1232;
    end
    if strcmp(kranionExports{i}, 'Patient39')
        maxTa = maxTa(1:5);
    end
    % Convert actual temps to delta above 37C
    deltaTemps = maxTa - 37;

    subplot(rows, cols, i);
    hold on;

    % Scatter plot of data
    plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'red');
    
    % Linear regression
    coeffs = polyfit(maxSimTemps, deltaTemps, 1);  % First-order polynomial fit
    xFit = linspace(min(maxSimTemps), max(maxSimTemps), 100);
    yFit = polyval(coeffs, xFit);
    %plot(xFit, yFit, 'b-', 'LineWidth', 1);  % Plot regression line in blue

    % Plot identity line
    xLimits = [0, 45];
    yLimits = [0, 45];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);  % Identity line
    
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Tempature Rise');
    ylabel('MR Thermometry Tempature Rise')
    cd('..');
end
%%
figure;
hold on;
legendHandles = [];      % To store plot handles
legendLabels = {};       % To store corresponding labels

for kk = 1:length(cValues)
    for i = 3:3
        KRXfile = kranionExports{2};
        cd(KRXfile);
        clear maxT
        clear maxTempsSimulation
        load(['maxTempsSimulationBoneSOS=', num2str(cValues(kk)), '.mat']);
        load('MRTIdata.mat');

        for j = 1:length(maxTempsSimulation)
            if isempty(maxTempsSimulation{j})
                maxTempsSimulation{j} = single(0);
            end
        end

        maxSimTemps = cell2mat(maxTempsSimulation);
        maxTa = cell2mat(maxT);
        

        % Convert actual temps to delta above 37C
        deltaTemps = maxTa - 37;
        

        maxSimTemps(4) = [];
        deltaTemps(4) = [];
        % Scatter plot of data
        h = plot(maxSimTemps, deltaTemps(1:5), 'o', 'MarkerSize', 8);  % Get handle
        legendHandles(end+1) = h;  % Store handle
        legendLabels{end+1} = ['SOS = ' num2str(cValues(kk))];  % Label for legend

        % Plot identity line
        xLimits = [0, 70];
        yLimits = [0, 70];
        plot(xLimits, xLimits, 'k-', 'LineWidth', 1);  % Identity line

        grid on;
        axis image;
        xlim(xLimits);
        ylim(yLimits);
        
        cd('..')
    end
end
% KRXfile = kranionExports{i};
% cd(KRXfile);
% load('maxTempsSimulationCSFPaperValuesGrayMatterAbsorbAtTarget.mat');
% maxSimTemps = cell2mat(maxTempsSimulation);
% h = plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8);
% legendHandles(end+1) = h;  % Store handle
% legendLabels{end+1} = 'Absorption = 0.012 Attenuation = 0.068';

% load('maxTempsSimulationBrainAtten=0.039Absorb=0.012.mat');
% maxSimTemps = cell2mat(maxTempsSimulation);
% h = plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8,'Color','b');
% legendHandles(end+1) = h; 
% legendLabels{end+1} = 'Absorption = 0.012 Attenuation = 0.039';

legend(legendHandles, legendLabels, 'Location', 'northwest','FontSize',16);
hold off;


%% Load data
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
clear all;
SonNum = 8;
cd('patient17')
load(['SonicationData',num2str(SonNum)]);
load(['SonicationParameters',num2str(SonNum)]);
load('Modl.mat');

%%
pout = sum(pressar,4);
Vhas = prepHAS_Model(Vhas,0,0,FocCTijk,0);
Foc_ijk = FocCTijk+1;

[maxforplot, i] = max(temps(:));
[xsim ysim zsim timeForFig] = ind2sub(size(temps), i);
sliceTemps = (temps(:,:,:,timeForFig));

%% HAS Model at focus
Vhas(Vhas == 3) = 2;
Vhas(Vhas == 5) = 3;
Vhas(Vhas == 6) = 4;
Vhas(Vhas == 7) = 5;
%%
figure;
colormap(parula(5));
imagesc(fliplr(flipud(Vhas(:,:,Foc_ijk(3)))));  

title("", 'FontSize', 25);  

ax = gca;  
ax.XAxis.FontSize = 15;  
ax.YAxis.FontSize = 15;  
axis image;

caxis([1 5]);

colorbar;

cb = colorbar;
cb.Ticks = [1, 2, 3, 4, 5];
cb.TickLabels = {'Water', 'Fat', 'Brain', 'Bone', 'Diploe'};
set(cb, 'FontSize', 30);

ylabel('P (mm)', 'FontSize', 25);  
xlabel('L (mm)', 'FontSize', 25);

% Set 4 ticks on the x and y axes
ax.XTick = [50 100 150 200] % 4 ticks evenly spaced across x-axis
ax.YTick = [50 100 150 200]; % 4 ticks evenly spaced across y-axis


%% HAS Model + Pressure at focus
figure;
imgWithPressure(fliplr(flipud((Vhas(:,:,Foc_ijk(3))))), fliplr(((abs(pout(:,:,Foc_ijk(3)))))));

%% HAS Model + Tempature at focus
figure;
imgWithTemp(fliplr(flipud(Vhas(:,:,Foc_ijk(3)))),fliplr(flipud(sliceTemps(:,:,Foc_ijk(3))+37)));
title("Simulation", 'FontSize',30);

%% Read in MR thermometry info
%09-40-25 = patient 30
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000030/20190821093528MR/Analysis/4208-2019-08-21-09-40-25');
%09-25-52 = patient 17
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000017/20190306092122MR/Analysis/4208-2019-03-06-09-25-52');
%10-02-55 = patient 22
cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000022/20190515081333MR/Analysis/4208-2019-05-15-10-02-55');
%15-07-28 = patient 23
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000023/20190605083053MR/Analysis/4208-2019-06-05-15-07-28');
%08-45-50 = patient 27
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000027/20190731100012MR/Analysis/4208-2019-07-31-08-45-50');
%09-21-53 = patient 15
%cd("/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000015/20190213084922MR/Analysis/4208-2019-02-13-09-21-53");
%19-52-40 = patient 24
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000024/20190703104520MR/Analysis/4208-2019-07-03-19-52-40');
%09-51-30 = patient 28
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000028/20190807095047MR/Analysis/4208-2019-08-07-09-51-30');
% patient30
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000030/20190821093528MR/Analysis/4208-2019-08-21-09-40-25');
% patient31
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000031/20190904083315MR/Analysis/4208-2019-09-04-09-17-21');
%patient37
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000037/20191106092257MR/Analysis/4208-2019-11-06-08-23-09');
%patient 39
%cd('/System/Volumes/Data/v/raid10/human_data/IRB00121352_InsightecET/000039/20191120081234MR/Analysis/4208-2019-11-20-08-41-41');
%patient 41
cd('/Users/zjohnson/Documents/MATLAB/KranionExports/4208-2019-12-04-09-13-09')
addpath('/v/raid2/hodeen/ScanData/ReconCode/_MRTIcode')
addpath('/v/raid2/hodeen/ScanData/ReconCode/CntrGUI/')
addpath('/v/raid2/hodeen/ScanData/ReconCode/')
clear tempsAll
clear tempsMRI
clear maxT
clear x
clear y
clear z
clear t
% Grab the folder names
SonicationFolders=grabFiles(pwd,'Sonic','d','d');
% Remove non-folders (here just checking if there's a '.' as in an extension...
for ii=1:length(SonicationFolders)
    if any(strfind(cell2mat(SonicationFolders(ii)),'.'))
        SonicationFolders(ii)=[];
    end
end
% Sort the folder names
for ii=1:length(SonicationFolders)
    SonicationFolders2(ii,1)=str2num(SonicationFolders{ii}(11:end));
end
SonicationFolders2=sort(SonicationFolders2);
for ii=1:length(SonicationFolders2)
    SonicationFolders{ii} = ['Sonication',num2str(SonicationFolders2(ii))];
end; clear ii SonicationFolders2


for sonication=1:length(SonicationFolders)
    cd(SonicationFolders{sonication}) % Go in to the folder
    pause(0.5)

    % Loop over and read in MRTI data
    % Grab all files containing MRTI data and open them
    files=grabFiles(pwd,'5-','raw',[1 2]);
    for ii=1:length(files)
        [tmp] = fopen(cell2mat(files(ii))); [tempsMRI(:,ii)] = fread(tmp,'float');
    end; clear ii tmp
    % Reshape
    tempsMRI = reshape(tempsMRI, [256 256 size(tempsMRI,2)]);

    % Save to structure
    tempsAll{sonication}=tempsMRI;

    clear tempsMRI files

    cd ..
end

for ii=1:numel(tempsAll)
    [x{ii} y{ii} z{ii} t{ii} maxT{ii}] = findHottestVoxel(permute(tempsAll{ii},[1 2 4 3]), [128-3 128+3 128-3 128+3 1 1 1 size(tempsAll{ii},3) ]);
end
asummary=readtable('TreatSummary.csv');

%%
figure;
maxTtable =table2array(asummary(:,25));
for i=1:numel(tempsAll)
    subplot(4,6,i);
    imagesc(tempsAll{i}(:,:,t{i}), [37, 57] );
end


%%
% MRI
figure;
thermometry = permute(tempsAll{SonNum}(:,:,t{SonNum}),[2,1,3]);
ax1 = axes;

thermometry_resized = imresize(thermometry, [size(Vhas, 1), size(Vhas, 2)],'bilinear');
imagesc(ax1, (thermometry) ,[37 57]); axis image, colormap hot;
set(ax1, 'Position',[.07 .11 .685 .815]);
cb = colorbar('Position',[.75 .11 .0675 .815], 'FontSize', 20);
caxis([37 55]);
ylabel(cb, [char(176),'C'], 'FontSize', 25);
ax1.XAxis.FontSize = 15;  % Set font size for x-axis ticks
ax1.YAxis.FontSize = 15;
xlabel(ax1, 'L (mm)', 'FontSize', 25);
ylabel(ax1, 'P (mm)', 'FontSize', 25);
ax1.XTick = [50 100 150 200]; % 4 ticks evenly spaced across x-axis
ax1.YTick = [50 100 150 200];
title("MR thermometry",'FontSize', 30)

%% Compare MR temp to simulated temp for 1 sonication. 
figure;
hold on;
cd('/Users/zjohnson/Documents/MATLAB/KranionExports')
cd('Patient17')
load('MRTIdata.mat')

acqTime = 3.456;
err = (acqTime/2)*ones(1,size(tempsAll{SonNum},3));
MRTItempOneSonication = zeros(1,size(tempsAll{SonNum},3));
acqTimeMean = zeros(1,size(tempsAll{SonNum},3));
for ii = 1:size(tempsAll{SonNum},3)
    columns = x{SonNum} - 1 : x{SonNum} + 1; 
    rows = y{SonNum} - 1 : y{SonNum} + 1;
    subRegion = tempsAll{SonNum}(columns,rows,ii);
    MRTItempOneSonication(ii) = mean(subRegion(:));
    acqTimeMean(ii) = ((ii-1) + 0.5) * (acqTime);
end 

load(['SonicationDataMay20',num2str(SonNum)]);

maxTempsSonicationX = zeros(size(temps,4),1);
for i = 1:length(maxTempsSonicationX)
    slice = temps(:,:,:,i);
    maxTempsSonicationX(i) = max(slice(:));
end 
timevec = timevec-1;

plot([timevec], [maxTempsSonicationX], '-', 'LineWidth', 2);  % Thicker line
hold on
plot([0, acqTimeMean(2:end)-acqTime], [0, MRTItempOneSonication(2:end)-37], 'o-', 'MarkerSize', 10, 'LineWidth', 2);  % Larger markers and thicker lines
ax = gca;
ax.FontSize = 16;
legend('Simulation', 'MRTI', 'FontSize', 24);
xlabel('Time (Seconds)', 'FontSize', 25);
ylabel(['Maximum Temperature', char(176),'C'], 'FontSize', 25);
xlim([0,26])
grid on;
title("Tempature rise over time", 'FontSize',30)
axis square;

                   
