%% Patients to plot
kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient41'};
sdr = {0.57, 0.39, 0.38,0.68,0.37,0.49,0.49,0.48,0};
%% plot measured and sim tempature at 3rd time point
figure;
targetTimePoint = 3;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationMay20Values.mat')

    targetTime = (acqTime/2) + ((targetTimePoint-2) * 3.456);
    maxTempsTimePoint = {};
    maxTempsTimePointMR = {};
    for jj = 1:length(maxTempsSimulation)

        load(['SonicationDataMay20',num2str(jj),'.mat']);        
        % Find the two surrounding time indices
        timevec = timevec-1;
        idx_before = find(timevec <= targetTime, 1, 'last');
        idx_after = find(timevec >= targetTime, 1, 'first');

        if idx_before == idx_after  % Exact match
            interpolated_temp = temps(:,:,:,idx_before);
        else
            t1 = timevec(idx_before);
            t2 = timevec(idx_after);
            w1 = (t2 - targetTime) / (t2 - t1);
            w2 = (targetTime - t1) / (t2 - t1);

            % Linear interpolation of the temperature volumes
            interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
        end

        % Find the maximum temperature in that interpolated volume
        interpolated_temp = interpolated_temp(:,:,30:end-30);
        max_temp = max(interpolated_temp(:));

        maxTempsTimePoint{jj} = max_temp;
        tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
        maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
    end
    for j = 1:length(maxTempsTimePoint)
        if isempty(maxTempsTimePoint{j})
            maxTempsTimePoint{j} = single(0);
        end
    end
    if(strcmp(KRXfil,'Patient31'))
        maxTempsTimePointMR(5) = [];
        maxTempsTimePoint(5) = [];
    end

    maxSimTemps = cell2mat(maxTempsTimePoint);
    maxTa = cell2mat(maxTempsTimePointMR);
    deltaTemps = maxTa - 37;

    subplot(rows, cols, i);
    hold on;
    plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'red');
    xLimits = [0, 30];
    yLimits = [0, 30];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Tempature Rise');
    ylabel('MR Thermometry Tempature Rise')
    title(KRXfil)
    cd('..');
end
%% plot mrti
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    for ii=1:numel(tempsAll)
        [x{ii} y{ii} z{ii} t{ii} maxT{ii}] = findHottestVoxel(permute(tempsAll{ii},[1 2 4 3]), [128-2 128+2 128-2 128+2 1 1 1 size(tempsAll{ii},3) ]);
    end

    for j=1:numel(tempsAll)
        subplot(4,6,j);
        imagesc(tempsAll{i}(:,:,t{j}), [37, 57] );
    end


    cd('..');
end
%% plot temp over time per sonication with marker for 3rd time point
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('maxTempsSimulationMay23SwappedDiploeValuesNewModl.mat');
    load('MRTIdata.mat');

    for ii = 1:numel(tempsAll)
        [x{ii}, y{ii}, z{ii}, t{ii}, maxT{ii}] = findHottestVoxel(permute(tempsAll{ii},[1 2 4 3]), [128-2 128+2 128-2 128+2 1 1 1 size(tempsAll{ii},3)]);
    end

    figure;
    n = length(maxTempsSimulation) - 1;
    colsSub = ceil(sqrt(n)) + 1;
    rowsSub = ceil(n / colsSub) + 1;

    acqTime = 3.456;
    targetTimePoint = 3 -1;
    targetTime = (acqTime / 2) + ((targetTimePoint - 1) * acqTime);  % same as in second block

    for SonNum = 1:length(maxTempsSimulation)
        err = (acqTime/2)*ones(1,size(tempsAll{SonNum},3));
        MRTItempOneSonication = zeros(1,size(tempsAll{SonNum},3));
        acqTimeMean = zeros(1,size(tempsAll{SonNum},3));
        for ii = 1:size(tempsAll{SonNum},3)
            columns = x{SonNum} - 1 : x{SonNum} + 1;
            rows = y{SonNum} - 1 : y{SonNum} + 1;
            subRegion = tempsAll{SonNum}(columns, rows, ii);
            MRTItempOneSonication(ii) = max(subRegion(:));
            acqTimeMean(ii) = ((ii - 1) + 0.5) * acqTime;
        end

        subplot(rowsSub, colsSub, SonNum);
        hold on;
        load(['SonicationDataMay23SwappedDiploeNewModl', num2str(SonNum)]);

        maxTempsSonicationX = zeros(size(temps,4), 1);
        for j = 1:length(maxTempsSonicationX)
            slice = temps(:,:,:,j);
            maxTempsSonicationX(j) = max(slice(:));
        end

        plot(timevec -1, maxTempsSonicationX, '-', 'LineWidth', 1.5);  % Sim temps
        plot([0, acqTimeMean(2:end) - acqTime], [0, MRTItempOneSonication(2:end) - 37], 'o-', 'MarkerSize', 10);  % MRTI temps

        % Add vertical line for 3rd acquisition time point
        xline(targetTime, 'r--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');

        title(['Sonication ', num2str(SonNum)]);
        xlabel('Time (s)');
        ylabel('Temperature Rise (°C)');
        grid on;
        hold off;
    end

    sgtitle(strrep(kranionExports{i}, '_', '\_'));  % Cleaner figure title
    cd('..');
end
%% look at power output
figure;
hold on;
for i = 1:11
    load(['son',num2str(i),'_var.mat']);
    plot(power_acou);
end
%% plot measured and sim tempature at 3rd time point with sonication markers
figure;
targetTimePoint = 3;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');
for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationMay23SwappedDiploeValuesNewModl.mat')

    targetTime = (acqTime / 2) + ((targetTimePoint-2) * 3.456);
    maxTempsTimePoint = {};
    maxTempsTimePointMR = {};
    for jj = 1:length(maxTempsSimulation)
        load(['SonicationDataMay23SwappedDiploeNewModl', num2str(jj), '.mat']);
        % Find the two surrounding time indices
        idx_before = find(timevec <= targetTime, 1, 'last');
        idx_after = find(timevec >= targetTime, 1, 'first');

        if idx_before == idx_after  % Exact match
            interpolated_temp = temps(:,:,:,idx_before);
        else
            t1 = timevec(idx_before);
            t2 = timevec(idx_after);
            w1 = (t2 - targetTime) / (t2 - t1);
            w2 = (targetTime - t1) / (t2 - t1);

            % Linear interpolation of the temperature volumes
            interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
        end

        % Find the maximum temperature in that interpolated volume
        max_temp = max(interpolated_temp(:));

        maxTempsTimePoint{jj} = max_temp;
        tempsAllSlice = tempsAll{jj}(128 - 2:128 + 2, 128 - 2:128 + 2, targetTimePoint);
        maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
    end
    for j = 1:length(maxTempsTimePoint)
        if isempty(maxTempsTimePoint{j})
            maxTempsTimePoint{j} = single(0);
        end
    end

    maxSimTemps = cell2mat(maxTempsTimePoint);
    maxTa = cell2mat(maxTempsTimePointMR);
    deltaTemps = maxTa - 37;

    subplot(rows, cols, i);
    hold on;
    h = plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'red');
    xLimits = [0, 30];
    yLimits = [0, 30];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Temperature Rise');
    ylabel('MR Thermometry Temperature Rise');
    title(KRXfil);

    % Add sonication number labels to each point
    for jj = 1:length(maxSimTemps)
        text(maxSimTemps(jj) + 0.5, deltaTemps(jj) + 0.5, num2str(jj), 'FontSize', 8, 'Color', 'blue');
    end
    
    cd('..');
end

%% first 3 time points may 20th values
figure;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationMay20Values.mat')

    all_maxSimTemps = [];
    all_deltaTemps = [];

    for targetTimePoint = 1:3
        targetTime = (acqTime/2) + ((targetTimePoint - 2) * acqTime);
        maxTempsTimePoint = {};
        maxTempsTimePointMR = {};

        for jj = 1:length(maxTempsSimulation)
            load(['SonicationDataMay20', num2str(jj), '.mat']);
            if(size(temps,4) == 0)
                maxTempsTimePointMR{jj} = 0;
                maxTempsTimePoint{jj} = 0;
                continue;
            end
            % Find the two surrounding time indices
            idx_before = find(timevec <= targetTime, 1, 'last');
            idx_after = find(timevec >= targetTime, 1, 'first');

            if isempty(idx_before) || isempty(idx_after)
                interpolated_temp = zeros(size(temps(:,:,:,1))); % fallback
            elseif idx_before == idx_after  % Exact match
                interpolated_temp = temps(:,:,:,idx_before);
            else
                t1 = timevec(idx_before);
                t2 = timevec(idx_after);
                w1 = (t2 - targetTime) / (t2 - t1);
                w2 = (targetTime - t1) / (t2 - t1);

                % Linear interpolation of the temperature volumes
                interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
            end

            max_temp = max(interpolated_temp(:));
            maxTempsTimePoint{jj} = max_temp;

            if targetTimePoint <= size(tempsAll{jj},3)
                tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
                maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
            else
                maxTempsTimePointMR{jj} = 37;  % default baseline if index exceeds
            end
        end

        % Clean up empty entries
        for j = 1:length(maxTempsTimePoint)
            val = maxTempsTimePoint{j};
            if isempty(val) || ~isnumeric(val)
                maxTempsTimePoint{j} = 0;
            else
                maxTempsTimePoint{j} = double(val);  % or single(val), if you're using single precision
            end
        end


        maxSimTemps = cell2mat(maxTempsTimePoint);
        maxTa = cell2mat(maxTempsTimePointMR);
        deltaTemps = maxTa - 37;

        all_maxSimTemps = [all_maxSimTemps; maxSimTemps(:)];
        all_deltaTemps = [all_deltaTemps; deltaTemps(:)];
    end

    subplot(rows, cols, i);
    hold on;
    plot(all_maxSimTemps, all_deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'red');
    xLimits = [0, 30];
    yLimits = [0, 30];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Temperature Rise');
    ylabel('MR Thermometry Temperature Rise');
    title(KRXfil, 'Interpreter', 'none');
    cd('..');
end
%% plot measured and sim temp for iterative adjustment
figure;
targetTimePoint = 4;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationIterativeAdjValues.mat');

    targetTime = (acqTime/2) + ((targetTimePoint-2) * 3.456);
    maxTempsTimePoint = {};
    maxTempsTimePointMR = {};
    pHAS_value = NaN;  % initialize

    for jj = 1:length(maxTempsSimulation)
        load(['SonicationDataIterativeAdj', num2str(jj), '.mat']);        
        if jj == 1
            % Save the value of pHAS.a_abs(5) from the first sonication
            if isfield(pHAS, 'a_abs') && length(pHAS.a_abs) >= 5
                pHAS_value = pHAS.a_abs(5);
            end
        end

        timevec = timevec - 1;
        idx_before = find(timevec <= targetTime, 1, 'last');
        idx_after = find(timevec >= targetTime, 1, 'first');

        if idx_before == idx_after
            interpolated_temp = temps(:,:,:,idx_before);
        else
            t1 = timevec(idx_before);
            t2 = timevec(idx_after);
            w1 = (t2 - targetTime) / (t2 - t1);
            w2 = (targetTime - t1) / (t2 - t1);
            interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
        end

        interpolated_temp = interpolated_temp(:,:,30:end-30);
        max_temp = max(interpolated_temp(:));

        maxTempsTimePoint{jj} = max_temp;
        tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
        maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
    end

    for j = 1:length(maxTempsTimePoint)
        if isempty(maxTempsTimePoint{j})
            maxTempsTimePoint{j} = single(0);
        end
    end

    if strcmp(KRXfil, 'Patient31')
        maxTempsTimePointMR(5) = [];
        maxTempsTimePoint(5) = [];
    end

    maxSimTemps = cell2mat(maxTempsTimePoint);
    maxTa = cell2mat(maxTempsTimePointMR);
    deltaTemps = maxTa - 37;

    subplot(rows, cols, i);
    hold on;
    plot(maxSimTemps, deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'blue');
    xLimits = [0, 40];
    yLimits = [0, 40];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Tempature Rise');
    ylabel('MR Thermometry Tempature Rise');
    title(['adjusted absoption: ' , num2str(pHAS.a_abs(5))]);

    cd('..');
end
%% first 3 time points iterative adjustment
figure;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationIterativeAdjValues.mat')

    all_maxSimTemps = [];
    all_deltaTemps = [];

    for targetTimePoint = 1:3
        targetTime = (acqTime/2) + ((targetTimePoint - 2) * acqTime);
        maxTempsTimePoint = {};
        maxTempsTimePointMR = {};

        for jj = 1:length(maxTempsSimulation)
            load(['SonicationDataIterativeAdj', num2str(jj), '.mat']);
            if(size(temps,4) == 0)
                maxTempsTimePointMR{jj} = 0;
                maxTempsTimePoint{jj} = 0;
                continue;
            end
            % Find the two surrounding time indices
            idx_before = find(timevec <= targetTime, 1, 'last');
            idx_after = find(timevec >= targetTime, 1, 'first');

            if isempty(idx_before) || isempty(idx_after)
                interpolated_temp = zeros(size(temps(:,:,:,1))); % fallback
            elseif idx_before == idx_after  % Exact match
                interpolated_temp = temps(:,:,:,idx_before);
            else
                t1 = timevec(idx_before);
                t2 = timevec(idx_after);
                w1 = (t2 - targetTime) / (t2 - t1);
                w2 = (targetTime - t1) / (t2 - t1);

                % Linear interpolation of the temperature volumes
                interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
            end

            max_temp = max(interpolated_temp(:));
            maxTempsTimePoint{jj} = max_temp;

            if targetTimePoint <= size(tempsAll{jj},3)
                tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
                maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
            else
                maxTempsTimePointMR{jj} = 37;  % default baseline if index exceeds
            end
        end

        % Clean up empty entries
        for j = 1:length(maxTempsTimePoint)
            val = maxTempsTimePoint{j};
            if isempty(val) || ~isnumeric(val)
                maxTempsTimePoint{j} = 0;
            else
                maxTempsTimePoint{j} = double(val);  % or single(val), if you're using single precision
            end
        end


        maxSimTemps = cell2mat(maxTempsTimePoint);
        maxTa = cell2mat(maxTempsTimePointMR);
        deltaTemps = maxTa - 37;

        all_maxSimTemps = [all_maxSimTemps; maxSimTemps(:)];
        all_deltaTemps = [all_deltaTemps; deltaTemps(:)];
    end

    subplot(rows, cols, i);
    hold on;
    plot(all_maxSimTemps, all_deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'blue');
    xLimits = [0, 30];
    yLimits = [0, 30];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Temperature Rise');
    ylabel('MR Thermometry Temperature Rise');
    title(['absorption: ', num2str(pHAS.a_abs(5))]);
    cd('..');
end



%% first 3 time points swapped diploe
figure;
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load('maxTempsSimulationMay23SwappedDiploeValuesNewModl.mat')

    all_maxSimTemps = [];
    all_deltaTemps = [];

    for targetTimePoint = 1:3
        targetTime = (acqTime/2) + ((targetTimePoint - 2) * acqTime);
        maxTempsTimePoint = {};
        maxTempsTimePointMR = {};

        for jj = 1:length(maxTempsSimulation)
            load(['SonicationDataMay23SwappedDiploeNewModl', num2str(jj), '.mat']);
            if(size(temps,4) == 0)
                maxTempsTimePointMR{jj} = 0;
                maxTempsTimePoint{jj} = 0;
                continue;
            end
            % Find the two surrounding time indices
            idx_before = find(timevec <= targetTime, 1, 'last');
            idx_after = find(timevec >= targetTime, 1, 'first');

            if isempty(idx_before) || isempty(idx_after)
                interpolated_temp = zeros(size(temps(:,:,:,1))); % fallback
            elseif idx_before == idx_after  % Exact match
                interpolated_temp = temps(:,:,:,idx_before);
            else
                t1 = timevec(idx_before);
                t2 = timevec(idx_after);
                w1 = (t2 - targetTime) / (t2 - t1);
                w2 = (targetTime - t1) / (t2 - t1);

                % Linear interpolation of the temperature volumes
                interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
            end

            max_temp = max(interpolated_temp(:));
            maxTempsTimePoint{jj} = max_temp;

            if targetTimePoint <= size(tempsAll{jj},3)
                tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
                maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
            else
                maxTempsTimePointMR{jj} = 37;  % default baseline if index exceeds
            end
        end

        % Clean up empty entries
        for j = 1:length(maxTempsTimePoint)
            val = maxTempsTimePoint{j};
            if isempty(val) || ~isnumeric(val)
                maxTempsTimePoint{j} = 0;
            else
                maxTempsTimePoint{j} = double(val);  % or single(val), if you're using single precision
            end
        end


        maxSimTemps = cell2mat(maxTempsTimePoint);
        maxTa = cell2mat(maxTempsTimePointMR);
        deltaTemps = maxTa - 37;

        all_maxSimTemps = [all_maxSimTemps; maxSimTemps(:)];
        all_deltaTemps = [all_deltaTemps; deltaTemps(:)];
    end

    subplot(rows, cols, i);
    hold on;
    plot(all_maxSimTemps, all_deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'blue');
    xLimits = [0, 40];
    yLimits = [0, 40];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel('Simulation Temperature Rise');
    ylabel('MR Thermometry Temperature Rise');
    title(KRXfil);
    cd('..');
end