function [] = plotMeasuredVsSimulated(maxTempFile,tempOverTimeFile,firstTimePoint,lastTimePoint)
kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient41'};
titles = {'Patient 1','Patient 2','Patient 3','Patient 4','Patient 5','Patient 6','Patient 7','Patient 8','Patient 9'};
sdr = {0.57, 0.39, 0.38,0.68,0.37,0.49,0.49,0.48,0.47};
beamIndex = {92, 27.81,26.51,137.79,105.09,74.52,50.96,59.04,47.94};
acqTime = 3.456;
n = length(kranionExports);
cols = ceil(sqrt(n));
rows = ceil(n / cols);
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

for i = 1:length(kranionExports)
    KRXfil = kranionExports{i};
    cd(KRXfil);
    load('MRTIdata.mat');
    load(maxTempFile)

    all_maxSimTemps = [];
    all_deltaTemps = [];

    for targetTimePoint = firstTimePoint:lastTimePoint
        targetTime = (acqTime/2) + ((targetTimePoint - 2) * acqTime);
        maxTempsTimePoint = {};
        maxTempsTimePointMR = {};

        for jj = 1:length(maxTempsSimulation)
            load([tempOverTimeFile, num2str(jj), '.mat']);
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

    % After your plotting commands inside the loop:
    subplot(rows, cols, i);
    hold on;
    if(strcmp(KRXfil,'Patient17'))
        all_maxSimTemps(4) = [];
        all_deltaTemps(4) = [];
    end
    plot(all_maxSimTemps, all_deltaTemps, 'o', 'MarkerSize', 8, 'Color', 'blue', 'LineWidth',  1.1);
    xLimits = [0, 40];
    yLimits = [0, 40];
    plot(xLimits, xLimits, 'k-', 'LineWidth', 1);
    grid on;
    axis image;
    xlim(xLimits);
    ylim(yLimits);
    xlabel({'Simulation', 'Temperature Rise ($^\circ$C)'}, 'Interpreter', 'latex');
    ylabel({'MR Thermometry', 'Temperature Rise ($^\circ$C)'}, 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 11;
    ax.XTick = [10 20 30 40]; 
    ax.YTick = [10 20 30 40];  
    % Add a textbox with SDR and Beam Index
    infoStr = sprintf('Beam Index = %.2f\nSDR = %.2f', beamIndex{i},sdr{i});
    text(xLimits(1)+1, yLimits(2)-5, infoStr,'Interpreter','latex','FontSize',13);
    
    title(titles{i}, 'FontSize', 17);

    cd('..');
end
end