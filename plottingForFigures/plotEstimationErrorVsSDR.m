function [] = plotEstimationErrorVsSDR(maxTempFile, tempOverTimeFile, firstTimePoint, lastTimePoint)

kranionExports = {'Patient17','Patient22','Patient24','Patient27','Patient28','Patient30','Patient31','Patient37','Patient41'};
sdr = [0.57, 0.39, 0.38, 0.68, 0.37, 0.49, 0.49, 0.48, 0.47];
beamIndex = [92, 27.81,26.51,137.79,105.09,74.52,50.96,59.04,47.94];
acqTime = 3.456;
cd('/Users/zjohnson/Documents/MATLAB/KranionExports');

meanErrors = zeros(1, length(kranionExports));

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

            idx_before = find(timevec <= targetTime, 1, 'last');
            idx_after = find(timevec >= targetTime, 1, 'first');

            if isempty(idx_before) || isempty(idx_after)
                interpolated_temp = zeros(size(temps(:,:,:,1)));
            elseif idx_before == idx_after
                interpolated_temp = temps(:,:,:,idx_before);
            else
                t1 = timevec(idx_before);
                t2 = timevec(idx_after);
                w1 = (t2 - targetTime) / (t2 - t1);
                w2 = (targetTime - t1) / (t2 - t1);
                interpolated_temp = w1 * temps(:,:,:,idx_before) + w2 * temps(:,:,:,idx_after);
            end

            max_temp = max(interpolated_temp(:));
            maxTempsTimePoint{jj} = max_temp;

            if targetTimePoint <= size(tempsAll{jj},3)
                tempsAllSlice = tempsAll{jj}(128-2:128+2,128-2:128+2,targetTimePoint);
                maxTempsTimePointMR{jj} = max(tempsAllSlice(:));
            else
                maxTempsTimePointMR{jj} = 37;
            end
        end

        for j = 1:length(maxTempsTimePoint)
            val = maxTempsTimePoint{j};
            if isempty(val) || ~isnumeric(val)
                maxTempsTimePoint{j} = 0;
            else
                maxTempsTimePoint{j} = double(val);
            end
        end

        maxSimTemps = cell2mat(maxTempsTimePoint);
        maxTa = cell2mat(maxTempsTimePointMR);
        deltaTemps = maxTa - 37;

        all_maxSimTemps = [all_maxSimTemps; maxSimTemps(:)];
        all_deltaTemps = [all_deltaTemps; deltaTemps(:)];
    end

    % Remove known bad data for Patient17
    if(strcmp(KRXfil,'Patient17'))
        all_maxSimTemps([4,25]) = [];
        all_deltaTemps([4,25]) = [];
    end

    % Compute mean signed error (positive = overestimate, negative = underestimate)
    errors = all_maxSimTemps - all_deltaTemps;
    meanErrors(i) = mean(errors);

    cd('..');
end

% Plot mean error vs SDR
figure;
sdr_and_beam = sdr + beamIndex;
sdr = sdr_and_beam;
scatter(sdr, meanErrors, 80, 'filled');
grid on;
xlabel('SDR', 'FontSize', 14);
ylabel('Mean Simulation Error (Sim - MR, °C)', 'FontSize', 14);
title('Correlation Between Estimation Error and SDR', 'FontSize', 16);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--'); % Reference line at zero
set(gca, 'FontSize', 12);

% Optional: Add correlation coefficient
R = corrcoef(sdr, meanErrors);
text(min(sdr)+0.01, max(meanErrors)-0.5, sprintf('r = %.2f', R(1,2)), 'FontSize', 13);

end
