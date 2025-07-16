function [] = plotTempOverTime(temps,tempsAll,timevec,SonNum,x,y)  
    hold on;
    acqTime = 3.456;
    err = (acqTime/2)*ones(1,size(tempsAll{SonNum},3));
    MRTItempOneSonication = zeros(1,size(tempsAll{SonNum},3));
    acqTimeMean = zeros(1,size(tempsAll{SonNum},3));
    for ii = 1:size(tempsAll{SonNum},3)
        columns = x{SonNum} - 1 : x{SonNum} + 1; 
        rows = y{SonNum} - 1 : y{SonNum} + 1;
        subRegion = tempsAll{SonNum}(columns,rows,ii);
        MRTItempOneSonication(ii) = max(subRegion(:));
        acqTimeMean(ii) = ((ii-1) + 0.5) * (acqTime);
    end 
    
    
    maxTempsSonicationX = zeros(size(temps,4),1);
    for i = 1:length(maxTempsSonicationX)
        slice = temps(:,:,:,i);
        maxTempsSonicationX(i) = max(slice(:));
    end 
    timevec = timevec-1;
    
    plot(timevec, maxTempsSonicationX, '-', 'LineWidth',1.1);  % Thicker line
    hold on
    plot([0, acqTimeMean(2:end)-acqTime], [0, MRTItempOneSonication(2:end)-37], 'o-');  % Larger markers and thicker lines
    ax = gca;
    legend('Simulation', 'MRTI','FontSize',11);
    xlabel('Time (Seconds)','FontSize',13);
    ylabel('Temperature rise ($^\circ$C)', 'FontSize',13);
    xlim([0,30])
    grid on;
    title('Tempature rise over time','FontSize', 17)
    axis square;
end