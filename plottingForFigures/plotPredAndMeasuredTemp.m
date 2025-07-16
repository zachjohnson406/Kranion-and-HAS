
function [] = plotPredAndMeasuredTemp(pred,measured)
    maxSimTemps = pred;
    maxTa = measured;
    
    hold on;
    plot(maxSimTemps, maxTa - 37, 'o', 'MarkerSize', 10, 'Color', 'blue');
    
    % Linear regression
    p = fitlm(maxSimTemps, maxTa - 37); % p(1) is the slope, p(2) is the intercept
    
    % Plot the linear regression line on the scatter plot
    xValues = linspace(0, 37, 100);  % Create 100 points for smooth line
    yValues = p.Coefficients.Estimate(1) + p.Coefficients.Estimate(2) * xValues;  % y = mx + b
    plot(xValues, yValues, 'r-', 'LineWidth', 2);  % Plot the regression line in red
    
    rsquared = p.Rsquared.Ordinary;
    pvalue = p.Coefficients.pValue(2);
    grid on;
    axis image;
    xLimits = [0, 30];
    yLimits = [0, 30];
    xlim(xLimits);
    ylim(yLimits);
    xlabel(['Simulation Temperature', char(176),'C'], 'FontSize',14);
    ylabel(['MR Thermometry Temperature', char(176), 'C'], 'FontSize', 14);
    xLine = [0 max(xLimits)];
    yLine = xLine;
    plot(xLine, yLine, 'k-', 'LineWidth', 1);
end