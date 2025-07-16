function [ jlast ] = stepimg( img,intscale,dostep )
%stepimg does loop to display
%   Detailed explanation goes here
%DLP 12/16/20 added output of current image displayed


minval = min(img(:));
maxval = max(img(:));
set(gca, 'FontSize', 16);
tickLabels = {'water', 'fat', 'brain', 'bone', 'diploe'};

% Adjust the colormap to 5 colors
colormap(parula(5));

% Plot the image
imagesc(img, [minval maxval]);
axis image;
axis off;
colorbar;

% Adjust the colorbar
cb = colorbar('Position', [.75 .11 .0675 .815]);
cb.TickLabels = tickLabels;

% Set the tick positions to match the number of categories
numTicks = length(tickLabels);
cb.Ticks = linspace(minval, maxval, numTicks);

% Update axis and title
set(gca, 'FontSize', 20, 'Position', [.07 .11 .685 .815]);
title('', 'FontSize', 16);
