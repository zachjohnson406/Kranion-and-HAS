% plotFocus is a helper function that plots a focus along the x, y, and z
% axes.

function [loc,loc_ijk,xAx,yAx,zAx] = plotFocus2(x,y,z,p,varargin)%,plotAxs,h,ax,ls,colorIdx)

%% Parse inputs
if mod(length(varargin),2)
    error('Expected type/argument pairs')
end

ls = '-';
plotAxs = [1,1,1];
for ii = 1:2:length(varargin)
    switch(lower(varargin{ii}))
        case 'plotaxs'
            plotAxs = varargin{ii+1};
        case 'h'
            h = varargin{ii+1};
        case 'ax'
            ax = varargin{ii+1};
        case 'linestyle'
            ls = varargin{ii+1};
        case 'color'
            colorIdx = varargin{ii+1};
        case 'window'
            window = varargin{ii+1};
        case 'focus'
            focus = varargin{ii+1};
        otherwise
            error([varargin{ii}, ' is not a recognized parameter'])
    end
end
%% Set defaults for values not set
if ~exist('h','var')
    h = figure;
end
figure(h);
%h.Position = [-666    46   667   976];
if ~exist('ax','var')
    ax = gca;
end
axes(ax);
if ~exist('colorIdx','var')
    colorIdx = ax.ColorOrderIndex;
end

if exist('window','var') && ~exist('focus','var')
    error('if you define a window you must define an intended focal location')
end

if exist('window','var')
    tmp = zeros(size(p));
    [~,fx] = min(abs(x-focus(1)));
    [~,fy] = min(abs(y-focus(2)));
    [~,fz] = min(abs(z-focus(3)));
    tmp((fy-window):(fy+window),(fx-window):(fx+window),(fz-window):(fz+window)) = ...
        p((fy-window):(fy+window),(fx-window):(fx+window),(fz-window):(fz+window));
    
    p = tmp;
end

[~,idx] = max(abs(p(:)));
[a,b,c] = ind2sub(size(p),idx);

disp(['<', num2str(y(a)), ',', num2str(x(b)), ',' num2str(z(c)),'>'])
loc = [y(a);x(b);z(c)];
loc_ijk = [a,b,c];

axes(ax);
hold on
if ~exist('colorIdx','var')
    colorIdx = 1;
end
ax.ColorOrderIndex = colorIdx;
labIdx = 1;
lblx = 0;
if ~exist('lbl','var')
    lbl = cell(1,sum(plotAxs));
    lblx = 1;
end
if plotAxs(1)
    xAx = squeeze(abs(p(a,:,c)));
    plot(x,xAx,ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'x';
        labIdx = labIdx+1;
    end
end
if plotAxs(2)
    yAx = squeeze(abs(p(:,b,c)));
    plot(y,yAx,ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'y';
        labIdx = labIdx+1;
    end
end
if plotAxs(3)
    zAx = squeeze(abs(p(a,b,:)));
    plot(z,zAx,ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'z';
    end
end
% legend(lbl)
%makeFigureBig(h)