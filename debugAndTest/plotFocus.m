% plotFocus is a helper function that plots a focus along the x, y, and z
% axes.

function [loc,loc_ijk] = plotFocus(x,y,z,p,plotAxs,h,ax,ls,colorIdx)

if ~exist('plotAxs','var')
    plotAxs = ones(1,3);
end

if ~exist('h','var')
    h = figure;
end

figure(h);
% h.Position = [1252          87         667         890];
%h.Position = [-666    46   667   976];
if ~exist('ax','var')
    ax = gca;
end
if ~exist('ls','var')
    ls = '-';
end

[~,idx] = max(abs(p(:)));
[a,b,c] = ind2sub(size(p),idx);

%disp(['<', num2str(y(b)), ',', num2str(x(a)), ',' num2str(z(c)),'>'])
loc = [y(b);x(a);z(c)];
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
    plot(x,squeeze(abs(p(a,:,c))),ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'x';
        labIdx = labIdx+1;
    end
end
if plotAxs(2)
    plot(y,squeeze(abs(p(:,b,c))),ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'y';
        labIdx = labIdx+1;
    end
end
if plotAxs(3)
    plot(z,squeeze(abs(p(a,b,:))),ls,'linewidth',2);
    if lblx
        lbl{labIdx} = 'z';
    end
end
legend(lbl)
makeFigureBig(h)