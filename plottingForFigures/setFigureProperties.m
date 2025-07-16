function [] = setFigureProperties(fig,figName)
    fname = figName;
    picturewidth = 25; % set this parameter and keep it forever
    hw_ratio = 0.65; % feel free to play with this ratio
    %set(findall(fig,'-property','FontSize'),'FontSize',17); % adjust fontsize to your document
    set(findall(fig,'-property','Box'),'Box','off'); % optional
    set(findall(fig,'-property','Interpreter'),'Interpreter','latex');
    set(findall(fig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex');
    set(fig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth]);
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);
    %cd('/Users/zjohnson/Documents/MATLAB/KranionFigures');
    %print(fig,fname,'-dpdf','-painters','-fillpage')
    %print(hfig,fname,'-dpng','-painters')
end 