    pOutCorrected = pressure; 


    pOutUncorrected = pressure2; 

    h = figure;
    ax1 = subplot(311);
    ax2 = subplot(312);
    ax3 = subplot(313);
    ii = 1;

    xLg{ii} = 'Uncorrected';
    xLg{ii+1} = 'Corrected';
    yLg{ii} = 'Uncorrected';
    yLg{ii+1} = 'Corrected';
    zLg{ii} = 'Uncorrected';
    zLg{ii+1} = 'Corrected';
    
    x = linspace(-size(pOutCorrected,2)/2,size(pOutCorrected,2)/2,size(pOutCorrected,2));
    y = linspace(-size(pOutCorrected,1)/2,size(pOutCorrected,1)/2,size(pOutCorrected,1));
    z = linspace(-size(pOutCorrected,3)/2,size(pOutCorrected,3)/2,size(pOutCorrected,3));
    fs = 0;
    fs = fs*1e3;
    z = z - z(Foc_ijk(3));

    plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','--','color',ii,'focus',fs);
    xlim([-50,50])
    plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','--','color',ii,'focus',fs);
    xlim([-50,50])
    plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','--','color',ii,'focus',fs);
    xlim([-50,50])

    % axes(ax1)
    % plot(x,curModl(loc(1),:,loc(3))*1e4,'k-');
    
    plotFocus2(x,y,z,pOutCorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','-','color',ii,'focus',fs);
    xlim([-50,50])
    plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','-','color',ii,'focus',fs);
    xlim([-50,50])
    plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','-','color',ii,'focus',fs);
    xlim([-50,50])
    axes(ax1);
    legend(xLg);
    ax1.ColorOrderIndex = ii;

    %plot([fs(1),fs(1)],[0,1e5],':','linewidth',2);
    title('All Plates');
    axes(ax2);
    legend(yLg);
    ax2.ColorOrderIndex = ii;
    %plot([fs(2),fs(2)],[0,1e5],':','linewidth',2);
    axes(ax3);
    legend(zLg);
    ax3.ColorOrderIndex = ii;

    xLg{ii} = 'Uncorrected';
    xLg{ii+1} = 'Corrected';
    yLg{ii} = 'Uncorrected';
    yLg{ii+1} = 'Corrected';
    zLg{ii} = 'Uncorrected';
    zLg{ii+1} = 'Corrected';
    

