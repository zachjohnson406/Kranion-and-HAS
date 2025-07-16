
    %pOutCorrected = pOutCorrected.pout;
    pOutUncorrected = abs(sum(pressarUnCo,4));
    pOutCorrected = abs(sum(pressar,4));
    %pOutUncorrected = pOutUncorrected.pout;
    %pOutCorrected = sum(pOutCorrected,4);
    %pOutUncorrected = sum(pOutUncorrected,4);

    
    h = figure;
    ax1 = subplot(311);
    ax2 = subplot(312);
    ax3 = subplot(313);
    ii = 1;
    jXducerSection = 1;
    x = linspace(-size(pOutCorrected,2)/2,size(pOutCorrected,2)/2,size(pOutCorrected,2));
    y = linspace(-size(pOutCorrected,1)/2,size(pOutCorrected,1)/2,size(pOutCorrected,1));
    z = linspace(-size(pOutCorrected,3)/2,size(pOutCorrected,3)/2,size(pOutCorrected,3));
    fs = 0;
    fs = fs*1e3;
    %z = z - z(Foc_ijk(3)+30);

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
    ax1.ColorOrderIndex = ii;
    %plot([fs(1),fs(1)],[0,1e5],':','linewidth',2);
    title(['Composite Pressure']);
    axes(ax2);
    ax2.ColorOrderIndex = ii;
    %plot([fs(2),fs(2)],[0,1e5],':','linewidth',2);
    axes(ax3);
    ax3.ColorOrderIndex = ii;
    %plot([fs(3),fs(3)],[0,1e5],':','linewidth',2);

    % xLg{ii} = ['xSteering = ', num2str(allSons.Sonication(sonNo).SteeredFoc(1))];
    % yLg{ii} = ['ySteering = ', num2str(allSons.Sonication(sonNo).SteeredFoc(2))];
    % zLg{ii} = ['zSteering = ', num2str(allSons.Sonication(sonNo).SteeredFoc(3))];

    xLg{ii} = 'Uncorrected';
    xLg{ii+1} = 'Corrected';
    yLg{ii} = 'Uncorrected';
    yLg{ii+1} = 'Corrected';
    zLg{ii} = 'Uncorrected';
    zLg{ii+1} = 'Corrected';

    % With time reversal correction
    % With correction
    % tmpPR = pR;
    % 
    % tmpPR.usePhaseCorrection = false;
    % tmpPR.usePhaseTimeRev = true;
    % curModl(xCentIdx,yCentIdx,Foc_ijk(3)) = 100;
    % [pOutCorrectedTr,~,angCorrectedTr] = pHAS_NoGUIFullfunc(curModl,pHAS,tmpPR,pERFA);
    % % pOutCorrectedTr = rotvolpivrecenterinterp(pOutCorrectedTr,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
    % 
    % plotFocus(x,y,z,pOutCorrectedTr,[1,0,0],h,ax1,':',ii);
    % plotFocus(x,y,z,pOutCorrectedTr,[0,1,0],h,ax2,':',ii);
    % plotFocus(x,y,z,pOutCorrectedTr,[0,0,1],h,ax3,':',ii);    

    % phV(ampV==0) = nan;
    % phV = phV-mean(phV,'omitnan');
    % h2 = plotPhaseSpatial(Seg(jXducerSection).ElemLoc_Cart(:,1),Seg(jXducerSection).ElemLoc_Cart(:,2),Seg(jXducerSection).ElemLoc_Cart(:,3),phV,[-pi,pi]);
    % h2.Position = [584    87   667   800];
    % title('Kranion Phases')

    % ang(ampV==0) = nan;
    % ang = ang-mean(ang,'omitnan');
    % h2 = plotPhaseSpatial(Seg(jXducerSection).ElemLoc_Cart(:,1),Seg(jXducerSection).ElemLoc_Cart(:,2),Seg(jXducerSection).ElemLoc_Cart(:,3),ang,[-pi,pi]);
    % h2.Position = [584    87   667   800];
    % title('Steering Phases')

%%

    
%     pOutCorrected_load = load('pressar_HASGPi-L.mat');
%     pOutCorrected_all = pOutCorrected_load.pressar;
%     
%     pOutUncorrected_load = load('pressar_HASGPi-LnoPhCo.mat');
%     pOutUncorrected_all  = pOutUncorrected_load.pressar;
% %%
%     
%     pOutUncorrected = sum(pOutUncorrected_all,4);
%     pOutCorrected = sum(pOutCorrected_all,4);
% 
%     h = figure; 
%     
%     lgnd = {'Corrected','Uncorrected','SteeredFocus'};
%     ax1 = subplot(311);
%     ax2 = subplot(312);
%     ax3 = subplot(313);
%     x = linspace(-size(pOutCorrected,2)/2,size(pOutCorrected,2)/2,size(pOutCorrected,2));
%     y = linspace(-size(pOutCorrected,1)/2,size(pOutCorrected,1)/2,size(pOutCorrected,1));
%     z = linspace(-size(pOutCorrected,3)/2,size(pOutCorrected,3)/2,size(pOutCorrected,3));
%     z = z - z(Foc_ijk(3));
%     
%     fs = [Sonication.SteeredFoc(1), Sonication.SteeredFoc(2),Sonication.SteeredFoc(3)];
%     %fs = [str2num(Sonication.NatFoc{2}), str2num(Sonication.NatFoc{1}),str2num(Sonication.NatFoc{3})];
%     %fs = fs*1e-3;
%     
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','-','color',1,'focus',fs);
%     %xlim([-50,50])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','-','color',1,'focus',fs);
%     xlim([-50,50])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','-','color',1,'focus',fs);
%     xlim([-50,50])
%     legend(lgnd);
%     
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','--','color',2,'focus',fs);
%     line([fs(1), fs(1)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %title("Export: HAStest05172024 . With Kranion bone values. ")
%     xlim([-50,50])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','--','color',2,'focus',fs);
%     line([fs(2) fs(2)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     xlim([-50,50])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','--','color',2,'focus',fs);
%     line([fs(3) fs(3)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     xlim([-50,50])
%     legend(lgnd);
% 
% %% 
%  for i = 1:7
% 
%     h = figure; 
%     pOutUncorrected = pOutUncorrected_all(:,:,:,i);
%     pOutCorrected = pOutCorrected_all(:,:,:,i);
%     
%     lgnd = {'Corrected','Uncorrected','SteeredFocus'};
%     ax1 = subplot(311);
%     ax2 = subplot(312);
%     ax3 = subplot(313);
%     x = linspace(-size(pOutCorrected,2)/2,size(pOutCorrected,2)/2,size(pOutCorrected,2));
%     y = linspace(-size(pOutCorrected,1)/2,size(pOutCorrected,1)/2,size(pOutCorrected,1));
%     z = linspace(-size(pOutCorrected,3)/2,size(pOutCorrected,3)/2,size(pOutCorrected,3));
%     z = z - z(Foc_ijk(3));
%     
%     fs = [Sonication.SteeredFoc(1), Sonication.SteeredFoc(2),Sonication.SteeredFoc(3)];
%     %fs = [str2num(Sonication.NatFoc{2}), str2num(Sonication.NatFoc{1}),str2num(Sonication.NatFoc{3})];
%     %fs = fs*1e-3;
%     
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])
%     legend(lgnd);
%     
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','--','color',2,'focus',fs);
%     line([fs(1), fs(1)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %xlim([-20,20])
%     legend(lgnd);
%     title(["Plate:", num2str(i)]);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','--','color',2,'focus',fs);
%     line([fs(2) fs(2)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','--','color',2,'focus',fs);
%     line([fs(3) fs(3)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %xlim([-20,20])
%     legend(lgnd);
% end
% 
% %%
%     
%     pOutUncorrected = pout;
% 
%     pR.usePhaseCorrection = true;
%     pR.usePhaseTimeRev = false;
%     [pOutCorrected,Z,angpgvect,phasematsv,pR,pHAS] = pHAS_NoGUIFullfunc(Modl,pHAS,pR,pERFA);
%     pOutCorrected=rotvolpivrecenterinterp(pOutCorrected,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
%     pR.usePhaseCorrection = false;
%     pR.usePhaseTimeRev = false; 
%     
% 
%     h = figure; 
%     
%     lgnd = {'Corrected','Uncorrected','SteeredFocus'};
%     ax1 = subplot(311);
%     ax2 = subplot(312);
%     ax3 = subplot(313);
%     x = linspace(-size(pOutCorrected,2)/2,size(pOutCorrected,2)/2,size(pOutCorrected,2));
%     y = linspace(-size(pOutCorrected,1)/2,size(pOutCorrected,1)/2,size(pOutCorrected,1));
%     z = linspace(-size(pOutCorrected,3)/2,size(pOutCorrected,3)/2,size(pOutCorrected,3));
%     z = z - z(Foc_ijk(3));
%     
%     fs = [Sonication.SteeredFoc(1), Sonication.SteeredFoc(2),Sonication.SteeredFoc(3)];
%     %fs = [str2num(Sonication.NatFoc{2}), str2num(Sonication.NatFoc{1}),str2num(Sonication.NatFoc{3})];
%     %fs = fs*1e-3;
%     
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])s
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutCorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','-','color',1,'focus',fs);
%     %xlim([-20,20])
%     legend(lgnd);
%     
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[1,0,0],'h',h,'ax',ax1,'linestyle','--','color',2,'focus',fs);
%     %line([fs(1), fs(1)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     title(["Plate: ", num2str(jXducerSection)]);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,1,0],'h',h,'ax',ax2,'linestyle','--','color',2,'focus',fs);
%     %line([fs(2) fs(2)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %xlim([-20,20])
%     legend(lgnd);
%     plotFocus2(x,y,z,pOutUncorrected,'plotAxs',[0,0,1],'h',h,'ax',ax3,'linestyle','--','color',2,'focus',fs);
%     %line([fs(3) fs(3)], ylim, 'Color', 'r', 'LineStyle', ':','LineWidth', 2);
%     %xlim([-20,20])
%     legend(lgnd);