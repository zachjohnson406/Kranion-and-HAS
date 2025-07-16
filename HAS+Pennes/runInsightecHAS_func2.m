function [Q_rel,pressar,phmat,AmpPhar,pHAS,pR] = runInsightecHAS_func2(fullfilepath, pHAS, Modl_in, pR)
%   runInsightecHAS_func runs the HAS calculation on each of the 7 Insightec
%   transducer sections
%
%   This function rotates the skull/head model to align with the 7
%   pre-computed ERFA planes.
%
% DLP: 12/11/20 changed to output final 1024 array of amplitudes and phases

load([fullfilepath filesep 'ERFA8.mat']);

% PhCoQ = questdlg('Phase Correction?', 'Correction', ...
%     'No Correction', 'Insightec Correction', 'Test Correction', 'No Correction');

test = 0;
% switch PhCoQ
%     case 'No Correction'
%         PhCo = 1;
%     case 'Insightec Correction'
%         PhCo = 2;
%     case 'Test Correction'
%         PhCo = 2; 
%         test = 1;
%     otherwise
%         error('Unexpected case');
% end
PhCo = 2;

Foc_ijk = pR.Foc_ijk;
padVal = 0;
Modl_in = padarray(Modl_in,[ones(1,3)]*padVal,1);
Foc_ijk = Foc_ijk + padVal;

Modl = Modl_in;
% [xd,yd,zd] = size(Modl);
% Modl = ones(xd,yd,zd);
% Modl_in = ones(xd,yd,zd);

pR.useRayTracing = false;
if PhCo==1 %no correction
    usePhaseCorrection = false;  % true/false
    usePhaseTimeRev = false;  % true/false
elseif PhCo==2 %Insightec phase correction
    usePhaseCorrection = true;  % true/false
    usePhaseTimeRev = false;  % true/false
else PhCo==3 % time reversal phase correction
    usePhaseCorrection = false;  % true/false
    usePhaseTimeRev = true;  % true/false
end

pR.usePhaseCorrection = usePhaseCorrection;
pR.usePhaseTimeRev = usePhaseTimeRev;
pHAS.padsize = padVal;

disp(['PhCo=',num2str(PhCo),' usePhaseCorrection=',num2str(usePhaseCorrection),' usePhaseTimeRev=',num2str(usePhaseTimeRev)]);

Powerv = 1;       
[xdnx,xdny,xdnz] = size(Modl);
xcent = xdnx/2; ycent = xdny/2;
phmat = 0;

a_abs = pHAS.a_abs;
%DLP 10/20/20 put put in the Modl for rotation purposes
if(PhCo == 3)
    sm = size(Modl);
    focusindv=round( (pHAS.geom(2) + pHAS.vmm)/pHAS.Dy + 0.5 + (sm(1)/2) );  % y index of new focus in padded Modlpd;
    focusindh=round( (pHAS.geom(1) + pHAS.hmm)/pHAS.Dx + 0.5 + (sm(2)/2) );  % distances in mm.
    % focusindz=round( (pHAS.geom(3) + pHAS.zmm)/pHAS.Dz);
    focusindz = Foc_ijk(3);
    Modl_in(focusindv,focusindh,focusindz) = pHAS.indtarget; %stick with brain properties at the focus 
end

pressar = single(zeros(xdnx,xdny,xdnz,7));
AmpPhFull = single(zeros(128,2,8));


for jXducerSection = 1:7
    %load ERFA file (jX)
    fname = fullfilepath; 

    pERFA.perfa = perfa_seg{jXducerSection};
    pERFA.fMHz = fMHz;
    pERFA.Len = Len;
    pERFA.sx = sx;
    pERFA.R = R;
    pERFA.isPA = isPA;
    pERFA.ElemLoc = Seg(jXducerSection).ElemLoc_Polar;
    pERFA.dxp = dxp;
    pERFA.dyp = dyp;
    pERFA.relem = relem;
    
    % Transducer power for sector
    if jXducerSection == 1
        Pr = Powerv/4;
    else 
        Pr = Powerv/8;
    end
   
    % load in each phase correction for each zone
    if usePhaseCorrection == true || usePhaseTimeRev == true
        fname = [fullfilepath filesep 'ERFA8.mat'];       
    end

    pR.phaseFullFile = fname;
    pR.jXducerSection = jXducerSection;
    pR.Pr = Pr;       
    
    TH_rot = -rad2deg(TH_rotation);
    PHI_rot = -rad2deg(PHI_rotation);
    
    Modl = prepHAS_Model(Modl_in,PHI_rot(jXducerSection), -TH_rot(jXducerSection),Foc_ijk,-pHAS.xtilt);

    %run HAS-nogui
    if(PhCo == 3)
        [vmax, locmax] = max(Modl(:));
%       Modl(locmax) = Modl(locmax) - 100;
        [focusindv,focusindh,focusindz] = ind2sub(size(Modl),locmax);           %This is why we put the target in here, as the maximum value
        focusindv=focusindv + pHAS.padsize ;  % y index of new focus in padded Modlpd;
        focusindh=focusindh + pHAS.padsize ;  % y index of new focus in padded Modlpd;
           disp(['focusindv=',num2str(focusindv),' focusindh=',num2str(focusindh),' focusindz=',num2str(focusindz)]);
    end

    pR.ERFA_load = 0;
    pR.modl_load = 0;
    pR.HAS_savefile = 0;
    pR.f = pERFA.fMHz * 1e6;

    
    pHAS.vmm = pHAS.SteeredFocus(2);
    pHAS.hmm = pHAS.SteeredFocus(1);
    pHAS.zmm = pHAS.SteeredFocus(3);

    pHAS.dmm = 150-Foc_ijk(3) -1;
    
    %********************************************************
    %pHAS_NoGUIFull;         % THIS IS THE MAIN HAS PROGRAM !!!
    [pout,Z,angpgvect,phasematsv,pR,pHAS] = pHAS_NoGUIFullfunc(Modl,pHAS,pR,pERFA);
    %********************************************************
    if test ==1
        testPhCo; 
    end
    
    if jXducerSection == 1
        Z_modl=Z;
        %absmodl=pHAS.a_abs(Modl_in)*1e2*fMHz;     % DLP added 10/12/20 a(i) is pressure total attenuation coefficient (assume linear freq dep).
        %absmodl_modl=a_abs(Modl)*1e2*pERFA.fMHz;    % aabs(i) is pressure absorption coefficient (no random variation in it now).
        %absmodl_modl=absmodl;
    end
    
    pout=rotvolpivrecenterinterp(pout,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
   if test == 1
    pout2 = rotvolpivrecenterinterp(pout2,[ceil(xcent),ceil(ycent),Foc_ijk(3)],1,1,1,PHI_rot(jXducerSection),TH_rot(jXducerSection),0,1);
   end 
   if jXducerSection == 1
        pressure = pout;
        if test == 1
            pressure2 = pout2;
        end 
            %DLP debug
        if(PhCo == 3)
            [npx,npy,nel] = size(phasematsv);
            AmpPhFull(:,:,1) = pR.AmpPh1(1:128,:);
            AmpPhFull(:,:,2) = pR.AmpPh1(129:256,:);
            phmat = single(zeros(npx,npy,128,8));
            phmat(:,:,:,1) = phasematsv(:,:,1:128);
            phmat(:,:,:,2) = phasematsv(:,:,123:250);
        end
    else
        if(PhCo == 3)
            [~,~,nel] = size(phasematsv);
            AmpPhFull(:,:,jXducerSection+1) = pR.AmpPh1;
            phmat(:,:,1:nel,1+jXducerSection) = phasematsv;
        end
        pressure = pressure + pout;
        if test == 1
            pressure2 = pressure2 + pout2;
        end 
    end

    pressar(:,:,:,jXducerSection) = pout;
    
    if(PhCo == 3)
        fname=['AmpPh2_Z',num2str(jXducerSection),'_s',num2str(pR.json),'_c',num2str(pHAS.patnum),'.mat'];
        AmpPh1 = pR.AmpPh1;
        save(fname,'AmpPh1');
    end
end

if(PhCo == 3)
   AmpPhar = [squeeze(AmpPhFull(:,:,1));squeeze(AmpPhFull(:,:,6));squeeze(AmpPhFull(:,:,2));squeeze(AmpPhFull(:,:,3)); ...
                squeeze(AmpPhFull(:,:,7));squeeze(AmpPhFull(:,:,8));squeeze(AmpPhFull(:,:,4));squeeze(AmpPhFull(:,:,5))];
            fname = (['AmpPhar-PatNum',num2str(pHAS.patnum),'son',num2str(pR.json),'vel',num2str(pHAS.c(7)),'.mat']);
            save(fname,'AmpPhar');
else
    AmpPhar = [];
end
%    zpa1 = [ampphar(1:128,:);ampphar(257:384,:)];
%     zpa2 = ampphar(385:512,:);
%     zpa3 = ampphar(769:896,:);
%     zpa4 = ampphar(897:1024,:);
%     zpa5 = ampphar(129:256,:);
%     zpa6 = ampphar(513:640,:);
%     zpa7 = ampphar(641:768,:);

if test ==1
    plotPhCo; 
end
%Q_rel=abs((pressure.*conj(pressure)).*absmodl_modl./Z_modl);        %Note the absmodl should be pressure attenuation DLP 10/12/20
Modl_Q = prepHAS_Model(Modl_in,0, 0,Foc_ijk,-pHAS.xtilt);
Q_rel=abs((pressure.*conj(pressure)).*a_abs(Modl_Q)*1e2*pERFA.fMHz./Z_modl);        %DLP 1/17/22 no need for absmodl_modl for just one use. save some memory
%Q_rel=abs((pressure.*absmodl_modl.*conj(pressure.*absmodl_modl))./Z_modl);        %Note the absmodl should be pressure attenuation DLP 10/12/20

end