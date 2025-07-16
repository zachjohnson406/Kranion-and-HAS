% FUSF_HAS_Prep
%
% A first-run function for a new user/site. This is currently designed
% for use with Insightec hardware / treatment exports, and does the
% following:
% - *Calls loadKranionScene() to load parsed Kranion export into memory
% - *Segments Insightec transducer into 7 sections for use with HAS
% - *Computes the Rayleigh-Sommerfeld integral for each segment
% - *Saves the "ERFA" planes and associated tilt angles required by HAS
%
% - *Generates the Modl file based on the CT volume
% - *Executes 7-segment HAS simulation
% - *Saves pressure and tempature field into Kranion-readable format
%
%
% @OUTPUTS
%   Q_rel: Q pressure pattern in Kranion coordinate space
%   pressar: Pressure in Kranion coordinate space
%   temps: Predicted tempature feild in Kranion coordinate space
%
%
%   Taylor Webb, Zach Johnson, and Matt Eames
%   June 19th, 2024

function [Q_rel,pressar,temps] = FUSF_HAS_Prep_loop_new_values(KRXfile,KRXpath)

%% Load Kranion Scene
[loc, Vhas, NumSon, toKranion] = loadKranionScene_loop(KRXfile, KRXpath);

fullfilepath=[KRXpath KRXfile];
ctinxd = Vhas;
%% Create tissue segmentation and 'phantom' model for HAS
Vhas = FUSF_Make_Modl(Vhas, fullfilepath);

%% ERFA Process, skip if saved already
if isfile([fullfilepath filesep 'ERFA8.mat'])
    disp('ERFA8.mat file found, skipping computations...')
    load([fullfilepath filesep 'ERFA8.mat'])
    disp('ERFA8.mat file finished loading');
else
    [Seg, TH_rotation, PHI_rotation] = segmentTransducerElements_taylor(loc);
    FUSF_ERFAMaker(Seg,TH_rotation,PHI_rotation,loc,fullfilepath);
end


[pHAS] = FUSF_LoadHASphantom(fullfilepath,Vhas);
pHAS.ctinxd = ctinxd;
Vhas = cleanHASModel(Vhas,pHAS);
fname = [fullfilepath filesep 'Modl.mat'];
save(fname,'Vhas');
maxTempsSimulation = {};
for i = 1:NumSon
    disp(['Simulating Sonication #',num2str(i)]);
    try
        load([fullfilepath filesep ['SonicationParameters',num2str(i),'.mat']]);
    catch
        disp(['no data for sonication', num2str(i)]);
        continue;
    end
    phase = Sonication.Phase;
    amplitude = Sonication.Amplitude;
    %%
    duration = Sonication.Duration;
    pHAS.SteeredFocus = [0,0,0];
    pHAS.xtilt = Sonication.XTilt;

    %% Create Phantom, run HAS, and compute tempature
    [pR] = matchPhaseToSections(Seg, phase,amplitude);
    pR.Foc_ijk = FocCTijk;

    pHAS.rho = [1000, 1007, 911, 1109, 1046, 1908, 1178];
    pHAS.c = [1500,1505,1440,1624,1546,2814,2118];
    %   water,  CSF,      fat,    skin,   brain,  bone,   diploe
    %%CSF paper values
    pHAS.a = [0,     0.068,  0.044,  0.21,    0.068,  .54553,    .47];        
    pHAS.a_abs = [0, 0.012,  0.044,  0.21,    0.068,  .54553,    .47];
    %%CSF paper values with average brain attenuation and absorption
    pHAS.a = [0,     0.068,  0.044,  0.21,    0.0395,  .54553,    .47];        
    pHAS.a_abs = [0, 0.012,  0.044,  0.21,    0.0395,  .54553,    .47];


%     load('maxTempsSimulationOriginalValues.mat');
%     load('MRTIdata');
%     tempDiff = maxTempMRTI(3) - maxTempSimulation(3);
%     percentDiff = ((tempDiff/(maxTempMRTI(3) +  maxTempMRTI(3)/2)));
%     pHAS.a_abs(5) = 0.012 * (1 + percentDiff); 
%     pHAS.a(5) = 0.012 * (1+percentDiff);
     [Q_rel,pressar,phmat,AmpPhar,pHAS,pR] = runInsightecHAS_func2(fullfilepath, pHAS, Vhas, pR);

    fname = ['son',num2str(i),'_var.mat'];
    try
        load(fname,'power_acou','timeDuration');
    catch
        continue;
    end
    [xdim,ydim,zdim] = size(Q_rel);
    xcent = round(xdim/2); ycent = round(ydim/2);
    [temps,timevec] = runPennes_func(Q_rel(xcent-10:xcent+10,ycent-10:ycent+10,:),power_acou,timeDuration,pHAS);

    %% Run Pennes, convert to Kranion coordinate space, save outputs
    %pout = sum(pressar,4);
    %cd(fullfilepath)
    %load("CT_Volume.mat",'unTransform');
    %load("CT_Volume.mat",'Size');
    %pressureInKranion = SaveToKranion(pout,pR.Foc_ijk,pHAS.xtilt,unTransform,Size(3),KRXfile,toKranion, 'pressure');
    %Q_relInKranion = SaveToKranion(Q_rel,pR.Foc_ijk,pHAS.xtilt,unTransform,Size(3),KRXfile,toKranion, 'Q_rel');
    [maxTempValue, index] = max(temps(:));
    [~,~,~,t] = ind2sub(size(temps), index);
    %tempInKranion = SaveToKranion(temps(:,:,:,t),pR.Foc_ijk,pHAS.xtilt,unTransform,Size(3),KRXfile,toKranion, 'temp');
    %save(['SonicationData',num2str(i),'.mat'],'temps', 'timevec', 'Modlinxd','pressar','Q_rel');
    maxTempsSimulation{i} = maxTempValue; 
   
end
maxSimTemps = cell2mat(maxTempsSimulation);
save('maxTempsSimulationCSFPaperValues.mat','maxTempsSimulation');

end


