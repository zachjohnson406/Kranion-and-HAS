function pHAS = FUSF_LoadHASPhantom(fullfilepath,Vhas)

load([fullfilepath filesep 'ERFA8.mat']);

c0 = 1500;
rho0 = 1000;

% Model voxel spacing
Dx = 1.0;       
Dy = 1.0;       
Dz = 1.0;     

zfocus = 0;

% Model size values 
[xdnx,xdny,xdnz] = size(Vhas);

%% 150 - back of transducer. 84 - top of modl in z direction  
dmm = 150 - ceil(xdnx); 

modlName = [fullfilepath filesep 'Modl.mat'];

% indices
pHAS.indwater = 1;pHAS.indCSF = 2; pHAS.indfat = 3; pHAS.indskin = 4; pHAS.indbrain = 5; pHAS.indbone = 6; pHAS.inddiploe = 7; pHAS.indgray = 8; pHAS.indwhite = 9; pHAS.indtarget = 10;
pR.addDQA = 1;
if (pR.addDQA == 0)
    %threshold values for the Modl:
    pHAS.thrsh_bone = 2000;
    pHAS.thrsh_diploe = 1550; %85;  %200;
    pHAS.thrsh_brain = 1000;
    pHAS.thrsh_fat = 800;
    pHAS.thrsh_skin = 0;
    %Values for the human brain modeling:
    %   water,  CSF,      fat,    skin,   brain,  bone,   diploe, gray, white, target
    %December values:
    %     c = [1500,  1504.5,   1480,   1600,   1546.3, 2653,   1200, 1546.3]; %2653,1600]; %water,gland,fat,skin,brain,bone,diploe
    %     a =  [  0,   0.001,  0.044,   0.22,    0.068,  1.3,    0.8,  0.068];        %DLP added an 8th for the focal point for time reversal
    %     a_abs = [0,  0.001,  0.044,   0.22,    0.012,  1.3,    1.0,  0.012];
    %     rho = [1000,  1007,    937,    937,     1046, 1738,   1300,   1046];
    %Oct values:
    %             water,  CSF,      fat,    skin,   brain,  bone,   diploe, gray,   white,   target
    CSF_vals.c = [1500,  1504.5,   1440,   1624,   1546.3, 2653,     2653,  1540,  1552.5   1546.3]; %2653,1600]; %water,CSF,fat,skin,brain,bone,diploe
    CSF_vals.a = [0,     0.001,  0.044,  0.21,    0.068,  .54553,    .47,  0.069, 0.0671, 0.068];        %DLP added an 8th for the focal point for time reversal
    CSF_vals.a_abs = [0, 0.001,  0.044,  0.21,    0.012,  .54553,    .47,  0.012, 0.0671, 0.012];
    CSF_vals.rho = [1000, 1007,    911,   1109,    1046,    1908,   1178,   1046,  1046,   1046];

    Brain_vals.c = [1500,  1546.3,   1440,   1624,   1546.3, 2653,   2653,  1540,  1552.5 , 1546.3]; %2653,1600]; %water,brain,fat,skin,brain,bone,diploe
    Brain_vals.a = [0,     0.068,  0.044,  0.21,    0.068,  .54553,    .47,  0.069, 0.0671, 0.068];        %DLP added an 8th for the focal point for time reversal
    Brain_vals.a_abs = [0, 0.012,  0.044,  0.21,    0.012,  .54553,    .47, 0.012, 0.0671,  0.012];
    Brain_vals.rho = [1000, 1046,    911,   1109,     1046, 1908,   1178, 1046, 1046, 1046];
    
    c = CSF_vals.c;       %[1500,  1500,   1480,   1480,   1546.3, 2653,   2653, 1546.3]; %2653,1600]; %water,gland,fat,skin,brain,bone,diploe
    a = CSF_vals.a;     %[0,     0.00133,  0.044,  0.22,    0.068,  1.3,    1.3, 0.068];        %DLP added an 8th for the focal point for time reversal
    a_abs = CSF_vals.a_abs;     %[0, 0.00133,  0.044,  0.22,    0.012,  1.3,    1.3, 0.012];
    rho = CSF_vals.rho;     %[1000, 1000,    937,   937,     1046, 1738,   1300, 1046];
    
    %Pennes values for human studies:
    k_mm = 0.527*10^-3;     % Thermal conductivity W/mm/°C
    rho_mm = 1045*10^-9;    % Density kg/mm3
    sp_h = 3600;            % Specific heat   J/kg/°C
    W = 5*10^-9;    %15*10^-9;         % Pennes perfusion parameter 9.3 x 10^-9 kg/mm3/s  *** Henrik 15 x 10^-9 kg/m3/s to conv to ml/kg/min 60s/min * 1000kg/m3 10^6 m3/mm3
else
    %threshold values for the Modl:
    pHAS.thrsh_bone = 2000;
    pHAS.thrsh_diploe = 1550; %85;  %200;
    pHAS.thrsh_brain = 1000;
    pHAS.thrsh_fat = 800;
    pHAS.thrsh_skin = 0;
    %Values for the phantom - with skull modeling (attenuation at 1MHz):
    %   water,  CSF,      fat,    skin,   DQA,  bone,   diploe,  target
    %c = [1500,  1504.5,   1480,   1600,   1560.0, 2653,   1800, 1560.0]; %2653,1600]; %water,gland,fat,skin,brain,bone,diploe
    c = [1500,  1504.5,   1480,   1600,   1560.0, 2652.0,   2652.0, 1560.0];
    a =  [  0,   0.001,  0.044,   0.22,    0.042,  1.3,    1.3,  0.042];        %DLP added an 8th for the focal point for time reversal
    a_abs = [0,  0.001,  0.044,   0.22,    0.042,  1.3,    1.3,  0.042];
    rho = [1000,  1007,    937,    937,     1046, 1738,   1300,   1046];
    pR.noSkull = 0;
    if(pR.noSkull ==1)
        c(6) = c(1); c(7) = c(1);
        a(6) = a(1); a(7) = a(1);
        a_abs(6) = a_abs(1); a_abs(7) = a_abs(1);
        rho(6) = rho(1); rho(7) = rho(1);
    end
    %Pennes values for phantom scans:
    k_mm = 0.55*10^-3;     % Thermal conductivity W/mm/°C
    rho_mm = 1045*10^-9;    % Density kg/mm3
    sp_h = 3400;            % Specific heat   J/kg/°C
    W = 0;         % Pennes perfusion parameter 9.3 x 10^-9 kg/mm3/s  *** Henrik 15 x 10^-9 kg/m3/s
% end
pHAS.k_mm = k_mm;
pHAS.rho_mm = rho_mm;
pHAS.sp_h = sp_h;
pHAS.W = W;

%a = [0,     0.133,  0.091,  0.086,  0.068,  0.3,      0.35];

%Pr = 100; %Transducer powers loaded from file
offsetxmm = 0; % initial values of offset of xducer axis to model axis in mm.
offsetymm = 0;
hmm = 0; % initial values of focus shift in mm.
vmm = 0;
zmm = 0;
Power=[288.75, 144.375]; %transucer power Z1 and Z2-7

randvc=zeros(1,length(a));
usePhaseCorrection = false;  % true/false
usePhaseTimeRev = false;  % true/false
useRayTracing = false;  % true/false
numrefl = 1; % (number of reflections)
% these are the full pathname and filename of the output .mat files (if empty, no .mat file will be written)
% poutFullFile = 'D:\HAS_Outputs\poutOut.mat';
% QFullFile = 'D:\HAS_Outputs\QOut.mat';
% ppFullFile = 'D:\HAS_Outputs\ppOut.mat';
% relgeomfocus = [21, 21, 21];
% padsize = 100;
% checkPhaseFileValidity = 0;

% dx = Dx;    %dlp klug - not sure where these should have been defined
% dy = Dy;
% dz = Dz;

pHAS.modlName = modlName;
pHAS.xdnx = xdnx;
pHAS.xdny = xdny;
pHAS.xdnz = xdnz;
pHAS.Dx = Dx;
pHAS.Dy = Dy;
pHAS.Dz = Dz;
pHAS.hmm = 0;
pHAS.vmm = 0;
pHAS.zmm = 0;
pHAS.dmm = dmm;
if (pR.addDQA == 0)
    pHAS.CSF_vals = CSF_vals;
    pHAS.Brain_vals = Brain_vals;
end
pHAS.c = c;
pHAS.a = a;
pHAS.a_abs = a_abs; %[0,     0.133,  0.044,  0.22,  0.012,  0.53,      0.53];
pHAS.rho = rho;
pHAS.randvc=randvc; %zeros(1,length(a));
pHAS.usePhaseCorrection = usePhaseCorrection;  % true/false
pHAS.usePhaseTimeRev = usePhaseTimeRev;  % true/false
pHAS.useRayTracing = useRayTracing;  % true/false
pHAS.numrefl = numrefl;
pHAS.offsetxmm = offsetxmm; % initial values of offset of xducer axis to model axis in mm.
pHAS.offsetymm = offsetymm;
pHAS.c0 = c0;
pHAS.rho0 = rho0;
pHAS.Qzstart = 151;
pHAS.Qzend = 450;
pHAS.Qdmm = pHAS.dmm + (pHAS.Qzstart - 1) * pHAS.Dz;
pHAS.relgeomfocus = [0, 0, 0];
pHAS.padsize = 100;
pHAS.checkPhaseFileValidity = 0;

pHAS.dx = Dx;    %dlp klug - not sure where these should have been defined
pHAS.dy = Dy;
pHAS.dz = Dz;
pHAS.zfocus = zfocus;
pHAS.geom = [0, 0, zfocus];
    pHAS.thrsh_bone = 2000; % Utah: 1000;
    pHAS.thrsh_diploe = 1550; % Utah: 200;
    pHAS.thrsh_brain = 1000; % Utah: -20;
    pHAS.thrsh_fat = 800; % Utah: -200;
%pHAS.Dx = Dx;	

end