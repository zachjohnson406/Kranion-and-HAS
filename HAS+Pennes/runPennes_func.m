function [temps,timevec] = runPennes_func(Q_rel_in,power_acou,timeDuration,pHAS)
% power_acou and timeDuration come from power output for each sonication
%path1= '/v/raid10/users/InsightecExablate/HAS_simulations/Emma/Pennes Bioheat/';
%cd (path1)
disp('running runPennes_func');
Q_rel=Q_rel_in/1000000000; %m3 to mm3
% *** GRAB FROM TREATMENT EXPORT FILE
% TEpower=965.8;
% TEdur=11.25;

Vox = [pHAS.Dx, pHAS.Dy, pHAS.Dz];          % Voxel size
dist_f = [0];           % Voxel spacing
k_mm = 0.527*10^-3;     % Thermal conductivity W/mm/°C
rho_mm = 1045*10^-9;    % Density kg/mm3
sp_h = 3600;            % Specific heat   J/kg/°C   
W = 0;         % Pennes perfusion parameter 9.3 x 10^-9 kg/mm3/s  *** Henrik 15 x 10^-9 kg/m3/s 
% k_mm = pHAS.k_mm;     % Thermal conductivity W/mm/°C
% rho_mm = pHAS.rho_mm;    % Density kg/mm3
% sp_h = pHAS.sp_h;            % Specific heat   J/kg/°C   
% W = pHAS.W;         % Pennes perfusion parameter 9.3 x 10^-9 kg/mm3/s  *** Henrik 15 x 10^-9 kg/m3/s 
Mask_Temps = 'n';       % If we want to just calc temps over some ROI
mask = 0;

%cd (path4)
[ts] = find_dt_size(Vox, k_mm, rho_mm, sp_h, W);    % Find what time step to use in calculation 

% baseLine=round(10/ts);    % in number of time points
% heating=round(TEdur/ts);     % Take from treatment Export
% cooling=round(20/ts);


% totTimePnts=baseLine+heating+cooling;
% power_vec=zeros(1,totTimePnts);
timevec = 1:ts:timeDuration(end);
power_vec=interp1(timeDuration, power_acou, timevec);
%power_vec(1,baseLine+1:baseLine+heating)=TEpower;   % Vector for when the power is on

% Run nummerical model
[temps] = Numerical_Model(Vox,dist_f, k_mm,rho_mm,sp_h,W,power_vec,Q_rel,ts,Mask_Temps,mask);
end

