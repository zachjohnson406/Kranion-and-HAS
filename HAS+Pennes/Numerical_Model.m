%*************************************************************
% Numerical_Model.m
%
% Program to Caclulate the model temperatures based on the 
%    numerical approach
%
% Nick Todd, University of Utah, January 2009
%
% E.g.
% temps_Mod = Numerical_Model(Vox, dist_f, k_mm, rho_mm, sp_h, W, power_vec, Q_rel, ts, Mask_Temps, mask);
%
%*************************************************************


function [temps_Mod] = Numerical_Model(Vox,dist_f, k_mm,rho_mm,sp_h,...
        W,power_vec,Q_rel,ts,Mask_Temps,mask)

[sx sy sz] = size(Q_rel);
st = length(power_vec);
temps_Mod = zeros(sx,sy,sz,st,'single');

%fprintf('\nModel iterations to go: ')

for ii=1:size(temps_Mod,4)-1
    
    if mod(ii,100)==0    
        disp(['Model iterations to go: ', num2str(size(temps_Mod,4)-ii)]),
        %fprintf('\b');
        %fprintf('%d', num2str(size(temps_Mod,4)-ii));
    end
    
    % Diffusion Term:
    [d2T] = Calc_Second_Deriv(temps_Mod(:,:,:,ii),Vox, dist_f);
    Diff_Term = d2T*k_mm/(rho_mm*sp_h);

    % Perfusion Term:
    Per_Term = -1/rho_mm*W.*temps_Mod(:,:,:,ii);

    % Heat Term:
    Heat_Term = 1/(rho_mm*sp_h)*power_vec(ii)*Q_rel;

    temps_Mod(:,:,:,ii+1) = temps_Mod(:,:,:,ii) + ...
        (Diff_Term + Per_Term + Heat_Term)*ts;

    if (Mask_Temps =='y')
        temps_Mod(:,:,:,ii+1) = temps_Mod(:,:,:,ii+1).*mask;
    end 
    
end

