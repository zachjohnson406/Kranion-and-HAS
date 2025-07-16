
% Code to calculate what temporal step size is needed for numerical
% simulations of Pennes.
% Need to work out the unitis of everything here...looks like if we input
% units according to Nick's Pennes solver this comes out ok, but those
% units are "interesting"...in /mm instead of /m.
%
% [dt] = find_dt_size(Vox, k_mm, rho_mm, sp_h, W)
% 


function [dt] = find_dt_size(Vox, k_mm, rho_mm, sp_h, W)
A=Vox(1)/Vox(2);
B=Vox(1)/Vox(3);
dt=1/(W/rho_mm+2*k_mm*(1+A^2+B^2)/(rho_mm*sp_h*Vox(1)^2)) ;
 
%dt_max=1/(w_max/rho_min+2*k_max*(1+A^2+B^2)/(rho_min*cp_min*dx^2)); 
% Code from Chris on 2017-11-28:
% % Calculate the maximum time step for stability of the thermal model
%     dx=str2num(get(handles.Props_dxEdit,'string'))/1000;
%     dy=str2num(get(handles.Props_dyEdit,'string'))/1000;
%     dz=str2num(get(handles.Props_dzEdit,'string'))/1000;
%     % Dimensionless Parameters
%         A=dx/dy;                       % Dimensionless increment
%         B=dx/dz;                       % Dimensionless increment
%     % Other Parameters
%         w_max=max(handles.Props_w);
%         rho_min=min(handles.Props_rho);
%         c_min=min(handles.Props_cp);
%         k_max=max(handles.Props_k);
%     dt_max=1/(w_max/rho_min+2*k_max*(1+A^2+B^2)/(rho_min*c_min*dx^2)) % Maximum allowable time step before iterations become unstable (s)