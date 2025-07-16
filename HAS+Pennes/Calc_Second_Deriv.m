%*************************************************************
% Calc_Second_Deriv.m
%
% Program to find the perfusion term using
%   temperature information from the cooling portion of the curve
% Based on solving for W in the Pennes bioheat equation and 
%   computing the spatial and time derivatives numerically
%
%
% Nick Todd, University of Utah, January 2009
%
% HOD:  2013-08-16
%       Modified the code to run faster
%       Mainly compressed it all into one calculation, rather than
%       calculating each direction in three steps and then adding all three
%       directions. 
%       An alternative to using this code at all would be to download
%       "DGradient" from Matlab Central and figure out how to compile
%       it...doesn't seem to be working on the Linux server platform?!
%
%*************************************************************

function [d2T] = Calc_Second_Deriv(temps,Vox, Gap)

% % HOD: Old code
% [sx, sy, sz] = size(temps);
% sxD = sx+2; syD = sy+2; szD = sz+2;
% 
% % HOD: Create one and copy and do it all at once, instead of creating each one separately. This is ~10x faster.
% dx_p = zeros(sxD,syD,szD,'single'); 
% dx_n = dx_p; dx_m = dx_p; dy_p = dx_p; dy_n = dx_p; dy_m = dx_p; dz_p = dx_p; dz_n = dx_p; dz_m = dx_p;
% 
% %---------------------------------------------------
% % Do Second spatial Derivative along each direction
% %---------------------------------------------------
% %dx_p = zeros(sxD,syD,szD,'single'); 
% %dx_n = zeros(sxD,syD,szD,'single');
% %dx_m = zeros(sxD,syD,szD,'single');
% dx_m(2:end-1,2:end-1,2:end-1)=temps(:,:,:,1);
% dx_p(2:end-1,3:end,2:end-1)  =temps(:,:,:,1);
% dx_n(2:end-1,1:end-2,2:end-1)=temps(:,:,:,1);
% 
% d2x = (dx_n-2*dx_m+dx_p)/Vox(1)^2;
% 
% %dy_p = zeros(sxD,syD,szD,'single');
% %dy_n = zeros(sxD,syD,szD,'single');
% %dy_m = zeros(sxD,syD,szD,'single');
% dy_m(2:end-1,2:end-1,2:end-1)=temps(:,:,:,1);
% dy_p(3:end,2:end-1,2:end-1)  =temps(:,:,:,1);
% dy_n(1:end-2,2:end-1,2:end-1)=temps(:,:,:,1);
% 
% d2y = (dy_n-2*dy_m+dy_p)/Vox(2)^2;
% 
% %dz_p = zeros(sxD,syD,szD,'single');
% %dz_n = zeros(sxD,syD,szD,'single');
% %dz_m = zeros(sxD,syD,szD,'single');
% dz_m(2:end-1,2:end-1,2:end-1)=temps(:,:,:,1);
% dz_p(2:end-1,2:end-1,3:end)  =temps(:,:,:,1);
% dz_n(2:end-1,2:end-1,1:end-2)=temps(:,:,:,1);
% 
% d2z = (dz_n-2*dz_m+dz_p)/(Vox(3)+Gap)^2;
% 
% d2T = d2x+d2y+d2z;
% 
% d2T = d2T(2:end-1,2:end-1,2:end-1);


%% HOD: New, FASTER, code
[sx, sy, sz] = size(temps);
%sxD = sx+2; syD = sy+2; szD = sz+2;
temps2=zeros(sx+2,sy+2,sz+2,'single');
temps2(2:end-1,2:end-1,2:end-1)=temps;

%temps2=padarray(temps,[1 1 1]); % This seems to be slightly slower (~0.5-1s slower over 1500 iteratons) than the above code

% d2x_2 = (temps2(2:end-1,1:end-2,2:end-1) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(2:end-1,3:end,2:end-1))/Vox(1)^2;
% d2y_2 = (temps2(1:end-2,2:end-1,2:end-1) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(3:end,2:end-1,2:end-1))/Vox(2)^2;
% d2z_2 = (temps2(2:end-1,2:end-1,1:end-2) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(2:end-1,2:end-1,3:end))/(Vox(3)+Gap)^2;
% 
% d2T_2 = d2x_2+d2y_2+d2z_2;

d2T = (temps2(2:end-1,1:end-2,2:end-1) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(2:end-1,3:end,2:end-1))/Vox(1)^2 ...        % x-dir
    + (temps2(1:end-2,2:end-1,2:end-1) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(3:end,2:end-1,2:end-1))/Vox(2)^2 ...        % y-dir
    + (temps2(2:end-1,2:end-1,1:end-2) - 2*temps2(2:end-1,2:end-1,2:end-1) + temps2(2:end-1,2:end-1,3:end))/(Vox(3)+Gap)^2 ;    % z-dir


