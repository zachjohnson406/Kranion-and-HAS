%% Function to find the hottest voxel
%
% Input:
% temps - 4D temp data set in format (x, y, z, t)
% region - Sub-region of data set to look at e.g.:
%                                                   [60 75 40 54 5  7  32 44] 
%                                                    x1 x2 y1 y2 z1 z2 t1 t2
% printRes - Volontary input. If set to 'n', the results won't be printed.
%
% Output:
% x, y, z, t - Coordinates of hottest voxel
% maxT = maximum temperature found
%
%
% Example:
%           [x y z t maxT] = findHottestVoxel(temps, [60 75 40 54 5 7 32 44])
%
% Henrik Odeen
% Univ. of Utah
% June 2012




function [x, y, z ,t, max_val] = findHottestVoxel(temps, region, printRes)

x1 = region(1); x2 = region(2);
y1 = region(3); y2 = region(4);
z1 = region(5); z2 = region(6);
t1 = region(7); t2 = region(8);


mask=zeros(size(temps));
mask(x1:x2,y1:y2,z1:z2,t1:t2)=1;
temps2=temps.*mask;

% Find the voxel with max temperature in temps, and its position
[max_val, position] = max(temps2(:)); 

% Transform the index to 4 indices, given the size of temps
[x,y,z,t] = ind2sub(size(temps),position);

if nargin==3 && printRes=='n'
    % Do nothing
else
    disp(['x=' ,int2str(x), ', y=' ,int2str(y), ', z=' ,int2str(z), ', t=' ,int2str(t), ', and maxT=' ,num2str(max_val), ' C'])
end
