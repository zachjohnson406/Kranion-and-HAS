% This function rotates a 3-D segmented volume, such as Modl, so the face of the volume
%   is parallel to the face of the transducer and the beam can be propagated using the
%   hybrid angular spectrum technique. The new volume has values that are picked by nearest
%   neighbor interpolation, except for points that fall outside the range of the old volume
%   (which will be set to outofrangeval). This function can also rotate back to the original
%   orientation a volume of field values, such as pressure or Q (but it doesn't linearly interpolate). 
%   Called by rotatemodl.m and rotateback.m; calls in turn rotcoordpiv.m
%
% function [new_vol] = rotvolpivrecenter(old_vol,pivind,dx,dy,dz,phi,theta,outofrangeval,isrotback)
%    old_vol    = a 3-D matrix of values on the original grid that is to be rotated.
%    new_vol    = a 3-D matrix of values on a grid that has been rotated.
%    pivind      = a 3x1 vector representing the INDEX point about which the rotation
%                 will be performed; order: [colpiv, rowpiv, pagepiv] or (x,y,z).
%    dx, dy, dz = physical spacing between the x, y, and z values respectively (usu. in m or mm units). 
%                 These are needed to keep the correct aspect ratio of the volume.
%    phi      = rotation about the y axis of the old_volume (in radians).
%    theta        = rotation about an intermediary x-axis (in radians).
%    outofrangeval	 = either 1(representing water) for Modl, or 0 for a field array.
%    isrotback  = flag that tells whether modl is being rotated back (1).


% Last updated by:
%   Doug Christensen 4/27/11 - Changed Scott Almquist's order of rotation to match coordinate system
%       of Modl: y (vertical) is in direction of increasing ROW index, x (horiz) is
%       in direction of increasing COLUMN index.  Thus phi now rotates first
%       around the y or row axis, and theta next rotates around the rotated x or column axis,
%       leaving the rotated x-axis in the x-z plane. Changed order and names of phi & theta in input
%       argument list. Also the pivot point vector is in the order [colpiv, rowpiv, pagepiv].
%   This function can be used to rotate forward a 3D array of media integers, as in a Modl, as well
%       as to rotate back to the original orientation an array of double precision field values
%       such as pressure or Q (in which case the flag isrotback is set to 1.)
%   Added outofrangeval (usually either 1 (=water) if volume being rotated is integer Modl, or 0
%       if volume is a field of pressure or Q values) to the input argument list. 
%       Changed the nested for loops (very slow for a 301x301x360 Modl) to array multiplications.
%   
%   Added shifts in xaxis and yaxis to place pivot point (which is also forced to be the geometric
%       focus) to be at x,y CENTER of new volume (this is done to allow beam
%       from transducer to be more centered on model and avoid clipping by edges of model). In turn,
%       this required that the geometric focus in gui be moved by the same amount; this was
%       accomplished by changes to rotatemodl.m and rotateback.m.
%       isrotback flag added to arguments to tell that volume is being rotated back, not forward.
%
%   4/15/17 - This function may NOT be needed for rotating BACK any volume if the Modl is not rotated
%       back, but rather replaced by OrigModl (see rotateback.m), and both Q and wdisp are rotated with
%       interpolation by rotvolpivrecenterinterp.m for improved resolution, not this function.
%
% Copyright D.A.Christensen 4/15/17
%
%   July 2020
%       Blake Zimmerman solved issue ("iss13" in GitHub). Replaced linear indexing with 3D
%       interpolation.
%
%   Feb 2021
%       Michelle Kline and Nick Richards solved iss18 in GitHub. After translation of pivot point
%       to 0,0,0, then rotation of the model, HAS expects the model to be recentered in x,y on the
%       pivot point (which is forced to be the geometric focus). Code was translating pivot point
%       back to original position, not to x,y center.

function [new_vol] = rotvolpivrecenterOrig(old_vol,pivind,phi,theta)

% Flip the theta and phi - This rotates the model, but we want that to be
% the inverse of how we are rotating the transducer. 
% theta = -theta;
% phi = -phi;

% Create the reference for the coordinate system
Rdefault = imref3d(size(old_vol));

% Create the rotation matrix about Y' axis by phi
yRotation = [cosd(phi) 0 sind(phi) 0; 
             0 1 0 0; 
             -sind(phi) 0 cosd(phi) 0;
             0 0 0 1];

% Create the rotation matrix about X' axis by theta
xRotation = [1 0 0 0;
             0 cosd(theta) sind(theta) 0;
             0 -sind(theta) cosd(theta) 0;
             0 0 0 1];

% Create the translation vector to make 0, 0, 0 be the center
tX = pivind(1);
tY = pivind(2);
tZ = pivind(3);
tTranslationToCenterAtOrigin = [1 0 0 0; 
                                0 1 0 0; 
                                0 0 1 0;
                                -tX -tY -tZ 1];

% Define center of model (for recentering on pivot point)
sizeDim = size(old_vol);
centerPt(1) = round(sizeDim(2)/2);
centerPt(2) = round(sizeDim(1)/2);
centerPt(3) = round(sizeDim(3)/2);
                            
% Create the translation vector to recenter the model in the Y/Z plane (move the pivot point to the center)
tTranslationBackToOriginalCenter = [1 0 0 0; 
                                    0 1 0 0; 
                                    0 0 1 0;
                                    centerPt(1) centerPt(2) tZ 1];  % iss18, was tX tY tZ 1
                                
% Create the affine object - it matters that x comes first
tformCenteredRotation = tTranslationToCenterAtOrigin*yRotation*xRotation*tTranslationBackToOriginalCenter;
tformCenteredRotation = affine3d(tformCenteredRotation);

% Now loop over the unique values in the array
unique_tissues = unique(old_vol(:));
% tissue_volumes = zeros([size(old_vol), length(unique_tissues)]);
current_max = zeros(size(old_vol));
new_vol = ones(size(old_vol));


for i=1:length(unique_tissues)
    temp = single(old_vol==unique_tissues(i));
    [temp, ~] = imwarp(temp, Rdefault, tformCenteredRotation, 'linear', 'OutputView', Rdefault);
    mask = temp > current_max;
    current_max(mask) = temp(mask);
    new_vol(mask) = unique_tissues(i);
end
