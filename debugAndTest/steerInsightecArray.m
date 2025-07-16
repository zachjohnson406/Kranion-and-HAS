% Steer Array
% 
% This function steers an ultrasound transducer array in a homogeneous
% medium
% 
% @INPUTS
%   loc: x,y,z location of elements (in m)
%   freq: frequency (in MHz)
%   c: speed of sound (in m/s)
%   desiredFocus (in m)
% 
% @OUTPUTS
%   ph: vector of phases corresponding to each element in loc (in radians)

function ph = steerInsightecArray(loc,freq,c,desiredFocus)

% Convert everything to standard units
freq = freq*1e6; % from MHz to Hz

% Find the distance between each element and the desired Focus
dist = sqrt((-loc(:,1)-desiredFocus(1)).^2+...
    (loc(:,2)-desiredFocus(2)).^2+...
    (loc(:,3)-desiredFocus(3)).^2);

% Compute the phase for each element
t = dist/c;
t = t-t(1);
ph = 2*pi*freq*t;