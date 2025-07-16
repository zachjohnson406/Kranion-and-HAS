function [thvect, phivect] = cart2sphd(x,y,z, R)
    % Goes from cartesian coordinates to, sperical coordinates (as defined
    % by Doug Chritensen for HAS). 
    %
    % Zach Johnson

    R = sqrt(x.^2+y.^2+z.^2);

    thvect = asin(y./R);
    den = R .* cos(thvect);

    % Force elements to be on sphere: This needs to change to account for
    % small deviations from the spherical surface
    arg = x./den;
    % arg(arg>1) = 1;
    % arg(arg<-1) = -1;

    phivect = asin(arg);

    if max(imag(phivect))
        keyboard
    end