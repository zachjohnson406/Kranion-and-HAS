function Modl = prepHAS_Model(Modl, ph, th, Foc_ijk, xtilt)

% Rotate to account for x-tilt of transducer
Modl = rotvolpivrecenterOrig(Modl, [Foc_ijk(1),Foc_ijk(2),Foc_ijk(3)],0,xtilt);


% Rotate to put plate in proper orientation
xCent = floor(size(Modl,1)/2);
yCent = floor(size(Modl,2)/2);
Modl = rotvolpivrecenterOrig(Modl, [xCent,yCent,Foc_ijk(3)],ph, th);

% Reverse x and y axes to match Kranion coordinates
Modl = Modl(end:-1:1,end:-1:1,:);
