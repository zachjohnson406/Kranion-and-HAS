function [outarray] = TriLinInterpV3(imageFrom,Afrom_to,imageTo)
%TriLinInterpV2 Summary of this function goes here
%   Use trilinear interpolation in image to determine value for each point
%   in outarray
%Input
% imageFrom     = image volume to interpolate from (sending image)
%            	 the array from which values are to be obtained
% Afrom_to      =   Affine transformation between 
% imageTo  Mapping of indices from receiving image (outarray)
%                to a floating point position (x, y, z) in image, the 
%                image to be interpolated
%Output
% outarray       the output array (interpolation result)
%

%TargToImageIndices has dimensions of the receiving image 
    Pos =  posMultInd2(Afrom_to,imageTo);
    Pos = Pos';
    %dctinct = TriLinInterpDLP(dynactimv,CTtoDCTindices);
[sx,sy,sz] = size(imageTo);       
sq = 4;
outarray = single(zeros(sx,sy,sz));

[nx,ny,nz] = size(imageFrom);
% 
% convert to a vector of addresses
% Pos = reshape(TargToImageIndices,[sx*sy*sz,sq]);
clear TargToImageIndices;
%Limit the position to the image boundaries
mask0 = Pos > 0;
mask1 = Pos(:,1) < nx;
mask2 = Pos(:,2) < ny;
mask3 = Pos(:,3) < nz;
mask = mask1.*mask2.*mask3 .* mask0(:,1) .*mask0(:,2) .* mask0(:,3);
Pos(:,1) = mask .* Pos(:,1);
Pos(:,2) = mask .* Pos(:,2);
Pos(:,3) = mask .* Pos(:,3);

clear('mask0','mask1','mask2','mask3','mask');
iPos = fix(Pos);
frac = Pos - iPos;      %a fraction between 0 and 1 for each position direction
clear Pos;
ctiadd = uint32(iPos(:,1)) + uint32(nx)*uint32(iPos(:,2)-1) + uint32(nx*ny)*uint32(iPos(:,3)-1);
% ctiadd = sub2ind(size(imageTo),iPos(:,1),iPos(:,2),iPos(:,3));

ctiadd = uint32(ctiadd);
ctiadd = max(1,ctiadd);
ctiadd = min((uint32(nx*ny*nz)),ctiadd);
% 
clear iPos;
 f1 = frac(:,1);
 f2 = frac(:,2);
 f3 = frac(:,3);
 clear frac;
outarray(:) =         (1-f1) .* (1-f2) .* (1-f3)   .* imageFrom(ctiadd);
outarray(:) = outarray(:)+    f1  .* (1-f2) .* (1-f3)   .* imageFrom(ctiadd + 1);
outarray(:) = outarray(:)+ (1-f1) .* (  f2) .* (1-f3)   .* imageFrom(ctiadd +     nx);
outarray(:) = outarray(:)+(   f1) .* (  f2) .* (1-f3)   .* imageFrom(ctiadd + 1 + nx);
outarray(:) = outarray(:)+ (1-f1) .* (1-f2) .* f3       .* imageFrom(ctiadd +        nx*ny);
outarray(:) = outarray(:)+ (  f1) .* (1-f2) .* f3       .* imageFrom(ctiadd + 1 +    nx*ny);
outarray(:) = outarray(:)+ (1-f1) .* (  f2) .* f3       .* imageFrom(ctiadd +     nx +nx*ny);
outarray(:) = outarray(:)+    f1  .*    f2  .* f3       .* imageFrom(ctiadd + 1 + nx +nx*ny);

end

