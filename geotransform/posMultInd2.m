function [PosFinal] = posMultInd2(A,imv)
%UNTITLED4 Summary of this function goes here
%do the 4x4 matrix multiply by a vector of 3D positions
%   Detailed explanation goes C:\Users\zaxpj\Downloads\TestKranionCode.mhere
%Input
% imv (image volume in the shape of the mesh grid to make
% A 4X4 matrix of the transformation to apply

%Create the mesh indices for the multiplication

    [NumRows,NumCols,NumSlices] = size(imv);
    [Pos(:,:,:,2),Pos(:,:,:,1),Pos(:,:,:,3)] = meshgrid(1:NumCols,1:NumRows,1:NumSlices);
    Pos(:,:,:,4) = ones(size(imv(:,:,:,1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AFFINE TRANSFORMATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PosFinal = A*double(reshape(Pos,[size(Pos,1)*size(Pos,2)*size(Pos,3),4])');

tmp1 = Pos(1:20:end,1:20:end,1:20:end,:);
tmp1 = double(reshape(tmp1,[size(tmp1,1)*size(tmp1,2)*size(tmp1,3),4])');

tmp2 = A*tmp1;

% figure(99)
% plot3(tmp1(1,:),tmp1(2,:),tmp1(3,:),'*')
% hold on
% plot3(tmp2(1,:),tmp2(2,:),tmp2(3,:),'o')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% % waitforbuttonpress
end

