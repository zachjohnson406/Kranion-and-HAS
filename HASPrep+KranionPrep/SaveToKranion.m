% Interpolates 3d array from HAS coordinate space to Kranion coordinate
% space. 
%
% @OUTPUTS
%   out: Interpolated array. 
%
% @INPUTS
%   array: Array to interpolate. 
%   Foc: Natural focus of transducer in model. 
%   Xtilt: Tilt of transducer along x-axis (degrees).
%   unTransform: Matrix going from Kranion voxel space to isotropic voxel
%   space. 
%   dz: Size of HAS model in the third dimension.    
%
%
% Zach Johnson
% June 28th, 2024


function [out] = SaveToKranion(array,Foc,Xtilt,unTransform,dz,export,toKranion, arrayType)        
   
    %variables to be made parameters
    padsize = 0; 
    %scale = 2.2069; 
    Foc = Foc + padsize; 

    [xdnx,xdny,xdnz] = size(array);
    xcent = ceil((xdnx)/2); ycent = ceil((xdny)/2); zcent = ceil((xdnz)/2);

    % Reverse x and y axes to match Kranion coordinates
    array = array(end:-1:1,end:-1:1,:);

    %get rid of Xtilt
    array=rotvolpivrecenterinterp(array,[xcent,ycent,Foc(3)],1,1,1,0,-Xtilt,1,1);

    %recenter
    array=rotvolpivrecenterinterp(array,[(xcent-(xcent-Foc(1))),ycent-(ycent-Foc(2)),ceil(zcent)],1,1,1,0,0,1,1);

    % padding
    array = array(1+padsize:end-padsize, 1+padsize:end-padsize, 1+padsize:end-padsize);

    array = array(1:end-1, 1:end-1, 1:end);


    %original CT was not isotropic 
    isotropic_scale = unTransform^(-1);


    out_xdim = 512;
    out_ydim = 512;
    [sim_xdim, sim_ydim, sim_zdim, sonication] = size(array);
    
    % Define the affine transformation matrix
    simToCT = [out_xdim/sim_xdim 0 0 0;
               0 out_ydim/sim_ydim 0 0;
               0 0 1 0;
               0 0 0 1];
    
    % Create the affine3d object
    affine = affine3d(simToCT);
    
    % Apply the affine transformation with nearest-neighbor interpolation
    out = imwarp(abs(array), affine, 'InterpolationMethod', 'nearest');
    %[out] = TriLinInterpV3(array,(isotropic_scale),zeros(512,512,dz)); 

    % row-major vs. column-major
    out = permute(out,[2 1 3]);
   
    % save image information to nifti. Set nifti header information.
    info = struct;
    voxelDimensions = [isotropic_scale(1,1), isotropic_scale(2,2), isotropic_scale(3,3)];
    info.ImageSize = size(out);
    info.PixelDimensions = voxelDimensions;
    info.Datatype = 'single';
    info.Description = 'HAS Pressure Magnitude';
    info.Version  = 'NIfTI1';
    info.BitsPerPixel= 64;
    info.SpaceUnits= 'Unknown';
    info.TimeUnits= 'None';
    info.AdditiveOffset= 0;
    info.MultiplicativeScaling= 0;
    info.TimeOffset= 0;
    info.SliceCode= 'Unknown';
    info.FrequencyDimension= 0;
    info.PhaseDimension= 0;
    info.SpatialDimension= 2;
    info.DisplayIntensityRange= [0 0];
    info.TransformName= 'Sform';
    toKranion(1,1) = toKranion(1,1) *  -1;
    toKranion(1,2) = toKranion(1,2) *  -1;
    toKranion(1,3) = toKranion(1,3) *  -1;
    toKranion(1,4) = toKranion(1,4) *  -1;


    toKranion(2,1) = toKranion(2,1) *  -1;
    toKranion(2,2) = toKranion(2,2) *  -1;
    toKranion(2,3) = toKranion(2,3) *  -1;
    toKranion(2,4) = toKranion(2,4) *  -1;
    affineMatrix = toKranion'; 

    info.Transform.T = affineMatrix;
    info.Qfactor = 1;
    info.raw= struct;
    niftiwrite(abs(out), [arrayType, export(1:end-4),'.nii'], info);

    