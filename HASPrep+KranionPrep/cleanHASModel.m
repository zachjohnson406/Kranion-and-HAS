function [new_Modl] = cleanHASModel(Modl,stop)
    % Removes partial volume voxels in the brain of Model. Changes brain that is outside of the skull to
    % fat. Removes CT holder. 
    % 
    % Zach Johnson
    [nx,ny,nz] = size(Modl);
    start_page = ceil(nz/2); 
    new_Modl = Modl(:,:,:); 
    mskBrain = zeros(size(new_Modl));
    mskBrain(new_Modl == 5) = 1;
    mskBrain = imfill(mskBrain,'holes');

    mskFat = zeros(size(new_Modl));
    mskFat(new_Modl == 3) = 1;
    mskFat = imfill(mskFat,'holes');

    mskSkull = zeros(size(new_Modl));
    mskSkull(new_Modl >= 6) = 1;
    mskSkull = imfill(mskSkull,'holes');
    for z = 1:size(mskSkull, 3)
        mskSkull(:, :, z) = imfill(mskSkull(:, :, z), 'holes');
    end

    new_Modl(mskBrain == 1) = 5;
    [nx, ny, nz] = size(mskBrain);
    start_row = ceil(nx / 2);
    start_col = ceil(ny / 2);
    
    connected_array = zeros(size(mskBrain));
    
    for z = 1:stop +10
        slice = mskBrain(:, :, z);
        CC = bwconncomp(slice, 4);
        start_idx = sub2ind([nx, ny], start_row, start_col);
        for j = 1:CC.NumObjects
            if ismember(start_idx, CC.PixelIdxList{j})
                connected_array(:, :, z) = false;
                connected_array(CC.PixelIdxList{j} + (z - 1) * nx * ny) = 1;
                break;
            end
        end
    end
    
    new_Modl((connected_array == 0) & (new_Modl == 5)) = 1;
    new_Modl(mskSkull == 1 & new_Modl ~= 6 & new_Modl ~= 5 & new_Modl ~= 3) = 7;
    new_Modl(Modl ~= 1 & new_Modl == 1) = 3;
    new_Modl(:,:,(stop+10):end) = Modl(:,:,(stop+10):end);
    new_Modl = removeCTholderZach(new_Modl);

    mskSkullFat = zeros(size(new_Modl));
    mskSkullFat(new_Modl >= 6 | new_Modl == 3) = 1;
    for z = 1:size(mskSkullFat, 3)
        mskSkullFat(:, :, z) = imfill(mskSkullFat(:, :, z), 'holes');
    end
    new_Modl(new_Modl == 7 & Modl == 1) = 1;


end
