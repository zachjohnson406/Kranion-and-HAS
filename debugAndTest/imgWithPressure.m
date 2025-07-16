function [] = imgWithPressure(img1, img2, orientation, ax, voxel_size)
    if nargin < 4
        ax = gca; % Use current axes if none provided
    end
    if nargin < 5
        voxel_size = [1, 1]; % Default voxel size [dx, dy] in mm
    end
    
    if exist('orientation','var') && ~isempty(orientation)
        switch lower(orientation)
            case 'coronal'
                img1 = permute(img1, [3,2,1]);
                img2 = permute(img2, [3,2,1]);
            case 'axial'
                img1 = permute(img1, [1,2,3]);
                img2 = permute(img2, [1,2,3]);
            case 'sagittal'
                img1 = permute(img1, [1,3,2]);
                img2 = permute(img2, [1,3,2]);
            otherwise
                error('Orientation must be "coronal", "axial", or "sagittal".');
        end
    end
    
    % Get slice dimensions
    [ny, nx] = size(img1); % Note: changed order to match image convention
    maxscale = max(img2(:));
    minscale = 0;
    
    % Define spatial resolution
    dx = voxel_size(1); % mm per pixel in x direction
    dy = voxel_size(2); % mm per pixel in y direction
    
    % Create coordinate vectors centered at origin
    x_coords = (1:nx) * dx - (nx/2 + 0.5) * dx;  % Center at 0
    y_coords = (1:ny) * dy - (ny/2 + 0.5) * dy;  % Center at 0
    
    % Convert grayscale background to RGB
    img1_rgb = repmat(mat2gray(img1), [1 1 3]); % scale to [0 1] and replicate to RGB
    
    % Plot RGB background image with proper coordinates
    image(ax, x_coords, y_coords, img1_rgb);
    hold(ax, 'on');
    
    % Overlay pressure data with colormap using proper coordinates
    hImg = imagesc(ax, x_coords, y_coords, img2, [minscale maxscale]);
    alphaData = ones(size(img2));
    alphaData(img2 <= 0) = 0; % transparent where pressure <= 0
    hImg.AlphaData = 0.7;
    colormap(ax, 'parula'); % only for the overlay
    
    % Set axis properties
    axis(ax, 'image');
    ylabel(ax, 'P/A (mm)', 'FontSize', 13);
    xlabel(ax, 'R/L (mm)', 'FontSize', 13);
    
    % Set symmetric tick marks
    max_extent_x = max(abs(x_coords));
    max_extent_y = max(abs(y_coords));
    
    % Create symmetric tick marks (adjust spacing as needed)
    tick_spacing = 50; % mm - adjust this value as needed
    x_ticks = [];
    y_ticks = [];
    
    % Generate tick marks that fit within the image bounds
    for i = -10:10 % Generous range
        x_tick = i * tick_spacing;
        y_tick = i * tick_spacing;
        if x_tick >= min(x_coords) && x_tick <= max(x_coords)
            x_ticks = [x_ticks, x_tick];
        end
        if y_tick >= min(y_coords) && y_tick <= max(y_coords)
            y_ticks = [y_ticks, y_tick];
        end
    end
    
    % Ensure we always include 0
    if ~ismember(0, x_ticks)
        x_ticks = sort([x_ticks, 0]);
    end
    if ~ismember(0, y_ticks)
        y_ticks = sort([y_ticks, 0]);
    end
    
    ax.XTick = x_ticks;
    ax.YTick = y_ticks;
    

    
    % Add colorbar
    cb = colorbar(ax);
    cb.Label.String = 'kPa';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 13;

end