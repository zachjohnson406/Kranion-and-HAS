function [] = stepimgWithTemp(img1,img2,orientation)
    if(exist('orientation','var'))
       if strcmp(orientation,'coronal')
            img1 = permute(img1,[3,2,1]);
            img2 = permute(img2,[3,2,1]);
       elseif strcmp(orientation,'axial')
           img1 = permute(img1,[1,2,3]);
           img2 = permute(img2,[1,2,3]); 
       elseif strcmp(orientation, 'sagittal')
           img1 = permute(img1,[1,3,2]);
           img2 = permute(img2,[1,3,2]);
       else
           print('Must enter "coronal", "axial", or "sagittal"');
       end
    end
    [~,~,z] = size(img1);
    [nx,ny,nz] = size(img1);
    maxscale = max((img2(:)));
    minscale = 3;
    j=1;
    inc = 1;
    while j <= nz
        %plot first data 
        % Create axes for the background image
        ax1 = axes;
        im = imagesc(ax1, img1(:,:,j));
        im.AlphaData = 1; % Set transparency for background image
        axis square;
        colormap(ax1, 'gray'); % Set colormap for background to grayscale
        ax1.Visible = 'off';
        ax1.XTick = [];
        ax1.YTick = [];
        
        hold on;
        
        % Create axes for the foreground image
        ax2 = axes;
        im1 = imagesc(ax2, img2(:,:,j), [minscale maxscale]);
        alphaData = ones(size(img2(:,:,j))) * 1; % Initialize transparency data
        alphaData(img2(:,:,j) >= 0 & img2(:,:,j) <= minscale) = 0; % Adjust as needed
        im1.AlphaData = alphaData;
        
        axis square;
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        
        % Link axes for synchronized zooming
        linkaxes([ax1, ax2]);
        
        % Optionally, adjust the position of the axes if necessary
        % set([ax1, ax2], 'Position', [0.1 0.1 0.8 0.8]);
        
        % Colormap adjustments (if needed, though not specified in the question)
        % colormap(ax1, 'gray');
        colormap(ax2, 'hot');
        %set the axes and colorbar position 
        set([ax1,ax2],'Position',[.07 .11 .685 .815]); 
        cb1 = colorbar(ax2,'Position',[.75 .11 .0675 .815]); 
        ylabel(cb1, '$^\circ$C', 'Interpreter', 'latex', 'FontSize', 16);
        jlast = j;
        newj = was('next image',j+inc);
        if(newj == 0 )
            inc = -inc;
        else
            j = newj;
        end
    end
end
