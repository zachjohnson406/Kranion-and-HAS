function [] = imgWithTemp(img1,img2,orientation)
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
        minscale = 0;

        %plot first data 
        % Create axes for the background image
        ax1 = axes;
        im = imagesc(ax1, img1(:,:,:));
        im.AlphaData = 1; % Set transparency for background image
        axis image;
        colormap(ax1, 'gray'); % Set colormap for background to grayscale
        
        hold on;
        
        % Create axes for the foreground image
        ax2 = axes;
        im1 = imagesc(ax2, img2(:,:,:), [minscale ceil(maxscale)]);
        alphaData = ones(size(img2(:,:,:))) * 1; % Initialize transparency data
        alphaData(img2(:,:,:) >= 0 & img2(:,:,:) <= minscale) = 0; % Adjust as needed
        im1.AlphaData = 0.5;
        
        axis image;
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        axis off; 
        
        % Link axes for synchronized zooming
        linkaxes([ax1, ax2]);
        
        % Optionally, adjust the position of the axes if necessary
        set([ax1, ax2], 'Position', [0.1 0.1 0.8 0.8]);
        
        % Colormap adjustments (if needed, though not specified in the question)
        % colormap(ax1, 'gray');
        colormap(ax2, 'hot');
        %set the axes and colorbar position 
        set([ax1,ax2],'Position',[.07 .11 .685 .815]); 
        cb1 = colorbar(ax2,'Position',[.75 .11 .0675 .815]); 
        %caxis([37 58]);
       
        ylabel(cb1, 'Temperature ($^\circ$C)', 'Interpreter', 'latex');

        xlabel(ax1, 'L (mm)');
        ylabel(ax1, 'P (mm)');
        ax1.XTick = [50 100 150 200]; % 4 ticks evenly spaced across x-axis
        ax1.YTick = [50 100 150 200];
        title('Simulation');
end
%%