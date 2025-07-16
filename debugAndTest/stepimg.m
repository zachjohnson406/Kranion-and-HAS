function [ jlast ] = stepimg( img,intscale,dostep )
%stepimg does loop to display
%   Detailed explanation goes here
%DLP 12/16/20 added output of current image displayed
if(exist('dostep','var'))
    if(dostep ~=1) 
        return; 
    end
end
[nx,ny,nz] = size(img);
maxscale = max((img(:)));
minscale = min(img(:));
inscale = 1.0;
minval = minscale;
maxval = maxscale;
if(exist('intscale','var'))
    if(length(intscale) == 2)
        minval = intscale(1);
        maxval = intscale(2);
    else
        inscale = intscale;
        if(intscale < 0)
            inscale = intscale * maxscale;
        end
        minval = minscale/inscale;
        maxval = maxscale/inscale;
    end
end
j=1;
inc = 1;
set(gca,'FontSize',16);
while j <= nz
    imagesc(img(:,:,j),[minval maxval]), axis image, axis off,colorbar;
    jlast = j;
    set(gca,'FontSize',16);
    title(int2str(j),'FontSize',16);
    newj = was('next image',j+inc);
    if(newj == 0 )
        inc = -inc;
    else
        j = newj;
    end
end