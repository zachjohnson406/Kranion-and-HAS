
function mask = reggrow(Modl,xv)
%reggrow - region growing algorithm
%input is 3D matrix, Modl and xv is the 3D location to start from 
[nx,ny,nz] = size(Modl);
if(nz == 1) %2D region grow
    xoff = [-1,1,-nx,nx];
    jxyz1 = sub2ind(size(Modl),xv(1),xv(2));
else
    xoff = [-1,1,-nx,nx,-nx*ny,nx*ny];
    jxyz1 = sub2ind(size(Modl),xv(1),xv(2),xv(3));
end
value = Modl(jxyz1);
group = (zeros(nx*ny,1));
mask = single(0*Modl);
mask(jxyz1) = 1;
group(1) = jxyz1;
jxyzmx = nx*ny*nz;
ng = 1;
while(ng > 0)
    jxyz1 = group(ng);
    ng = ng-1;
    for jj=1:length(xoff)
        jxyzt = jxyz1+xoff(jj);
        if(jxyzt <jxyzmx && jxyzt>0)
            if(mask(jxyzt) == 0)    %not tested yet
                if(Modl(jxyzt) == value)                 %>5)
                    %mask(jxyzt) = 2;
                    ng = ng+1;
                    mask(jxyzt) = 1;
                    group(ng) = jxyzt;
                end
            end
        end
    end
end
    
%     %remove from the list
%     jxyztv = jxyz1 + xoff;
%     tstadd1 = jxyztv(jxyztv < jxyzmx & jxyztv > 0);
%     if(~isempty(tstadd1))
%         jxyztv = tstadd1;
%         tstv = mask(jxyztv)==0 & Modl(jxyztv)<6;
%         tstadd = jxyztv(tstv==1);
%         doneadd = jxyztv(tstv==0);
%         if(~isempty(doneadd))
%             mask(doneadd) = 2;
%         end
%         if(~isempty(tstadd))
%             jxyzadd = tstadd;
%             nadd = length(jxyzadd);
%             mask(jxyzadd) = 1;
%             group(ng+1:ng+nadd) = jxyzadd;
%             ng = ng+nadd;
%         end
%     end
% end
% 
