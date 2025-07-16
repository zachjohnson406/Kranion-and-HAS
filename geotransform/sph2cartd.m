function [locs] = sph2cartd(ElemLoc, R)
    % a function to convert from spherical coords to cartisian 
    thvect=ElemLoc(:,1); phivect=ElemLoc(:,2);
    for g=1:height(ElemLoc)
      
        y = R*sin(thvect(g)) ;
    
        ss2=R*cos(thvect(g));
        ss3=sin(phivect(g));
    
        x=ss2.*ss3	;				
    
        aa3=cos(phivect(g));
        z=-ss2.*aa3;

        %b=R-y;
        locs(g,:) = [x,y,z];
                                                                        % sign in h.                                                           
    end