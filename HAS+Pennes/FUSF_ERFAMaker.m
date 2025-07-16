%  This is a modified version of D.A. Christiensen's ERFA maker. 
%  Creates 7 efra planes, one for each section of
%  a hemispherical trasducer. Geometry loaded from Kranion. 
%
% @INPUTS
%   Seg: A stucture containing transducer element locations divied into 7
%   sections. 
%   TH_rotation: 7 rotation angles about the X axis to bring each of the 7
%   transducer segments to align with the axis of propagation. 
%   PHI_rotation: 7 rotation angles about the Y axis to bring each of the
%   7 transducer segments to align with the axis of propagation. 
%   fullarray: A maxrix of all trasducer element locations.
%   fullfilepath: Path to the Kranion export (.krx) file in use. 
%
%   Zach Johnson and Taylor Webb 
%   June 5, 2024
    
function FUSF_ERFAMaker(Seg,TH_rotation,PHI_rotation,fullarray,fullfilepath)
        %       Set parameters for making the ERFA for the InSightec array
        %       Later this can be a swtiched section based on the transducer geometry

        %       Parameter file for ERFA calculations of fields on an ERFA plane a distance sx
        %		away from a spherical transducer.  Used by program ERFAMaker4.
        %       The size and resolution of the ERFA plane should be appropriate for the Modl (see below).
        %       
        %
        %	    Copyright D.A. Christensen 2017.
        %	    Jan 20, 2017.

        %       Best to think of transducer axis aimed horizontally, as opposed to
        %       vertically

        % fMHz=.940;	% frequency (in MHz).
        % 
        % Dv=0.1538;	% max ARC length of source in theta (vert) direction (in m).
        % Dh=0.1538;	% max ARC length of source in phi (horiz) direction (in m).

        % imax=88;	% 201 number of theta (elevation or vert) angle increments to source.
        % kmax=88;	% 201 number of phi (azimuth or horiz) angle increments to source.
                % Note: imax and kmax should be ODD to keep symmetry in planes. The spacing referred
                % 	to the transducer surface should be smaller than 1/2 wavelength.
    
        % R=0.10;	% radius of curvature of transducer (in m).
                % Note: R can be large (to simulate planar transducers) but not inf.
    
        isPA=1; % this is a pamhased array.
    
        % relem=0.002;   % radius of circular element of array transducer (in m) (if applicable).
        % activearea=0.23; % ratio of active (radiating) area to total area of array (if applicable).
    
        % d=0.05;	% distance from furthest point on transducer curve to ERFA plane (in m).
        c0=1.5e3;	% speed of sound in medium between transducer and plane, e.g., water (in m/s).
        rho0=1e3;	% density of medium between transducer and plane, e.g., water (in kg/m^3).
    
        % Lv=0.10;	% length of ERFA plane in vert direction (in m).
        % Lh=0.10;	% length of ERFA plane in horiz direction (in m).
                % Note: The generalized ERFA technique assumes a square ERFA plane, so Lv=Lh.  Also, to
                %   allow for offset between the xducer and the Modl, the ERFA plane should be somewhat
                %   larger than the transducer or Modl (maybe twice the size of the model).
    
        % lmax=401;	% number of space increments in vert direction in ERFA plane; should be ODD.
        % mmax=401;	% number of space increments in horiz direction in ERFA plane; should be ODD.
                % Note: The incremental spacing in the ERFA plane should not be significantly larger than the
                %   spacing of the Modl so later interpolation is accurate.
    
        fMHz = 0.6700;
        Dv = 0.22;
        Dh=Dv;
        imax = 37;
        kmax = imax;
        R = 0.1500;
        relem = 0.0055;
        d = 0.0532;
        Lv = 0.1882;
        Lh = Lv;
        lmax = 401;
        mmax = lmax;
    
        dTh = asin(relem/R);
        Th_total = 2*asin((Dv/2)/R);
        imax = 3*ceil(Th_total/dTh);
        if ~mod(imax,2)
            imax = imax+1;
        end
        kmax = imax;


        
        %% Launch ERFA Maker for each array section
        %Run ERFA for each segment, rotating each by:
        % [TH_centers_rot, PHI_centers_rot]

        disp('No precomputed ERFA file found...computing ERFA now.')
        Dim=[Dv,Dh]; anginc=[imax,kmax]; Len=[Lv,Lh]; planinc=[lmax,mmax]; % combine some params.

        prompt={'Frequency of transducer (in MHz):',...
            'Overall ARC dimensions of transducer (vert, horiz) in m:',...
            'Number of angle increments to transducer (vert, horiz):',...
            'Radius of curvature of transducer face (in m):',...
            'Is this a phased-array transducer? (1=PA, 0=solid)',...
            ' >Radius OR 1/2 width/height of array elements (in m):',...
            ' >Reserved for future use',...
            'Distance from back of xducer to ERFA plane (in m):',...
            'Size of ERFA plane area (vert, horiz in m):',...
            'Number of increments in ERFA plane (vert, horiz):',''};
        commentstr='MAKE CHANGES if any, TO VALUES, THEN PRESS OK';
        titl='Parameters for ERFA Calculation'; lines=1;
        initans={num2str(fMHz),num2str(Dim),num2str(anginc),num2str(R),...
            num2str(isPA),num2str(relem),'0',num2str(d),...
            num2str(Len),num2str(planinc),commentstr};
        params=inputdlg(prompt,titl,lines,initans);
        if isempty(params); return; end

        fMHz=str2double(char(params(1)));
        Dim=str2num(char(params(2))); Dv=Dim(1); Dh=Dim(2);
            if Dv~=Dh; hdb=warndlg(['The transducer outline is NOT CIRCULAR. The Solid tranducer will be changed to '...
                'Circular. Continue if that is okay.'], 'Transducer Not Circular', 'modal'); uiwait(hdb); end 
        anginc=str2num(char(params(3))); imax=anginc(1); kmax=anginc(2);
        R=str2double(char(params(4)));
        isPA=str2double(char(params(5)));
        relem=str2double(char(params(6)));
        d=str2double(char(params(8)));
        Len=str2num(char(params(9))); Lv=Len(1); Lh=Len(2);
            if Lv~=Lh; hdb=warndlg('The ERFA plane is NOT SQUARE. Continue only if that is okay.', 'ERFA Not Square',...
                 'modal'); uiwait(hdb); end
        planinc=str2num(char(params(10))); lmax=planinc(1); mmax=planinc(2);
            if 2*round(lmax/2)==lmax || 2*round(mmax/2)==mmax;
              hdb=warndlg('The number of ERFA plane increments should be ODD integers. Continue only if even is okay.',...
                 'Not Odd', 'modal'); uiwait(hdb); end

        f=fMHz*1e6;	% convert to Hz.

        i=1:imax; k=1:kmax; l=1:lmax; m=1:mmax;     % set up indices.  Careful: i used as index here, j is imag number.
        dth=single(Dv/(R*imax)); dphi=single(Dh/(R*kmax));		% incremental size of source angle (in radians).
        dyp=single(Lv/(lmax-1)); dxp=single(Lh/(mmax-1));	% incremental size of steps in ERFA plane (in m).
        th1=dth*(i-round(imax/2));      % angle row vector, centered; imax and kmax should be odd for symmetry.
        phi1=dphi*(k-round(kmax/2));	% angle row vector, centered.
        thmesh1=repmat(th1',1,kmax); % imax x kmax matrix of theta values, 'meshgrid' style, over entire rectangular
                                        % angular area that encompasses xducer.
        phimesh1=repmat(phi1,imax,1);% imax x kmax matrix of phi values, 'meshgrid' style, over entire xducer area.

        Zm=rho0*c0;     % impedance of medium that waves radiate into; rho0 and c0 read in the parameter file.

        button1=questdlg('Is this a phased array or solid transducer?','TRANSDUCER FORMAT','Phased Array',...
            'Solid','Cancel','Cancel');

        switch button1  %  could also switch on isPA, but this is more flexible.
          case 'Cancel'; return
           % ----------------    
          case 'Solid';     % The next steps model a solid CIRCULAR transducer:
            hs=waitbar(0,'Evaluating Rayleigh-Sommerfeld integral...');
            Arc=min(Dv,Dh);   % diameter along arc of circular transducer. NOTE: Only valid for circular transducers.
            distfromcent=R*acos(cos(phimesh1).*cos(thmesh1));   % arc distance from center; see notes.
            pointmap=zeros(imax,kmax,'single');     % initialize map of valid points inside xducer circle.
            pointmap(distfromcent <= Arc/2) = 1;    % it's inside the xducer radius, so set point = 1; otherwise 0.
            cth=cos(thmesh1);					%  important term: cos(theta) matrix needed for spherical integrations.
            XducerArea=R*R*dth*dphi*sum(sum(pointmap.*cth));     % add up the valid point areas to get overall area.
            ptunif=single(sqrt(2*Zm*1/XducerArea));	% peak pressure for 1 W of total radiated power, if uniform and solid.
            pt=pointmap*ptunif;       % now set up the transducer pressure matrix.

            %--- Rayleigh-Summerfeld integral done next with three implicit for loops ---
            ss5=R*(cth.*sin(phimesh1));         % 2D matrix: imax x kmax.
            s=repmat(ss5,[1,1,lmax]);           % now 3D array.
            xp=dxp*(m-round(mmax/2));           % vector of x points along (with vers. 8) the horizontal axis.

            aa5=R*(cth.*cos(phimesh1));
            a=repmat(aa5,[1,1,lmax]);           % 3D array: imax x kmax x lmax.
            b=R-a;  clear a         % to save memory.

            tt= R*sin(thmesh1); 
            t=repmat(tt,[1,1,lmax]);            % 3D array of t.
            yyp(1,1,:)=dyp*(l-round(lmax/2));	% turn y points vector into a 'page' vector, lmax pages long.
            yp=repmat(yyp,[imax,kmax,1]);		% make 3D array.
            term1=(t-yp).^2;	clear t yp        % to save memory.		
            term3=(d-b).^2;		clear b          % to save memory.	

            ppc=pt.*cth;                        % pt is circ xducer source matrix; cth is cos for spherical integration.
            pc=repmat(ppc,[1,1,lmax]);          % product of pressure and cos theta now 3D array.
            kk=2*pi*f/c0;
            ppi=zeros(kmax,lmax,mmax,'single');          % preallocate storage for ppi.

            for mi=1:mmax
                % Rayleigh-Sommerfeld integral explicit over xp only here:
                r=sqrt(term1 + term3 + (s-xp(mi)).^2);	% r is 3D array for each value of xp.
                ppi(:,:,mi)=(f*R*R*dth*dphi/c0)*sum(pc.*exp(1j*(-kk*r + (pi/2)))./r);	
                waitbar(mi/mmax)												
            end
            close(hs);
            perfa=shiftdim(sum(ppi));			% perfa is pressure on ERFA plane, an lmax x mmax matrix.
            ElemLoc=[0 0];      % a solid transducer, so only one 'element'.
            button2='solid';

            hs1=imagesc(phi1,th1,pointmap); axis image; axis xy;
            title('Angle plot of circular transducer - face view'); xlabel('phi (radians)'); ylabel('theta (radians)');    
            figure; hs2=imagesc(xp,squeeze(yyp),abs(perfa)); axis image; axis xy;
            title('Space plot of ERFA pressure from transducer - face view'); xlabel('horizontal (m)'); ylabel('vertical (m)');

            perfa=fliplr(perfa);    % flip perfa so the x-axis (col) matches the direction of the x-axis in HAS.
          % ----------------  
          case 'Phased Array'; % The next steps model an array of elements, each either round with radius=relem OR
                               % square with width=height=2*relem:

            button2=questdlg('Does the phased array have round or square elements?','ELEMENT SHAPE','Round',...
                'Square','Cancel','Cancel');
            if strcmp(button2, 'Cancel'); return; end

            for ns = 1:length(Seg)
                % Rotate array segment into correct orientation, then save spherical coords
                % into thvect and phivect.

                Ry =        [cos(PHI_rotation(ns)) 0 sin(PHI_rotation(ns)); 
                            0 1 0; 
                            -sin(PHI_rotation(ns)) 0 cos(PHI_rotation(ns));];

            % Create the rotation matrix about X' axis by psi
                 Rx =        [1 0 0;
                             0 cos(TH_rotation(ns)) sin(TH_rotation(ns));
                              0 -sin(TH_rotation(ns)) cos(TH_rotation(ns));];
               
                a = Seg(ns).ElemLoc_Cart(:,1);
                b = Seg(ns).ElemLoc_Cart(:,2);
                c = Seg(ns).ElemLoc_Cart(:,3);
                rotated = [Rx*Ry*[a b c]']';
                [thvect, phivect] = cart2sphd(rotated(:,1),rotated(:,2),rotated(:,3), R);
                h = figure; h.Position = [1252          87         667         890];
                %plotTx(a,b,c,h)
                %hold on
                %plotTx(rotated(:,1),rotated(:,2),rotated(:,3),h);
                %view([1,0,0]);
                
                close(h);
                % a = Seg(ns).ElemLoc_Cart(:,1);
                % b = Seg(ns).ElemLoc_Cart(:,2);
                % c = Seg(ns).ElemLoc_Cart(:,3);
                % rotated = [Ry*Rz*[a b c]']';
                % [phivect, thvect,~] = cart2sph(rotated(:,1),rotated(:,2),rotated(:,3));
                
                    close all
                    plot3(fullarray(:,1),fullarray(:,2),fullarray(:,3),'g*');
                    axis image, grid on, hold on, xlabel('x'),ylabel('y'),zlabel('z')
                    plot3(a,b,c,'ko');
                    plot3(rotated(:,1),rotated(:,2),rotated(:,3),'r*');
                h1=waitbar(0,'Evaluating segment ERFA pressure pattern...');
                numelem=size(thvect,1);
                perfa=zeros(lmax,mmax,numelem,'single');     % preallocation for speed.         
                pointmap=zeros(imax,kmax,'single'); % initialize map of valid points inside xducer elements.


                % --- Loop through the elements (of the current segment) starting here ---
                for g=1:numelem       % cycle through elements, finding a page of ERFA for each element.
                    
                    distfromelem=R*acos((cos(thvect(g))*cos(thmesh1)).*cos(phivect(g)-phimesh1)...
                        +sin(thvect(g))*sin(thmesh1));   % great circle distance from center of element g.

                    indelemlin= find(distfromelem<=relem);  % find linear indices inside round element g.
                    [indelemi, indelemk]=ind2sub([imax,kmax],indelemlin);   % convert to subscripts.
                    is=min(indelemi);  ie=max(indelemi);    % find start and end indices in area around element g.
                    ks=min(indelemk);  ke=max(indelemk);
                    th2=th1(is:ie);     % angle row vector only around element.
                    phi2=phi1(ks:ke);     % angle row vector only around element.
                    isize=size(th2,2);      % size of rectangle enclosing element g only, so only integrate over element.
                    ksize=size(phi2,2);
    %                 if isize==1 || isize==0 || ksize==1 || ksize==0; 
    %                     errordlg(['Only 0 or 1 sample points inside at least one element. Increase number'...
    %                         ' of angle increments to transducer.']); return; end
                    thmesh2=repmat(th2',1,ksize); % isize x ksize meshgrid of theta values, only around element.
                    phimesh2=repmat(phi2,isize,1);  % isize x ksize meshgrid of phi values, only around element.
                    pt=zeros(imax,kmax,'single');        % start with blank pressure over entire field.  
                    if strcmp(button2,'Round')  % elements are round.
                        %TotElemArea=(pi*relem^2)*numelem;   % total area of all elements combined (legacy).                      
                        pt(indelemlin)=1;       % pt = 1 only where there are valid points inside element.                      
                    else                        % elements are square.
                        %TotElemArea=(4*relem^2)*numelem;   % total area of all elements combined (legacy).
                        pt(is:ie,ks:ke)=1;      % pt = 1 only where there are valid points inside element.                      
                    end

                    pt2=pt(is:ie,ks:ke);    % rectangular pressure matrix just around the element.
                    cth=cos(thmesh2);
                    ss5=R*(cth.*sin(phimesh2));
                    s=repmat(ss5,[1,1,lmax]);		% 3D array: isize x ksize x lmax.
                    xp=dxp*(m-round(mmax/2));		% vector of x points along (with vers. 8) the horizontal axis.

                    aa5=R*(cth.*cos(phimesh2));
                    a=repmat(aa5,[1,1,lmax]);		% 3D array: isize x ksize x lmax.
                    b=R-a;  clear a     % to save memory

                    tt= R*sin(thmesh2);
                    t=repmat(tt,[1,1,lmax]);		% 3D array of t.
                    yyp(1,1,:)=dyp*(l-round(lmax/2));	% turn y points into a 'page' vector, lmax pages long.
                    yp=repmat(yyp,[isize,ksize,1]);		% make 3D array.
                    term1=(t-yp).^2; clear t yp      % to save memory
                    term3=(d-b).^2;  clear b         % to save memory

                    ppc=pt2.*cth;					% ppc is isize x ksize matrix; cth is for spherical integration.
                    pc=repmat(ppc,[1,1,lmax]);		% product of pressure and cos theta now 3D array.
                    kk=2*pi*f/c0;
                    ppi=zeros(ksize,lmax,mmax,'single');     % prealocation for speed.

                    for mi=1:mmax
                         %--- Rayleigh-Summerfeld integral done next with three implicit for loops ---
                        r=sqrt(term1 + term3 + (s-xp(mi)).^2);	% r is 3D array for each value of xp.
                        ppi(:,:,mi)=(f*R*R*dth*dphi/c0)*sum(pc.*exp(1j*(-kk*r + (pi/2)))./r);	
                    end
                    perfa(:,:,g)=shiftdim(sum(ppi));

                    waitbar(g/numelem)
    %                 if g==1
    %                     hh1=imagesc(phi1,th1,pt); axis image; axis xy;
    %                     title('Angle plot of element 1 - face view'); xlabel('phi (radians)'); ylabel('theta (radians)');
    %                     figure; hh2=imagesc(xp,squeeze(yyp),squeeze(abs(perfa(:,:,g)))); axis image; axis xy;
    %                     title('Space plot of ERFA pressure from element 1 - face view'); 
    %                     xlabel('horizontal (m)'); ylabel('vertical (m)');
    %                     button6=questdlg('Do you want to continue with Rayleigh-Sommerfeld calculations?',...
    %                         'CONTINUE','Yes','Quit','Yes');
    %                     if strcmp(button6,'Quit'); close all; close(h1); return; end
    %                 end
                    pointmap=pointmap+pt;

                end         % end of g loop.

                TotalElemArea=R*R*dth*dphi*sum(sum(pointmap.*cos(thmesh1)));     % Add up areas of all points.
                ptelem=sqrt(2*Zm*1/TotalElemArea);      % pressure corresponding to 1 W total.
                perfa=perfa*ptelem;     % multiply normalized perfa by the pressure corresponding to 1 W total.

                close all; close(h1)
                hh3=imagesc(phi1,th1,pointmap); axis image; axis xy
                title('Angle plot of all elements - face view'); xlabel('phi (radians)'); ylabel('theta (radians)');

                
                phiInElem=phimesh1(pointmap>0);    % use logical indexing to find angles to points inside element.
                thetaInElem=thmesh1(pointmap>0);
                Xplot=R*cos(thetaInElem).*sin(phiInElem);
                Yplot=R*sin(thetaInElem);
                Zplot=R - R*cos(thetaInElem).*cos(phiInElem);  % put Z in ordinary direction coming out of xducer.
                hh4=plot3(Xplot,Yplot,Zplot,'.'); axis image;  axis xy % make 3D plot of points found in elements.
                title('3D space plot of all elements - face view'); xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); grid on;

                figure; hh5=imagesc(xp,squeeze(yyp),abs(sum(perfa,3))); axis image; axis xy
                title('Space plot of summed ERFA pressure from transducer - face view'); 
                xlabel('horizontal (m)'); ylabel('vertical (m)');

                perfa=fliplr(perfa);    % flip perfa so the x-axis (col) matches the direction of the x-axis in HAS.
                perfa_seg{ns} = perfa;
            end % end of s (segment) loop.
        end % end of switch

   
        button5 = questdlg('Do you wish to save the results?','Save?');
        if strcmp(button5,'Yes')
    %         [newfile,newpath] = uiputfile('ERFA8_.mat','Save binary ERFA_.mat file:');
    %         if newfile~=0;
                sx=d;
               % Save key parameters in .mat filese
                disp('Saving ERFA8.mat file...');
                save([fullfilepath filesep 'ERFA8.mat'],'perfa_seg','fMHz','Len','sx','R','isPA','fullarray','dxp','dyp','relem','Seg','TH_rotation','PHI_rotation'); 
                disp('Finished saving');
    %         end
        end
end   