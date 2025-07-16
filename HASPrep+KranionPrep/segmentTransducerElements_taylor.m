% Segment Array
% 
% This function segments an Insigtec transducer 7 segments. Uses distance
% beteween elements to determine 7 element clusters. 
% 
% @INPUTS
%   loc: x,y,z location of elements (in m)
% 
% @OUTPUTS
%   Seg: Transducer elements divided into 7 sections
%   TH_centers: Center of each segment theta coord
%   PHI_centers: Center of each sement phi coord
%   TH_rotation: Rotation of elements about y axis before running HAS
%   PHI_rotation: Rotation of elements about x axis before running HAS
% 
%   Taylor Web and Zach Johnson 
%   June 19th, 2024 

function [Seg, TH_rotation, PHI_rotation] = segmentTransducerElements_taylor(loc)

        sorted = sort(loc(:,3));
        dsorted = sorted(2:end)-sorted(1:end-1);
        I = find(dsorted == max(dsorted(1:end-100)));
        zThresh = sorted(I,1); %#ok<FNDSB>

        [el_num,I] = find(loc(:,3)<=zThresh);

        % To rotate into HAS +x propagation orientation:
        % y = y; x = -z; and z = x
        [TH,PHI,~] = cart2sph(-loc(el_num,3),loc(el_num,2),loc(el_num,1));
        Seg(1).ElemLoc_Polar = [TH, PHI];
        Seg(1).ElemLoc_Cart = [loc(el_num,1),loc(el_num,2),loc(el_num,3),loc(el_num,4)];
        
        % Specify Seg2-SegN, the outer "ring" of segments
        % First remove the elements in Seg1
        I = find(loc(:,3)>zThresh);
        loc = loc(I,:); %discard current segment of loc
        [TH,PHI,R] = cart2sph(loc(:,1),loc(:,2),loc(:,3));
 
        % Specify Theta angles for dividing the array into segments
        % Note: this is for use when the array is in "salad bowl" mode,
        % proagation along +z axis. This facilitates division into sectors
        % around Theta (azimuth).
        TH_div = (pi/6:pi/3:2*pi)-pi;
        TH_centers = [0 -pi:pi/3:2*pi/3];
        PHI_centers = [-pi/2 asin(zThresh/max(R))/2*ones(1, length(TH_centers)-1)];

        % Now we set a point at the center of each of the segments, which will
        % allow us to determine the Theta Phi rotations after rotating the
        % transducer into HAS propagation orientation (along minus-x axis)
        [x,y,z] = sph2cart(TH_centers,PHI_centers,150*ones(1,length(TH_centers)));
        [TH_centers_rot, PHI_centers_rot] = cart2sphd(x,y,z,150); %rotate centerpoints and compute centers in spherical coords.
        TH_rotation = TH_centers_rot; %Theta = pi represents no rotation (i.e. center of the array) in HAS orientation, so tip-minus-tail.
        PHI_rotation = PHI_centers_rot; % Phi = pi represents no rotation (ie. center of the array) in HAS orientation, so tip-minus-tail
        Segment_Centers = [TH_rotation', PHI_rotation'];
        for i = 1:7
            Seg(i).Center_Polar = Segment_Centers(i,:);
            Seg(i).Center_Cart = [-z(i),y(i),x(i)];
        end

        
        % Color key for segment plotting
        c = 'bgrcmkb';
        figure, axis image, grid on, hold on, xlabel('x'),ylabel('y'),zlabel('z') 
        xlim([-.150 .150]),ylim([-.150 .150]),zlim([-.150 .150])
        title('Loc before storing in Seg struct')

        % Loop through outer ring to set those array segments for ERFAmaker
        for i = 1:6
            I = find(TH < TH_div(i)); %Indices for current segment of elements.
            if i == 1
                I2 = find(TH > TH_div(end));
                I = [I; I2];
            end
            [T, P,~] = cart2sph(-loc(I,3),loc(I,2),loc(I,1));   % X and Z swapped to put transducer in HAS minus-x propagtion orientation
            Seg(i+1).ElemLoc_Polar = [T P];
            Seg(i+1).ElemLoc_Cart = [loc(I,1),loc(I,2),loc(I,3),loc(I,4)];             
            plot3(loc(I,1),loc(I,2),loc(I,3),[c(i) 'o']);

            % verify segment
            I = find(TH >= TH_div(i));
            TH = TH(I); PHI = PHI(I);%discard current segment of TH and PHI.
            loc = loc(I,:);   
        end

        % Check segments if needed - as of here, array is located on +x axis,
        % aimed in -x direction (at origin)
        figure, axis image, grid on, hold on, xlabel('x'),ylabel('y'),zlabel('z') 
        xlim([-.150 .150]),ylim([-.150 .150]),zlim([-.150 .150])
        title('Elem Loc loaded from Seg structure')
        for i = 1:7  
             plot3(Seg(i).ElemLoc_Cart(:,1), Seg(i).ElemLoc_Cart(:,2), Seg(i).ElemLoc_Cart(:,3),[c(i) 'o']) 
             % Seg(i).ElemLoc_Polar = fliplr(Seg(i).ElemLoc_Polar);
             [tmp1,tmp2] = cart2sphd(Seg(i).ElemLoc_Cart(:,1),Seg(i).ElemLoc_Cart(:,2),Seg(i).ElemLoc_Cart(:,3),0.15);
             Seg(i).ElemLoc_Polar = [tmp1,tmp2];
        end
              