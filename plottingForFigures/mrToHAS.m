function [mrFocsInHAS] = mrToHAS(KRXfile, KRXpath, attributeData)
        KRXsubfolder = KRXfile(1:end-4);
        %mrFocsInHAS = {};
        
        disp('Parsing KRX data.')
        
        %unzip([KRXpath KRXfile], KRXsubfolder);
        cd([KRXpath KRXsubfolder]);
        

        %% Open XML file and extract CT features
        FID = fopen('model.xml');
        TEXT = fscanf(FID,'%s');
        k1 = strfind(TEXT,'<CT__Image>');
        k2 = strfind(TEXT,'</CT__Image>');
        s1 = strfind(TEXT,'<Sonications');
        s2 = strfind(TEXT,'</Sonications>');
        loc1 = s2+1;
        loc2 = strfind(TEXT,'</Transducer>');
        CTxml = TEXT(k1:k2); clear k1 k2;
        SONxml = TEXT(s1:s2); clear s1 s2;
        Txdrxml = TEXT(loc1:loc2); clear loc1 loc2; fclose(FID); clear FID

        % Data within Dimensions:
        % Extract CT size
        CTsize = regexp(CTxml,'-*\d*\.*\d*(?=</dimensionSize>)'); 
        for i = 1:length(CTsize)
            Plus = strfind(CTxml(CTsize(i):CTsize(i)+15),'</dim')-1;
            Plus = CTsize(i)+Plus(1)-1;
            Size(i)=cell2mat(textscan(CTxml(CTsize(i):Plus),'%d'));
        end; clear CTsize; 

        % Extract CT spacing
        spacing = regexp(CTxml,'-*\d*\.*\d*(?=</regularSpacing>)');
        for i = 1:length(spacing)
            Plus = strfind(CTxml(spacing(i):spacing(i)+15),'</reg')-1;
            Plus = spacing(i)+Plus(1)-1;
            Spacing(i)=cell2mat(textscan(CTxml(spacing(i):Plus),'%f'));
        end; clear spacing

        % Extract CT position
        startF = '-*\d*\.*\d*(?=</float>)'; % regex: all signed floats followed by '</float>'
        endF = '</float>';                  % regex: ending '</float>'

        pos     = strfind(CTxml,'"ImagePosition'); 
        posPlus = strfind(CTxml(pos:end),'</Attribute>'); posPlus = pos+posPlus(1);
        pos1 = pos+regexp(CTxml(pos:posPlus),startF)-1;
        pos2 = pos+regexp(CTxml(pos:posPlus),endF)-1;
        for i=1:length(pos1)
            Pos(i)=cell2mat(textscan(CTxml(pos1(i):pos2(i)),'%f'));
        end; clear pos*

        % Extract CT orientation
        orient     = strfind(CTxml,'"ImageOrientation"'); 
        oriPlus = strfind(CTxml(orient:end),'</Attribute>'); oriPlus = orient+oriPlus(1);
        ori1 = orient+regexp(CTxml(orient:oriPlus),startF)-1;
        ori2 = orient+regexp(CTxml(orient:oriPlus),endF)-1;
        for i=1:length(ori1)
            Orient(i)=cell2mat(textscan(CTxml(ori1(i):ori2(i)),'%f'));
        end; clear ori* startF endF Plus i ans FID

        % Extract MR transform matrix "mm"
        startF = '<m\d{2}>';
        endF = '<\/m\d{2}>';

        MRrot   = strfind(CTxml,'ImageTransformMatExt');
        rotStarts = MRrot+regexp(CTxml(MRrot:end),startF)+4;
        rotEnds = MRrot+regexp(CTxml(MRrot:end),endF)-2;
        try
            for i=1:4
                for j=1:4
                    n = (i-1)*4+j;
                    mm(j,i) = cell2mat(textscan(CTxml(rotStarts(n):rotEnds(n)),'%f'));
                end
            end; clear rot* n j
        catch 
            mm = [ ];
        end 

        % Extract FUS element locations for HAS
        % Trans = strfind(,'<Transducername="');
        disp('Collecting FUS Per-Element Info')
        startF = '(?<=([xyz]=")).\d*.\d*';% Matches the numerical contents following any [xyz]= until the ending quotation mark. This includes nx=", so we have to search for those specifically to exclude them...
        tranStarts = regexp(Txdrxml,startF);
        nexp = '(?<=(n[xyz]=")).\d*.\d*'; % Matches the nx=", ny=", nz=" to exclude from regexp results.
        tranStartsNxyz = regexp(Txdrxml,nexp);
        tranStarts = setxor(tranStarts, tranStartsNxyz); % Find entries that only occur in one set (so exclude the entries that occur in both sets)      
        loc = zeros(1,3);
        for i = 0:length(tranStarts)/3-1
            start = tranStarts(i*3+1);            
            stop = strfind(Txdrxml(tranStarts(i*3+1):tranStarts(i*3+1)+12),'"'); stop = stop(1)-2;
            x = Txdrxml(start:start+stop);    
            loc(i+1,1) = str2num(x);

            start = tranStarts(i*3+2);
            stop = strfind(Txdrxml(tranStarts(i*3+2):tranStarts(i*3+2)+12),'"'); stop = stop(1)-2;
            y = Txdrxml(start:start+stop);    
            loc(i+1,2) = str2num(y);

            start = tranStarts(i*3+3);
            stop = strfind(Txdrxml(tranStarts(i*3+3):tranStarts(i*3+3)+12),'"'); stop = stop(1)-2;
            z = Txdrxml(start:start+stop);    
            loc(i+1,3) = str2num(z);     

            loc(i+1, 4) = i+1;
        end
        loc = loc/1000; % convert to meters   
        
        %% Load CT volume and compute focal spot (no longer: transform to correct coordinate space)
        disp('Loading CT Volume from Kranion Export')
        FID = fopen([KRXpath KRXsubfolder '/imageData/imageChannel0']);
        img = fread(FID,'ushort'); fclose(FID); clear FID;
        V = reshape(img,[Size(1) Size(2) Size(3)]);       
        B = permute(V,[2 1 3]); 
        ScaleFactor = 1; %(1 gives 1mm voxels)
        Scale = [...
                Spacing(1)/Spacing(3)/ScaleFactor 0 0 0;...
                0 Spacing(2)/Spacing(3)/ScaleFactor 0 0;...
                0 0 Spacing(3)/Spacing(3)/ScaleFactor 0;...
                0 0 0 1];       
        netRot = affine3d(Scale);
        Vhas = imwarp(B, netRot); 
        %[xdim,ydim,zdim] = size(B);
        %Vhas = TriLinInterpV3(B,Scale^(-1),zeros(ceil(xdim*Spacing(1)),ceil(ydim*Spacing(2)),ceil(zdim*Spacing(3))));
        unTransform = inv(Scale);
        


        %% Parse XML for sonication details (XTilt, Power, NaturalFocus,
        % FocusSteering, per-element Phase, per-element Amplitude
        SonNum = regexp(SONxml,'(?<=selectedSonication=").{1,2}','match');
        SonNum = cell2mat(textscan(SonNum{1},'%f'));

        tokens = regexp(SONxml, 'Sonicationscount="(\d+)"', 'tokens');
        NumSon = str2double(tokens{1}{1});
        
        for j = 0:(NumSon-1)
            %if isfile([[KRXpath KRXfile(1:end-4)] filesep ['SonicationParameters',num2str(j+1),'.mat']])
               %continue;
            %end
            SonNum = j;
            Sonication = struct([]);
            power = '(?<=power__w=")-?\d*.\d*';
            duration = '(?<=duration__s=")-?\d*.\d*';
            natural = '(?<=<.>)-*\d*.\d*(?=<\/.>)';
            xtilt = '(?<=transducerXTilt.*>)-*\d*.\d*(?=<\/Attribute>)';
            phase = '(?<=phase=")-*\d*\.?\d*(?:[eE][+-]?\d+)?(?=")';
            amplitude = '(?<=amplitude=")-*\d*\.?\d*(?:[eE][+-]?\d+)?(?=")';
    
            start = strfind(SONxml,['Sonicationnum="' num2str(SonNum) '"']);
            stop = strfind(SONxml(start:end),'</Sonication>')+start;
        
            Sonication(1).xml = SONxml(start:stop(1));
            Sonication.Power = str2double(cell2mat(regexp(Sonication.xml,power,'match')));
            Sonication.Duration = str2double(cell2mat(regexp(Sonication.xml,duration,'match')));
            Sonication.NatFoc = regexp(Sonication.xml,natural,'match'); % First three entries are x,y,and z, next three are x,y and z of steering... remainder are nonsense.
            Sonication.SteeredFoc(1) = str2double(Sonication.NatFoc{4});
            Sonication.SteeredFoc(2) = str2double(Sonication.NatFoc{5});
            Sonication.SteeredFoc(3) = str2double(Sonication.NatFoc{6});
    
    %        %%%%%%% NOT NEEDED (because I have per-element phase) %%%%%%%% Sonication.FocusSteering = ; %FocusSteering
            Sonication.XTilt = regexp(Sonication.xml,xtilt,'match'); % First entry is the tilt for this sonication...remainder are nonsense.
            try
                Sonication.XTilt = str2double(Sonication.XTilt{1});
            catch
                Sonication.XTilt = (input('Enter xtilt>>'));
            end
            Sonication.Phase = regexp(Sonication.xml,phase,'match');
            Sonication.Amplitude = regexp(Sonication.xml,amplitude,'match');
            
            Sonication.Phase = str2double(Sonication.Phase)';
            Sonication.Amplitude = str2double(Sonication.Amplitude)';
            try
                Sonication.Phase(:,2) = 1:1024;
                Sonication.Amplitude(:,2) = 1:1024;
            catch 
                continue;
            end 
                % end
            
            % Set dimension vectors x,y,z based on matrix size and CT spacing
            x = 0:Size(1)-1; 
            x = x*Spacing(1);
            y = x;
            z = 0:Size(3)-1; 
            z = z*Spacing(3);
            
            [X,Y,Z] = meshgrid(x,y,z);
            
            fx = str2num(Sonication.NatFoc{1});
            fy = str2num(Sonication.NatFoc{2});
            fz = str2num(Sonication.NatFoc{3});
    
            % Create a vector defining the natural focus of the transducer in world
            % coordinates (in millimeters)
            Foc = [fx fy fz 1]';
            Sonication.NatFoc = Foc;

            % Variable mm is the transform from imported CT coordinates to
            % the "World".
            % Use the inverse of mm to convert Foc (in "world") to imported coords.
            try
                mrFoc = attributeData(j+1,:);
            catch
                continue;
            end
            if sign(mrFoc(1)) ~= sign(Foc(1))
                mrFoc(1) = -mrFoc(1);
            end
            if sign(mrFoc(2)) ~= sign(Foc(2))
                mrFoc(2) = -mrFoc(2);
            end
            if sign(mrFoc(3)) ~= sign(Foc(3))
                mrFoc(3) = -mrFoc(3);
            end
            
            FocCT = inv(mm) * [mrFoc 1]'; %% FocCT is in millimeters 
            
            xscal = x(1:length(Vhas(1,:,1)))/Scale(1,1);
            yscal = y(1:length(Vhas(:,1,1)))/Scale(2,2);
            zscal = z(1:length(Vhas(1,1,:)))/Scale(3,3);
           
            [~, FocCTijk(1)] = min(abs(xscal - round(FocCT(1))));
            [~, FocCTijk(2)] = min(abs(yscal - round(FocCT(2))));
            [~, FocCTijk(3)] = min(abs(zscal - round(FocCT(3))));


            mrFocsInHAS{j+1} = FocCTijk;
            
         
        end 
        sm=size(Vhas); % sm (in y,x,z order) is used extensively throughout the gui.
        if 2*round(sm(1)/2)==sm(1) 
            Vhas=padarray(Vhas,[1 0 ],1,'post');
        end
        if 2*round(sm(2)/2)==sm(2)
            Vhas=padarray(Vhas,[0 1 ],1,'post');
        end
end
