function [maxT,tempsAll, t,x,y] = getThermometryInfo(pathTreat,patnum)
cd(pathTreat)
addpath('/v/raid2/hodeen/ScanData/ReconCode/_MRTIcode')
addpath('/v/raid2/hodeen/ScanData/ReconCode/CntrGUI/')
addpath('/v/raid2/hodeen/ScanData/ReconCode/')
clear tempsAll
clear temps
clear maxT
clear x
clear y
clear z
clear t
% Grab the folder names
SonicationFolders=grabFiles(pwd,'Sonic','d','d');
% Remove non-folders (here just checking if there's a '.' as in an extension...
SonicationFolders = SonicationFolders(~contains(SonicationFolders, '.'));

% Sort the folder names
for ii=1:length(SonicationFolders)
    SonicationFolders2(ii,1)=str2num(SonicationFolders{ii}(11:end));
end
SonicationFolders2=sort(SonicationFolders2);
for ii=1:length(SonicationFolders2)
    SonicationFolders{ii} = ['Sonication',num2str(SonicationFolders2(ii))];
end; clear ii SonicationFolders2

if(patnum == 23 || patnum == 217 || patnum == 152 || patnum == 162 || patnum ==196 || patnum==201 ||patnum==204)
    numberFolders = length(SonicationFolders) -1;
else
    numberFolders = length(SonicationFolders);
end
for sonication=1:numberFolders
    cd(SonicationFolders{sonication}) % Go in to the folder
    pause(0.5)

    % Loop over and read in MRTI data
    % Grab all files containing MRTI data and open them
    files=grabFiles(pwd,'5-','raw',[1 2]);
    for ii=1:length(files)
        [tmp] = fopen(cell2mat(files(ii))); [temps(:,ii)] = fread(tmp,'float');
    end; clear ii tmp
    % Reshape
    temps = reshape(temps, [256 256 size(temps,2)]);

    % Save to structure
    tempsAll{sonication}=temps;

    clear temps files

    cd ..
end

for ii=1:numel(tempsAll)
    [x{ii} y{ii} z{ii} t{ii} maxT{ii}] = findHottestVoxel(permute(tempsAll{ii},[1 2 4 3]), [128-3 128+3 128-3 128+3 1 1 1 size(tempsAll{ii},3) ]);
end
end