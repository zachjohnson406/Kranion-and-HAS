% Based on CTtoHASmodl - from dcmCTtoXDbyA - reads in CT and MR volume images and uses the Insightec
% determined afine transformation to map CT to MR and to map CT and MR to
% the Insightec transducer coordinates
%
% @INPUTS
%   Vhas: CT from Kranion.
%   fullfilepath: Path to Kranion export file in use
% 
% @OUTPUTS
%   Modl: Segmented HAS model. 
% 
%   Zach Johnson and Matt Eames
%   June 19th, 2024 

function Modl = FUSF_Make_Modl(Vhas,fullfilepath)
    %% Segment CT volume into bone, diploe, brain, fat, and water volumes
   
    % Assign tissue types to CT value ranges
    thrsh_bone = 2000; % Utah: 1000;
    thrsh_diploe = 1550; % Utah: 200;
    thrsh_brain = 1000; % Utah: -20;
    thrsh_fat = 800; % Utah: -200;
   
    % Threshold volume
    mskbone = Vhas >= thrsh_bone;
    mskdiploe = (thrsh_bone > Vhas) & (Vhas >= thrsh_diploe);
    mskbrain = (thrsh_diploe > Vhas) & (Vhas >= thrsh_brain);
    mskfat = (thrsh_brain > Vhas) & (Vhas >= thrsh_fat);
    mskwater = Vhas < thrsh_fat;       %for Insightec, air is replaced by water
    
    % Assign segment flag/ID
    ctbone = mskbone .* 6;
    ctdiploe = mskdiploe .* 7;
    ctbrain = mskbrain .*5;
    ctfat = mskfat .* 3;
    ctwater = mskwater .* 1;
    
    % Combine segments into unified Modl variable
    Modl = ctwater + ctfat + ctbrain + ctdiploe + ctbone;

    %% Refine the Modl - Remove partial volume voxels.
    new_Modl = Modl;


    %% Update Modl and save
    Modl = new_Modl;


end

