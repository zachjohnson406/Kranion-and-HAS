function  pERFA = FUSF_loadERFA(fullfilepath,segment)

%loadERFA  Summary of this function goes here
%load the ERFA file saved by ERFAMaker8.m

%   Detailed explanation goes here
%            save(fullfile,'perfa','fMHz','Len','sx','R','isPA','ElemLoc','pfilename','dxp','dyp','relem') 
load([fullfilepath filesep 'ERFA8.mat']);
pERFA.perfa = perfa_seg{segment};
pERFA.fMHz = fMHz;
pERFA.Len = Len;
pERFA.sx = sx;
pERFA.R = R;
pERFA.isPA = isPA;
pERFA.ElemLoc = Seg(segment).ElemLoc_Polar;
% pERFA.pfilename = pfilename;
pERFA.dxp = dxp;
pERFA.dyp = dyp;
pERFA.relem = relem;

end

