    pOutCorrected = pout; 

    usePhaseCorrection = false;  
    pR.usePhaseCorrection = false; 

    [pout2,Z,angpgvect,phasematsv,pR,pHAS] = pHAS_NoGUIFullfunc(Modl,pHAS,pR,pERFA);

    usePhaseCorrection = true;  
    pR.usePhaseCorrection = true; 

   