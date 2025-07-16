% Matches element location to corresponding phase and amplitude
% or performs steering if given a target location. 
%
% @OUTPUTS
%   pR: A struct containing infomation relating to element location, phase
%   and amplitude. 
%
% @INPUTS
%   segmentedElementLocations: A stucture with each of the 7 segments of
%   the transducer. 
%   phase: NumberOfElements x 2 matrix with element phase in first column 
%   and element number in the second. 
%   amplitude: NumberOfElements x 2 matrix with element amplitude in first column 
%   and element number in the second. 
%   target: optional parameter, steers focus to target by altering phase, if entered. 
%
%
% Zach Johnson, Matt Eames, and Taylor Webb
% June 28th, 2024

function pR = matchPhaseToSections(segmentedElementLocations,phase,amplitude,target)
      
        %% Steering
        if exist('target','var')
            for ii = 1:length(segmentedElementLocations)
                phaseAlt = SteeringPhasesPA8(target(2),target(1),target(3),0.15,segmentedElementLocations(ii).ElemLoc_Polar,670000,1500);
                pR.Phases(ii).AmpPh1 = [ones(size(phaseAlt)),phaseAlt];

            end
        else
            phase = phase(:,1);
            amplitude = amplitude(:,1);
            for ii = 1:length(segmentedElementLocations)
                idx = round([segmentedElementLocations(ii).ElemLoc_Cart(:,4)]*1e3);
                phaseAlt = phase(idx);
                ampAlt = amplitude(idx)/norm(amplitude(idx));
                pR.Phases(ii).AmpPh1 = [ampAlt,phaseAlt];
            end
        end

        