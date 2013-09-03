function simParams = constructSimulationParameters(theParams,staticParams)
% simParams = constructSimulationParameters(theParams,staticParams)
%
% Deal out the parameters structures to create a struct
% array of parameters, one entry of which will be used
% in each (potentially parallel) simulation run.
%
% See also setParameters.  Changes here may require modifications there.
%
% 8/24/13  dhb  Pulled this out.

%% Construct list of conditions
%
% These are strung out so that we can chunk through
% them in a big parfor loop below.
cdAngles = linspace(0,staticParams.dirAngleMax,staticParams.nColorDirections+1);
cdAngles = cdAngles(1:end-1);
testLevels = linspace(0,1,staticParams.nTestLevels);
paramIndex = 1;
for os = 1:length(theParams.OBSERVER_STATES)
    for ct = 1:length(theParams.DO_TAFC_CLASSIFIER_STATES)
        for m = 1:length(theParams.macularPigmentDensityAdjustments)
            for cdi = 1:staticParams.nColorDirections
                simParams(paramIndex).OBSERVER_STATE = theParams.OBSERVER_STATES{os};
                simParams(paramIndex).DO_TAFC_CLASSIFIER = theParams.DO_TAFC_CLASSIFIER_STATES(ct);
                simParams(paramIndex).macularPigmentDensityAdjust = theParams.macularPigmentDensityAdjustments(m);
                simParams(paramIndex).cdAngle = cdAngles(cdi);
                
                % Nuisance parameters
                
                % Set fixed params.  These can be made variable by adding a loop in this routine
                simParams(paramIndex).nTestLevels = staticParams.nTestLevels;
                simParams(paramIndex).noiseType = theParams.noiseType;
                simParams(paramIndex).surroundType = theParams.surroundType;
                simParams(paramIndex).surroundSize = theParams.surroundSize;
                simParams(paramIndex).surroundWeight = theParams.surroundWeight;
                simParams(paramIndex).integrationArea = theParams.integrationArea;
                simParams(paramIndex).secondSiteFanoFactor = theParams.secondSiteFanoFactor;
                
                % Kluge for now to select subregion of total cones
                simParams(paramIndex).fractionUse = 0.005;
                
                % Bump counter
                paramIndex = paramIndex+1;
            end
        end
    end
end