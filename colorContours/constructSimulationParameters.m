function simParams = constructSimulationParameters(theParams,staticParams)
% simParams = constructSimulationParameters(theParams,staticParams)
%
% Deal out the parameters structures to create a struct
% array of parameters, one entry of which will be used
% in each (potentially parallel) simulation run.
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
                for t = 1:length(testLevels)
                    simParams(paramIndex).OBSERVER_STATE = theParams.OBSERVER_STATES{os};
                    simParams(paramIndex).DO_TAFC_CLASSIFIER = theParams.DO_TAFC_CLASSIFIER_STATES(ct);
                    simParams(paramIndex).macularPigmentDensityAdjust = theParams.macularPigmentDensityAdjustments(m);
                    simParams(paramIndex).cdAngle = cdAngles(cdi);
                    simParams(paramIndex).testLevel = testLevels(t);
                    
                    % Set test contrast maximum length.
                    %
                    % The best values depend on observer state
                    % and are currently set manually based on
                    % experience.
                    switch (simParams(paramIndex).OBSERVER_STATE)
                        case 'LMandS'
                            simParams(paramIndex).testContrastLengthMax = 0.3;
                        case 'LSonly'
                            simParams(paramIndex).testContrastLengthMax = 1;
                        case 'MSonly'
                            simParams(paramIndex).testContrastLengthMax = 1;
                        otherwise
                            error('Unknown dichromat/trichromat type specified');
                    end
                    
                    % Nuisance parameters
                    
                    % Set fixed params.  These can be made variable by adding a loop in this routine
                    simParams(paramIndex).noiseType = theParams.noiseType;
                    simParams(paramIndex).surroundType = theParams.surroundType;
                    simParams(paramIndex).surroundSize = theParams.surroundSize;
                    simParams(paramIndex).surroundWeight = theParams.surroundWeight;
                    simParams(paramIndex).integrationArea = theParams.integrationArea;
                    simParams(paramIndex).opponentLevelNoiseSd = theParams.opponentLevelNoiseSd;
                    
                    % Kluge for now to select subregion of total cones
                    simParams(paramIndex).fractionUse = 0.005;
                   
                    % Bump counter
                    paramIndex = paramIndex+1;
                end
                
            end
        end
    end
end