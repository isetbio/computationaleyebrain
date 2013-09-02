function [theParams,staticParams] = setParameters(parameterPreset)
%  [theParams,staticParams] = setParameters([parameterPreset])
%
% Manage simulation parameters, return as two structs.  One (theParams)
% is variables that do or might someday vary within a simulation,
% the other (staticParams) is variables that remain constant throughout a single
% simulation.
%
% The division between theParams and staticParams is not completely rational,
% and carries some inertia from before this split was set up.  The split is
% important because some variables that get passed around are really big and
% we don't want multiple copies.
%
% Sets of parameters may be defined and invoked by the
% passed name.  Preset options are:
%   'BasicNoSurround'
%   'BasicRDrawSurround'
%   'BasicDetermSurround'
%   'BasicDetermSurroundWithNoise'
%   'MacularPigmentVary'
%   'QuickTest' [Default]
%
% See also constructSimulationParameters.  Changes here may require modifications there.
%
% 9/24/13  dhb  Pulled this out from main routine.
                            
%% Parameter section

% Visual system related
staticParams.integrationTimeSecs = 0.050;                    % Temporal integration time for detecting mechanisms.
staticParams.fieldOfViewDegrees = 2;                         % Field of view specified for the scenes.
staticParams.scenePixels = 64;                               % Size of scenes in pixels
staticParams.pupilDiameterMm = 3;                            % Pupil diameter.  Used explicitly in the PSF calc.
                                                             % [** Need to check that this is carried through to 
                                                             % the absorption calculations.  We might be using an isetbio
                                                             % default rather than the value set here.]
staticParams.coneProportions = [0.1 .6 .2 .1];               % Proportions of cone types in the mosaic, order: empty, L,M,S
staticParams.isetSensorConeSlots = [2 3 4];                  % Indices for LMS cones in iset sensor returns.   These run 2-4 because
                                                             % of the empty pixels
staticParams.coneApertureMeters = [sqrt(4.1) sqrt(4.1)]*1e-6;% Size of (rectangular) cone apertures, in meters.
                                                             % The choice of 4.1 matches the area of a 2.3 micron diameter IS diameter,
                                                             % and that is PTB's default.
                                                
% Stimulus related
staticParams.stimulus.type = 'rgb_uniform';                  % Set to allow different types of stimulus specification.
staticParams.stimulus.monitorName = 'LCD-Apple';             % Monitor spectrum comes from this file
staticParams.stimulus.backRGBValue = 0.5;                    % Define background for experment in monitor RGB
staticParams.stimulus.isetGammaValue = 2.2;                  % Needed to deal with gamma correction done on image read by iset.

staticParams.nColorDirections = 16;                          % Number of color directions for contour.
staticParams.dirAngleMax = 2*pi;                             % Use pi for sampling directions from hemicircle, 2*pi for whole circle
staticParams.nTestLevels = 10;                               % Number of test levels to simulate for each test direction psychometric function.
staticParams.nDrawsPerTestStimulus = 400;                    % Number of noise draws used in the simulations, per test stimulus
staticParams.criterionCorrect = 0.82;                        % Fraction correct for definition of threshold in TAFC simulations.

% Data management parameters
staticParams.theContourPlotLim = 0.5;                        % Axis limit for contour plots.
staticParams.outputRoot = 'out';                             % Plots get dumped a directory with this root name, but with additional
                                                             % characters to identify parameters of the run tacked on below.
staticParams.parameterPreset = parameterPreset;              % Name of preset used to determine entries of theParams.

% Convenience parameters
staticParams.nSensorClasses = length(staticParams.isetSensorConeSlots);
  
% Preset parameters, vary with preset name.
if (nargin < 1 || isempty(parameterPreset))
    parameterPreset = 'QuickTest';
end
switch (parameterPreset)
                            
    case 'NoSurroundNoSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1];         % Numbers of each cone class to use in the classifier.

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; % Simulate various tri and dichromats
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             % Can be true, false, or [true false]
        theParams.macularPigmentDensityAdjustments = [0];         % Amount to adjust macular pigment density for cone fundamentals of simulated observer.
                                                                  % Note that stimuli are computed for a nominal (no adjustment) observer.
                                                        
        theParams.noiseType = 1;                                  % Type of photoreceptor noise.  1 -> Poisson.  0 -> none.
        theParams.surroundType = 'none';                          % Define type of surround calc to implement
        theParams.surroundSize = 0;                               % Parameter defining surround size.
        theParams.surroundWeight = 0;                             % Parameter defining surround weight.
        theParams.integrationArea = 0;                            % Stimulus integration area.  NOT YET IMPLEMENTED.
        theParams.secondSiteFanoFactor = 0;                       % Noise added after opponent recombination, if any added.
                                                                  % Expressed as a fraction of Poisson variance to use.
 
    case 'NoSurroundWithSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [0];         
        
        theParams.noiseType = 1;  
        theParams.surroundType = 'none';                          
        theParams.surroundSize = 0;                               
        theParams.surroundWeight = 0;                             
        theParams.integrationArea = 0;                            
        theParams.secondSiteFanoFactor = 1;                                                                        
                                                       
    case 'RandomSurroundNoSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [0];
        
        theParams.noiseType = 1;
        theParams.surroundType = 'random_wiring';                         
        theParams.surroundSize = 10;                             
        theParams.surroundWeight = 0.7;                        
        theParams.integrationArea = 0;                            
        theParams.secondSiteFanoFactor = 0;
        
    case 'RandomSurroundWithSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [0]; 
        
        theParams.noiseType = 1;
        theParams.surroundType = 'random_wiring';                         
        theParams.surroundSize = 10;                             
        theParams.surroundWeight = 0.7;                        
        theParams.integrationArea = 0;                            
        theParams.secondSiteFanoFactor = 1;
        
     case 'SelectiveSurroundNoSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [0];
        
        theParams.noiseType = 1;
        theParams.surroundType = 'cone_specific';                         
        theParams.surroundSize = 10;                             
        theParams.surroundWeight = 0.7;                        
        theParams.integrationArea = 0;                            
        theParams.secondSiteFanoFactor = 0;
        
    case 'SelectiveSurroundWithSecondSiteNoise'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [0]; 
        
        theParams.noiseType = 1;
        theParams.surroundType = 'cone_specific';                         
        theParams.surroundSize = 10;                             
        theParams.surroundWeight = 0.7;                        
        theParams.integrationArea = 0;                            
        theParams.secondSiteFanoFactor = 1;

    case 'MacularPigmentVary'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        theParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];             
        theParams.macularPigmentDensityAdjustments = [-0.3 0 0.3];
        
        theParams.noiseType = 1;
        theParams.surroundType = 'none';                          
        theParams.surroundSize = 0;                              
        theParams.surroundWeight = 0;                             
        theParams.integrationArea = 0;                           
        theParams.secondSiteFanoFactor = 0;
    
    case 'QuickTest'
        staticParams.stimulus.coneNumbersToUse = [4 2 1]; 

        staticParams.nColorDirections = 4;
        staticParams.dirAngleMax = pi;
        staticParams.nTestLevels = 4;
        staticParams.nDrawsPerTestStimulus = 100;
        
        theParams.OBSERVER_STATES = {'LMandS'};
        theParams.DO_TAFC_CLASSIFIER_STATES = [true];
        theParams.macularPigmentDensityAdjustments = [0];
        
        theParams.noiseType = 1;
        theParams.surroundType = 'none';                        
        theParams.surroundSize = 0;                              
        theParams.surroundWeight = 0;                            
        theParams.integrationArea = 0;                         
        theParams.secondSiteFanoFactor = 1;                      
                                                       
    otherwise
        error('Unknown parameter presets');
end

