function [tParams,sParams] = setParameters(parameterPreset)
% [tParams, sParams] = setParameters([parameterPreset])
%
% Manage simulation parameters, return as two structs.  One (theParams) is
% variables that do or might someday vary within a simulation, the other
% (staticParams) is variables that remain constant throughout a single
% simulation.
%
% The division between theParams and staticParams is not completely
% rational, and carries some inertia from before this split was set up.
% The split is important because some variables that get passed around are
% really big and we don't want multiple copies. (But maybe Matlab doesn't
% create the extra copy unless you touch it ... we will see).
%
% Sets of parameters may be defined and invoked by the passed name.  Preset
% options are (you can have spaces and upper lower.  Everything squeezed by
% ieParamFormat)
%
%   'Basic No Surround'
%   'Basic RDraw Surround'
%   'Basic Determ Surround'
%   'Basic Determ Surround With Noise'
%   'Macular Pigment Vary'
%   'Quick Test' [Default]
%
% See also constructSimulationParameters.  Changes here may require
% modifications there.
%
% 9/24/13  dhb  Pulled this out from main routine.
                            
%% Parameter initialization

if ieNotDefined('parameterPreset'), parameterPreset = 'Quick Test'; end 

% Visual system related
sParams.integrationTimeSecs = 0.050;                    % Temporal integration time for detecting mechanisms.
sParams.fieldOfViewDegrees = 2;                         % Field of view specified for the scenes.
sParams.scenePixels = 64;                               % Size of scenes in pixels
sParams.pupilDiameterMm = 3;                            % Pupil diameter.  Used explicitly in the PSF calc.
                                                             % [** Need to check that this is carried through to 
                                                             % the absorption calculations.  We might be using an isetbio
                                                             % default rather than the value set here.]
sParams.coneProportions = [0.1 .6 .2 .1];               % Proportions of cone types in the mosaic, order: empty, L,M,S
sParams.isetSensorConeSlots = [2 3 4];                  % Indices for LMS cones in iset sensor returns.   These run 2-4 because
                                                             % of the empty pixels
sParams.coneApertureMeters = [sqrt(4.1) sqrt(4.1)]*1e-6;% Size of (rectangular) cone apertures, in meters.
                                                             % The choice of 4.1 matches the area of a 2.3 micron diameter IS diameter,
                                                             % and that is PTB's default.
                                                
% Stimulus related
sParams.stimulus.type = 'rgb_uniform';                  % Set to allow different types of stimulus specification.
sParams.stimulus.monitorName = 'LCD-Apple';             % Monitor spectrum comes from this file
sParams.stimulus.backRGBValue = 0.5;                    % Define background for experment in monitor RGB
sParams.stimulus.isetGammaValue = 2.2;                  % Needed to deal with gamma correction done on image read by iset.

sParams.nColorDirections = 16;                          % Number of color directions for contour.
sParams.dirAngleMax = 2*pi;                             % Use pi for sampling directions from hemicircle, 2*pi for whole circle
sParams.nTestLevels = 10;                               % Number of test levels to simulate for each test direction psychometric function.
sParams.nDrawsPerTestStimulus = 400;                    % Number of noise draws used in the simulations, per test stimulus
sParams.criterionCorrect = 0.82;                        % Fraction correct for definition of threshold in TAFC simulations.

% Data management parameters
sParams.outputRoot = 'out';                             % Plots get dumped a directory with this root name, but with additional
                                                             % characters to identify parameters of the run tacked on below.
sParams.parameterPreset = parameterPreset;              % Name of preset used to determine entries of theParams.

% Convenience parameters
sParams.nSensorClasses = length(sParams.isetSensorConeSlots);

%%
parameterPreset = ieParamFormat(parameterPreset);
switch (parameterPreset)
                            
    case 'nosurroundnosecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1];         % Numbers of each cone class to use in the classifier.

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; % Simulate various tri and dichromats
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             % Can be true, false, or [true false]
        tParams.macularPigmentDensityAdjustments = 0;         % Amount to adjust macular pigment density for cone fundamentals of simulated observer.
                                                                  % Note that stimuli are computed for a nominal (no adjustment) observer.
                                                        
        tParams.noiseType = 1;                                  % Type of photoreceptor noise.  1 -> Poisson.  0 -> none.
        tParams.surroundType = 'none';                          % Define type of surround calc to implement
        tParams.surroundSize = 0;                               % Parameter defining surround size.
        tParams.surroundWeight = 0;                             % Parameter defining surround weight.
        tParams.integrationArea = 0;                            % Stimulus integration area.  NOT YET IMPLEMENTED.
        tParams.secondSiteFanoFactor = 0;                       % Noise added after opponent recombination, if any added.
                                                                  % Expressed as a fraction of Poisson variance to use.
 
    case 'nosurroundwithsecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = 0;         
        
        tParams.noiseType = 1;  
        tParams.surroundType = 'none';                          
        tParams.surroundSize = 0;                               
        tParams.surroundWeight = 0;                             
        tParams.integrationArea = 0;                            
        tParams.secondSiteFanoFactor = 1;                                                                        
                                                       
    case 'randomsurroundnosecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = 0;
        
        tParams.noiseType = 1;
        tParams.surroundType = 'random_wiring';                         
        tParams.surroundSize = 10;                             
        tParams.surroundWeight = 0.7;                        
        tParams.integrationArea = 0;                            
        tParams.secondSiteFanoFactor = 0;
        
    case 'randomsurroundwithsecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = 0; 
        
        tParams.noiseType = 1;
        tParams.surroundType = 'random_wiring';                         
        tParams.surroundSize = 10;                             
        tParams.surroundWeight = 0.7;                        
        tParams.integrationArea = 0;                            
        tParams.secondSiteFanoFactor = 4;
        
     case 'selectivesurroundnosecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = 0;
        
        tParams.noiseType = 1;
        tParams.surroundType = 'cone_specific';                         
        tParams.surroundSize = 10;                             
        tParams.surroundWeight = 0.7;                        
        tParams.integrationArea = 0;                            
        tParams.secondSiteFanoFactor = 0;
        
    case 'selectivesurroundwithsecondsitenoise'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = 0; 
        
        tParams.noiseType = 1;
        tParams.surroundType = 'cone_specific';                         
        tParams.surroundSize = 10;                             
        tParams.surroundWeight = 0.7;                        
        tParams.integrationArea = 0;                            
        tParams.secondSiteFanoFactor = 4;

    case 'macularpigmentvary'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        tParams.OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; 
        tParams.DO_TAFC_CLASSIFIER_STATES = true;             
        tParams.macularPigmentDensityAdjustments = [-0.3 0 0.3];
        
        tParams.noiseType = 1;
        tParams.surroundType = 'none';                          
        tParams.surroundSize = 0;                              
        tParams.surroundWeight = 0;                             
        tParams.integrationArea = 0;                           
        tParams.secondSiteFanoFactor = 0;
    
    case 'quicktest'
        sParams.stimulus.coneNumbersToUse = [4 2 1]; 

        sParams.nColorDirections = 16;
        sParams.dirAngleMax = 2*pi;
        sParams.nTestLevels = 8;
        sParams.nDrawsPerTestStimulus = 100;
        
        tParams.OBSERVER_STATES = {'LMandS'};
        tParams.DO_TAFC_CLASSIFIER_STATES = true;
        tParams.macularPigmentDensityAdjustments = 0;
        
        tParams.noiseType = 1;
        tParams.surroundType = 'none';                        
        tParams.surroundSize = 0;                              
        tParams.surroundWeight = 0;                            
        tParams.integrationArea = 0;                         
        tParams.secondSiteFanoFactor = 1;                      
                                                       
    otherwise
        error('Unknown parameter presets');
end

end
