% colorContours
%
% This script will eventually compute predicted ideal observer
% isodetection contours for various psychophyiscal discrimination
% experiments.
%
% At present, this is Phase 1 (in which Doris gets her oats).
%
% Steps currently implemented:
%  * Create isetbio scenes based on realistic display primaries.
%  * Model the transformations from stimulus to cone absorptions.
%  * Specify a direction in color space and use isetbio machinery
%    to compute distributions of mosaic responses for signal and blank
%    trials at various levels of steps in the target color direction.
%  * Use svmlib classifier to predict percent correct responses at
%    each level.
%  * Fit predicted psychometric function and find psychophysical threshold.
%  * Loop over different observer types and macular pigment density adjustments.
%  * can do both YN and TAFC svm based observers.
%
% To do:
%  * Check cone absorption computations against PTB calculations.
%     - ISET works by taking the cone quantal efficiences computed elsewhere
%       and using these directly.
%     - This code is currently sticking the PTB computed versions into
%       the ISET structure.
%     - With that, there are slight differences between PTB and iset
%       computation, due to different conversion factor for computing
%       irradiance from radiance.  The differences are small and stem
%       from the fact that iset takes magnification and paraxial rays
%       into account, whereas PTB does not.
%     - We want to compare Hiroshi's isomerization code with PTBs.
%     - Neither set of code incorporates Stiles-Crawford effect.
%  * Add more sensible code to control spatial integration area.
%  * Add eye movements.
%  * Add simple foveal midget ganglion cell model.  This would be
%    pretty easy if we continue to ignore spatial structure.
%  * Break this big long script into sensible subfunctions.
%  * Clusterize.
%  * Shouldn't fit ellipses to dichromatic data -- want pairs of lines.
%  * I don't think the ellipse fitting enforces a center of zero.  It should.
%  * Need better system for tracking output.  Currently filenames tell you
%    something but not enough.  This is going to start to matter soon.  May
%    also want to save data as well as plots.
%
% Known bugs:
%  * 8/14/13 - Runs out of memory when run on Penn cluster head node.  Perhaps
%              this isn't really a problem since we should only debug on head node.

%  Some specific and minor things to patch up are indicated with comments
%  starting with [**] below, where they apply.
%
% Requires:
%   isetbio
%     - Available on gitHub as https://github.com/wandell/isetbio.git
%   PsychophysicsToolbox
%     - Available on gitHub as https://github.com/Psychtoolbox-3/Psychtoolbox-3.git
%
% 8/2/13  DHB/BW Our excellent adventure commences.
% 8/3/13  DHB/BW Tune up and add classifier.
% 8/4/13  DHB/BW Check irradiance calcs between ISET and PTB.
% 8/6/13  DHB/BW Got contour going.  End of this DHB visit west.
% 8/11/13 DHB    Added ellipse fitting, macular pigment adjustment
%         DHB    Added TAFC option.  Not fully tested.
% 8/12/13 DHB    Added dichromats, not tested at all.
% 8/13/13 DHB    Sometimes Matlab's svmtrain lives in the bioinfo toolbox.  Remove that too.
%         DHB    A little work on memory management.  Tweak params to leave running overnight.
% 8/16/13 DHB    Working on parallization.  In a broken state right now but gotta run.

%% Clear out the junk.  Remember where you are.
%
% Sometime, but not always s_initISET clears debugger stop points.
% So I (DHB) commented it out pending a better understanding.
clear; close all;  %s_initISET
talkD = pwd;
saveFlag = 0;

%% Parameters
integrationTimeSecs = 0.050;                    % Temporal integration time for detecting mechanisms.
fieldOfViewDegrees = 2;                         % Field of view specified for the scenes.
scenePixels = 64;                               % Size of scenes in pixels

monitorName = 'LCD-Apple';                      % Monitor spectrum comes from this file
backRGBValue = 0.5;                             % Define background for experment in monitor RGB

pupilDiameterMm = 3;                            % Pupil diameter.  Used explicitly in the PSF calc.
% Need to carry this through to the absorption calculations.

coneProportions = [0.1 .6 .2 .1];               % Proportions of cone types in the mosaic, order: empty, L,M,S
coneApertureMeters = [sqrt(4.1) sqrt(4.1)]*1e-6;% Size of (rectangular) cone apertures, in meters.
% The choice of 4.1 matches the area of a 2.3 micron diameter IS diameter,
% and that is PTB's default.
isetSensorConeSlots = [2 3 4];                  % Indices for LMS cones in iset sensor returns.   These run 2-4 because
% of the empty pixels
nSensorClasses = length(isetSensorConeSlots);   % For convenience, specify the number of sensor classes.

nColorDirections = 16;                          % Number of color directions for contour.
dirAngleMax = 2*pi;                             % Use pi for sampling directions from hemicircle, 2*pi for whole circle

nTestLevels = 8;                                % Number of test levels to simulate for each test direction psychometric function.
nDrawsPerTestStimulus = 400;                    % Number of noise draws used in the simulations, per test stimulus
noiseType = 1;                                  % Noise type passed to isetbio routines.  1 -> Poisson.

criterionCorrect = 0.82;                        % Fraction correct for definition of threshold in TAFC simulations.
testContrastLengthMax = 0.5;                    % Default maximum contrast lenght of test color vectors used in each color direction.
% Setting this helps make the sampling of the psychometric functions more efficient.
% This value can be overridden in a switch statement on OBSERVER_STATE in a loop below.

outputRoot = 'output';                          % Plots get dumped a directory with this root name, but additoinal
% characters to identify parameters.
psychoPlotDir = 'psychometricFcnPlots';

surroundSize = 0;                               % Parameter defining surround size.  NOT YET IMPLEMENTED.
surroundWeight = 0;                             % Parameter defining surround weight.  NOT YET IMPLEMENTED.
integrationArea = 0;                            % Stimulus integration area.  NOT YET IMPLEMENTED.

macularPigmentDensityAdjustments = [-0.3 0 0.3]; % Amount to adjust macular pigment density for cone fundamentals of simulated observer.
% Note that stimuli are computed for a nominal (no adjustment) observer.
DO_TAFC_CLASSIFIER_STATES = [true false];       % Can be true, false, or [true false]
OBSERVER_STATES = {'MSonly' 'LSonly'};          % Simulate various tri and dichromats

QUICK_TEST_PARAMS = true;                       % Set to true to override parameters with a small number of trials for debugging.

%% Process quick test option
if (QUICK_TEST_PARAMS)
    nColorDirections = 4;
    dirAngleMax = pi;
    nTestLevels = 4;
    nDrawsPerTestStimulus = 100;
    macularPigmentDensityAdjustments = [-0.3 0 0.3];
    DO_TAFC_CLASSIFIER_STATES = [true false];
end

%% Make output directories if they doesn't exist.
%
% Try to make name tell us a lot about static conditions.
if (dirAngleMax == 2*pi)
    outputDir = sprintf('%s_fullCircle_%d_%d_%d  ',outputRoot,nColorDirections,nTestLevels,nDrawsPerTestStimulus,...
        round(100*surroundSize),round(100*surroundWeight),round(100*integrationArea));
else
    outputDir = sprintf('%s_halfCircle_%d_%d_%d  ',outputRoot,nColorDirections,nTestLevels,nDrawsPerTestStimulus,...
        round(100*surroundSize),round(100*surroundWeight),round(100*integrationArea));
end

if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
if (~exist(fullfile(outputDir,psychoPlotDir,''),'file'))
    mkdir(fullfile(outputDir,psychoPlotDir,''));
end

%% Make sure random number generator seed is different each run.
ClockRandSeed;

%% Get stats path off of the path, temporarily.
%
% This avoids a name space conflict with the libsvn
% routines that we use here.
if (exist('RemoveMatchingPaths','file'))
    path(RemoveMatchingPaths(path,'stats'));
    path(RemoveMatchingPaths(path,'bioinfo'));
end


%%  Get display spectra
d = displayCreate(monitorName);
displaySpd = displayGet(d,'spd');
wavelengthsNm = displayGet(d,'wave');
% vcNewGraphWin; plot(w,displaySpd)

%% Set background
backRGB = [backRGBValue backRGBValue backRGBValue]';
backSpd = displaySpd*backRGB;

%% Create a scene with the background spectrum
%
% To patch into scene creation routines, we create an image
% on disk with spatially uniform RGB values.  These are then treated
% as frame buffer values, and will be raised to the 2.2 (gamma uncorrected)
% by the scene creation routine.
%
% [**] May want to find a way to do this (and the corresponding version
% for the test) that does not involve writing an image to disk and
% that skips the gamma correction stuff.
gammaValue = 2.2;
backImg = ones(scenePixels,scenePixels,3)*(backRGBValue)^(1/gammaValue);
imwrite(backImg,'backFile.png','png');
sceneB = sceneFromFile('backFile.png','rgb',[],'LCD-Apple.mat',wavelengthsNm);
sceneB = sceneSet(sceneB,'name','background');
sceneB = sceneSet(sceneB,'fov',fieldOfViewDegrees);
vcAddAndSelectObject(sceneB);
%sceneWindow;

%% Create standard human polychromatic PSF
% using wavefront tools, and make it an
% iset OI thingy.
wvf = wvfCreate('wave',wavelengthsNm);
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf = wvfSet(wvf,'zcoeffs',sample_mean);
wvf = wvfComputePSF(wvf);
oiD = wvf2oi(wvf,'shift invariant');
optics = oiGet(oiD,'optics');
focalLengthMm = opticsGet(optics,'focal length','mm');
vcAddAndSelectObject(oiD);
clear wvf
% oiWindow;
% vcNewGraphWin; plotOI(oiD,'psf')

%% Construct list of conditions
%
% These are strung out so that we can chunk through
% them in a big parfor loop below.
cdAngles = linspace(0,dirAngleMax,nColorDirections+1);
cdAngles = cdAngles(1:end-1);
testLevels = linspace(0,1,nTestLevels);
paramIndex = 1;
for os = 1:length(OBSERVER_STATES)
    for ct = 1:length(DO_TAFC_CLASSIFIER_STATES)
        for m = 1:length(macularPigmentDensityAdjustments)
            for cdi = 1:nColorDirections
                for t = 1:length(testLevels)
                    params(paramIndex).OBSERVER_STATE = OBSERVER_STATES{os};
                    params(paramIndex).DO_TAFC_CLASSIFIER = DO_TAFC_CLASSIFIER_STATES(ct);
                    params(paramIndex).macularPigmentDensityAdjust = macularPigmentDensityAdjustments(m);
                    params(paramIndex).cdAngle = cdAngles(cdi);
                    params(paramIndex).testLevel = testLevels(t);
                    
                    % Set test contrast maximum length.
                    %
                    % The best values depend on observer state
                    % and are currently set manually based on
                    % experience.
                    switch (params(paramIndex).OBSERVER_STATE)
                        case 'LMandS'
                            params(paramIndex).testContrastLengthMax = 0.3;
                        case 'LSonly'
                            params(paramIndex).testContrastLengthMax = 1;
                        case 'MSonly'
                            params(paramIndex).testContrastLengthMax = 1;
                        otherwise
                            error('Unknown dichromat/trichromat type specified');
                    end
                    
                    % Nuisance parameters
                    
                    % Set fixed params.  These can be made variable by
                    % adding a loop here.
                    params(paramIndex).sceneB = sceneB;
                    params(paramIndex).gammaValue = gammaValue;
                    params(paramIndex).oiD = oiD;
                    params(paramIndex).wavelengthsNm = wavelengthsNm;
                    params(paramIndex).backSpd = backSpd;
                    params(paramIndex).displaySpd = displaySpd;
                    params(paramIndex).pupilDiameterMm = pupilDiameterMm;
                    params(paramIndex).focalLengthMm = focalLengthMm;
                    params(paramIndex).integrationTimeSecs = integrationTimeSecs;
                    params(paramIndex).coneProportions = coneProportions;
                    params(paramIndex).coneApertureMeters = coneApertureMeters;
                    params(paramIndex).fieldOfViewDegrees = fieldOfViewDegrees;
                    params(paramIndex).backRGB = backRGB;
                    params(paramIndex).isetSensorConeSlots = isetSensorConeSlots;
                    params(paramIndex).nSensorClasses = nSensorClasses;
                    params(paramIndex).scenePixels = scenePixels;
                    params(paramIndex).nDrawsPerTestStimulus = nDrawsPerTestStimulus;
                    params(paramIndex).noiseType = noiseType;
                    
                    % Control diagnostics
                    % 
                    % These are useful when debugging but messy
                    % when running in parallel.  We probably don't
                    % need this fine level of control.
                    params(paramIndex).PLOT_COMPARE_IRRADIANCE = 0;
                    params(paramIndex).PLOT_TRAINING_TEST = 0;
                    params(paramIndex).PLOT_COMPARE_CONEQE = 0;
                    params(paramIndex).PRINT_OUT_PHOTORECEPTORS = 0;
                    params(paramIndex).VERBOSE = 0;

                    % Kluge for now to select subregion of total cones
                    params(paramIndex).fractionUse = 0.001;


                    % Bump counter
                    paramIndex = paramIndex+1;
                end
                
            end
        end
    end
end
nParams = length(params);
clear sceneB optics oiD

%% Loop over all the simulations in one big parfor loop.
if (exist('IsCluster','file') && IsCluster)
    parfor p = 1:nParams
        params(p).results = DoOneSimulation(params(p));      
    end
else
    for p = 1:nParams
        fprintf('\tSimulation %d of %d\n',p,nParams);
        fprintf('\tCalculations for observer state %s\n',params(p).OBSERVER_STATE);
        fprintf('\tTAFC state %d\n',params(p).DO_TAFC_CLASSIFIER);
        fprintf('\tMacular pigment density adjust %0.2f\n',params(p).macularPigmentDensityAdjust);
        fprintf('\tColor direction %0.3f\n',params(p).cdAngle);
        fprintf('\tTest level %0.3f\n',params(p).testLevel);
        
        params(p).results = DoOneSimulation(params(p));
        %                    
        
    end
end

%% Plot and fit psychometric function, extract threshold
%
% Plot data
if (~exist('psychoFig','var'))
    psychoFig = figure; clf; hold on
else
    figure(psychoFig); clf; hold on
end
plot(testLevels,fractionCorrect,'ro','MarkerSize',10,'MarkerFaceColor','r');

% Fit and find threshold, and add fit to plot
%
% Be setting USE_PALAMEDES to 1, you can use the Palamedes toolbox routines
% to do the fitting.  But, you'll need to install the Palamedes toolbox to
% do so.  By default, we use fitting routines in the Psychtoolbox.
testLevelsInterp = linspace(testLevels(1),testLevels(end),100);
USE_PALAMEDES = 0;
if (USE_PALAMEDES)
    % Palamedes toolbox (for fitting psychometric functions)
    %     - Available at http://www.palamedestoolbox.org, and described in Kingdom
    %       & Prins, Psychophysics: A Practical Introduction.
    %     - DHB's lab runs a version of this toolbox with some local modifications.
    %       Not sure if the code here relies on those.  The Palamedes distribution has
    %       been updated since DHB snagged his copy, and it looks like some of the features
    %       he put in are now available in the distribution.
    PF = @PAL_Weibull;                  % Alternatives: PAL_Gumbel, PAL_Weibull, PAL_CumulativeNormal, PAL_HyperbolicSecant
    PFI = @PAL_inverseWeibull;
    paramsFree = [1 1 0 0];             % 1: free parameter, 0: fixed parameter
    paramsValues0 = [mean(testLevels) 1/2 0.5 0];
    options = optimset('fminsearch');   % Type help optimset
    options.TolFun = 1e-09;             % Increase required precision on LL
    options.Display = 'off';            % Suppress fminsearch messages
    lapseLimits = [0 1];                % Limit range for lambda
    [paramsValues] = PAL_PFML_Fit(...
        testLevels',nCorrectResponses',nTotalResponses', ...
        paramsValues0,paramsFree,PF,'searchOptions',options, ...
        'lapseLimits',lapseLimits);
    probCorrInterp = PF(paramsValues,testLevelsInterp);
    thresholdEst = PFI(paramsValues,criterionCorrect);
else
    [alpha,beta] = FitWeibTAFC(testLevels,nCorrectResponses,nTotalResponses-nCorrectResponses,[],1/2);
    thresholdEst = FindThreshWeibTAFC(criterionCorrect,alpha,beta);
    probCorrInterp = ComputeWeibTAFC(testLevelsInterp,alpha,beta);
end
plot(testLevelsInterp,probCorrInterp,'r');
plot([thresholdEst thresholdEst],[0.5 criterionCorrect],'g');
plot([testLevels(1) thresholdEst],[criterionCorrect criterionCorrect],'g');
xlim([0 1]);
ylim([0.5 1]);
drawnow;

% Save individual psychometric function figures, mainly for debugging level choice
if (DO_TAFC_CLASSIFIER)
    outName = sprintf('psychoFig_%s_TAFC_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),cd);
else
    outName = sprintf('psychoFig_%s_YN_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),cd);
end
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(psychoFig,fullfile(outputDir,psychoPlotDir,outName),'png');

% Print threshold
fprintf('%d%% correct threshold is %0.1f\n',round(100*criterionCorrect),thresholdEst);

% Store results for this color direction
results(cd).thresholdLevel = thresholdEst;
results(cd).testLMSContrast = testLMSContrast;
results(cd).thresholdLMSContrast = thresholdEst*testLMSContrast;
results(cd).backgroundLMS = backLMS;
results(cd).testLMSGamut = testLMSGamut;

%% Fit an ellipse to the thresholds
LContourPoints = zeros(nColorDirections,1);
MContourPoints = zeros(nColorDirections,1);
for cd = 1:nColorDirections
    LContourPoints(cd) = results(cd).thresholdLMSContrast(1);
    MContourPoints(cd) = results(cd).thresholdLMSContrast(2);
end
if (dirAngleMax == pi)
    LContourPoints = [LContourPoints ; -LContourPoints];
    MContourPoints = [MContourPoints ; -MContourPoints];
end

%% Make a plot of the thresholds
contourFig = figure; clf; hold on
theContourPlotLim = 0.2;
plot(LContourPoints,MContourPoints,'ro','MarkerFaceColor','r','MarkerSize',8);
if (length(LContourPoints) > 6)
    % Sometimes the ellipse fitting routine throws an error if the data aren't close enough
    % to an ellipse.  We could diagnose this more, but it generally happens only for cases
    % where the number of parameters is set small to test something, and the try/catch keeps
    % the program from crashing out.
    try
        [ellipseZ, ellipseA, ellipseB, ellipseAlpha] = fitellipse([LContourPoints' ; MContourPoints']);
        plotellipse(ellipseZ,ellipseA,ellipseB,ellipseAlpha,'r');
    catch
        fprintf('Ellipse fit failed, skipping and moving on\n');
    end
end
plot([-theContourPlotLim theContourPlotLim],[0 0],'k:');
plot([0 0],[-theContourPlotLim theContourPlotLim],'k:');
xlim([-theContourPlotLim theContourPlotLim]);
ylim([-theContourPlotLim theContourPlotLim]);
axis('square');
xlabel('Nominal L cone contrast');
ylabel('Nominal M cone contrast');
title('Ideal observer threshold contours');
if (DO_TAFC_CLASSIFIER)
    outName = sprintf('colorContour_%s_TAFC_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust));
else
    outName = sprintf('colorContour_%s_YN_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust));
end
drawnow;
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(contourFig,fullfile(outputDir,outName),'png');

% These files are sort of big, so don't always save.
%
% Probably could figure out what not to save and get
% a useful piece of the data if we wanted to.
%save(outName);

% Close windows
close all

% Clear some big variables and pack
clear backSensorVals testSensorVals




