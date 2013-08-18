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
%  * Written so the slow part will run in parallel on a DCE cluster.
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
% 8/18/13 DHB    I think cluster stuff is working now.

%% Clear out the junk.  Remember where you are.
%
% Sometime, but not always s_initISET clears debugger stop points.
% So I (DHB) commented it out pending a better understanding.
clear; close all;  %s_initISET
talkD = pwd;
saveFlag = 0;

%% Parameter section

integrationTimeSecs = 0.050;                    % Temporal integration time for detecting mechanisms.
fieldOfViewDegrees = 2;                         % Field of view specified for the scenes.
scenePixels = 64;                               % Size of scenes in pixels

monitorName = 'LCD-Apple';                      % Monitor spectrum comes from this file
backRGBValue = 0.5;                             % Define background for experment in monitor RGB

pupilDiameterMm = 3;                            % Pupil diameter.  Used explicitly in the PSF calc.
                                                % [** Need to check that this is carried through to 
                                                % the absorption calculations.  We might be using an isetbio
                                                % default rather than the value set here.]

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

outputRoot = 'output';                          % Plots get dumped a directory with this root name, but with additional
                                                % characters to identify parameters of the run tacked on below.
psychoPlotDir = 'psychometricFcnPlots';         % Subdir for dumping psychometric function plots.

OBSERVER_STATES = {'LMandS' 'MSonly' 'LSonly'}; % Simulate various tri and dichromats
DO_TAFC_CLASSIFIER_STATES = [true];             % Can be true, false, or [true false]
macularPigmentDensityAdjustments = [0];         % Amount to adjust macular pigment density for cone fundamentals of simulated observer.
                                                % Note that stimuli are computed for a nominal (no adjustment) observer.

surroundSize = 10;                              % Parameter defining surround size.  NOT YET IMPLEMENTED.
surroundWeight = 0.7;                           % Parameter defining surround weight.  NOT YET IMPLEMENTED.
integrationArea = 0;                            % Stimulus integration area.  NOT YET IMPLEMENTED.
opponentLevelNoiseSd = 0;                       % Noise added after opponent recombination.  NOT YET IMPLEMENTED.

                                               
QUICK_TEST_PARAMS = false;                      % Set to true to override parameters with a small number of trials for debugging.

COMPUTE = true;                                 % Compute?
ANALYZE = true;                                 % Analyze

%% Initialization for running on the cluster
%
% Make sure we are in our directory.  This does not
% happen automatically when launched on the cluster.
try
    cd(fileparts(mfilename('fullpath'))); %#ok<MCCD>
    
    % Open Matlab pool if it hasn't been opened for us.
    %
    % If we didn't open it, we don't close it.
    if (exist('IsCluster','file') && IsCluster && exist('matlabpool','file') && matlabpool('size') == 0)
        desiredPoolSize = 15;
        matlabpool(desiredPoolSize);
        NEEDTOCLOSEPOOL = 1;
    else
        NEEDTOCLOSEPOOL = 0;
    end
    
    %% Process quick test option
    if (QUICK_TEST_PARAMS)
        nColorDirections = 4;
        dirAngleMax = pi;
        nTestLevels = 4;
        nDrawsPerTestStimulus = 100;
        macularPigmentDensityAdjustments = [0];
        DO_TAFC_CLASSIFIER_STATES = [false];
        OBSERVER_STATES = {'LMandS'};
    end
    
    %% Make output directories if they doesn't exist.
    %
    % Try to make name tell us a lot about static conditions,
    % so that we can keep separate what happened in different runs.
    if (dirAngleMax == 2*pi)
        outputDir = sprintf('%s_fullCircle_%d_%d_%d_%d_%d_%d_%d',outputRoot,nColorDirections,nTestLevels,nDrawsPerTestStimulus,...
            round(100*surroundSize),round(100*surroundWeight),round(100*integrationArea),round(100*opponentLevelNoiseSd));
    else
        outputDir = sprintf('%s_halfCircle_%d_%d_%d_%d_%d_%d_%d',outputRoot,nColorDirections,nTestLevels,nDrawsPerTestStimulus,...
            round(100*surroundSize),round(100*surroundWeight),round(100*integrationArea),round(100*opponentLevelNoiseSd));
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
    
    %% **************
    % Compute Section
    %% **************
    if (COMPUTE)
        
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
                            params(paramIndex).gammaValue = gammaValue;
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
                            params(paramIndex).surroundSize = surroundSize;
                            params(paramIndex).surroundWeight = surroundWeight;
                            params(paramIndex).integrationArea = integrationArea;
                            params(paramIndex).opponentLevelNoiseSd = opponentLevelNoiseSd;
                            
                            % Kluge for now to select subregion of total cones
                            params(paramIndex).fractionUse = 0.005;
                            
                            % Control diagnostics
                            %
                            % These are useful when debugging but messy
                            % when running in parallel.  We probably don't
                            % need this fine level of control.
                            params(paramIndex).PLOT_COMPARE_IRRADIANCE = false;
                            params(paramIndex).PLOT_TRAINING_TEST = false;
                            params(paramIndex).PLOT_COMPARE_CONEQE = false;
                            params(paramIndex).PRINT_OUT_PHOTORECEPTORS = false;
                            params(paramIndex).VERBOSE = false;
                            params(paramIndex).SVM_QUIET = true;
                            params(paramIndex).DO_PSYCHO_PLOTS = false;
                            
                            
                            % Bump counter
                            paramIndex = paramIndex+1;
                        end
                        
                    end
                end
            end
        end
        nParams = length(params);
        
        %% Static params
        %
        % We only keep one copy of these, because they are big
        staticParams.oiD = oiD;
        staticParams.sceneB = sceneB;
        clear sceneB optics oiD
        
        %% Make/clear output directory
        if (~exist(outputDir,'dir'))
            mkdir(outputDir);
        else
            unix(['rm -rf ' fullfile(outputDir,'*') ';']);
        end
        
        %% Loop over all the simulations in one big parfor loop.
        %
        % This is the long slow part.
        if (exist('IsCluster','file') && IsCluster)
            mkdir(fullfile(outputDir,'clusterLogFiles',''));
            parfor p = 1:nParams
                simResults(p) = DoOneSimulation(params(p),staticParams);
                
                % Write a little log file so we can track what's happening from afar
                fid = fopen(fullfile(outputDir,'clusterLogFiles',['done.' num2str(p) '_' num2str(nParams)]),'wt');
                fprintf(fid,'\n\tSimulation %d of %d\n',p,nParams);
                fprintf(fid,'\tCalculations for observer state %s\n',params(p).OBSERVER_STATE);
                fprintf(fid,'\tTAFC state %d\n',params(p).DO_TAFC_CLASSIFIER);
                fprintf(fid,'\tMacular pigment density adjust %0.2f\n',params(p).macularPigmentDensityAdjust);
                fprintf(fid,'\tColor direction %0.3f\n',params(p).cdAngle);
                fprintf(fid,'\tTest level %0.3f\n',params(p).testLevel);
                fprintf(fid,'\tFraction correct %0.2f\n',simResults(p).fractionCorrect);
                fclose(fid);
            end
        else
            for p = 1:nParams
                fprintf('\n\tSimulation %d of %d\n',p,nParams);
                fprintf('\tCalculations for observer state %s\n',params(p).OBSERVER_STATE);
                fprintf('\tTAFC state %d\n',params(p).DO_TAFC_CLASSIFIER);
                fprintf('\tMacular pigment density adjust %0.2f\n',params(p).macularPigmentDensityAdjust);
                fprintf('\tColor direction %0.3f\n',params(p).cdAngle);
                fprintf('\tTest level %0.3f\n',params(p).testLevel);
                
                simResults(p) = DoOneSimulation(params(p),staticParams);
                fprintf('\tFraction correct %0.2f\n',simResults(p).fractionCorrect);
            end
        end
        
        %% Save the results
        %
        % This lets us reload
        % and analyze/plot away from the cluster.
        save(fullfile(outputDir,'simResults.mat'),'params','simResults');
    end
    
    %% **************
    % Analyze Section
    %% **************
    if (ANALYZE)
        %% Load
        theData = load(fullfile(outputDir,'simResults'),'params','simResults');
        
        %% Figure out what was run
        %
        % If all is working right, we already have this from the parameters
        % at the top, but it seems wise to recreate from the data.
        %
        % [** Could implement a check that these match what is listed at the top.]
        OBSERVER_STATES_LIST = {theData.params.OBSERVER_STATE};
        THE_OBSERVER_STATES = unique(OBSERVER_STATES_LIST);
        DO_TAFC_CLASSIFIER_STATES_LIST = [theData.params.DO_TAFC_CLASSIFIER];
        THE_DO_TAFC_CLASSIFIER_STATES = unique(DO_TAFC_CLASSIFIER_STATES_LIST);
        macularPigmentDensityAdjustments_List = [theData.params.macularPigmentDensityAdjust];
        the_macularPigmentDensityAdjustments = unique(macularPigmentDensityAdjustments_List);
        DO_PSYCHO_PLOTS = any([theData.params.DO_PSYCHO_PLOTS]);
        VERBOSE = any([theData.params.VERBOSE]);
        
        %% Make psychometric function output dir, if necessary
        if (DO_PSYCHO_PLOTS)
            if (~exist(fullfile(outputDir,psychoPlotDir,''),'file'))
                mkdir(fullfile(outputDir,psychoPlotDir,''));
            end
        end
        
        %% We loop over the outer variables and get a contour for each, and
        % make and save a plot.
        for os = 1:length(THE_OBSERVER_STATES)
            OBSERVER_STATE = THE_OBSERVER_STATES{os};
            for ct = 1:length(THE_DO_TAFC_CLASSIFIER_STATES)
                DO_TAFC_CLASSIFIER_STATE = THE_DO_TAFC_CLASSIFIER_STATES(ct);
                for m = 1:length(the_macularPigmentDensityAdjustments)
                    macularPigmentDensityAdjust = the_macularPigmentDensityAdjustments(m);
                    index0 = find(strcmp(OBSERVER_STATE,OBSERVER_STATES_LIST)  & ...
                        DO_TAFC_CLASSIFIER_STATE == DO_TAFC_CLASSIFIER_STATES_LIST & ...
                        macularPigmentDensityAdjust == macularPigmentDensityAdjustments_List);
                    useParams0 = theData.params(index0);
                    useResults0 = theData.simResults(index0);
                    
                    cdAngles_List = [useParams0.cdAngle];
                    the_cdAngles = unique(cdAngles_List);
                    for cdi = 1:length(the_cdAngles);
                        cdAngle = the_cdAngles(cdi);
                        index1 = find(cdAngle == cdAngles_List);
                        useParams1 = useParams0(index1);
                        useResults1 = useResults0(index1);
                        testLevels_List = [useParams1.testLevel];
                        the_testLevels = unique(testLevels_List);
                        
                        for t = 1:length(the_testLevels)
                            the_fractionCorrects(t) = useResults1(t).fractionCorrect;
                            the_nCorrectResponses(t) = useResults1(t).nCorrectResponses;
                            the_nTotalResponses(t) = useResults1(t).nTotalResponses;
                            the_testLMSContrast(:,t) = useResults1(t).testLMSContrast;
                            the_backLMS(:,t) = useResults1(t).backgroundLMS;
                            the_testLMSGamut(:,t) = useResults1(t).testLMSGamut;
                        end
                        
                        % Plot data
                        if (DO_PSYCHO_PLOTS)
                            if (~exist('psychoFig','var'))
                                psychoFig = figure; clf; hold on
                            else
                                figure(psychoFig); clf; hold on
                            end
                            plot(the_testLevels,the_fractionCorrects,'ro','MarkerSize',10,'MarkerFaceColor','r');
                        end
                        
                        % Fit and find threshold, and add fit to plot
                        %
                        % Be setting USE_PALAMEDES to 1, you can use the Palamedes toolbox routines
                        % to do the fitting.  But, you'll need to install the Palamedes toolbox to
                        % do so.  By default, we use fitting routines in the Psychtoolbox.
                        testLevelsInterp = linspace(the_testLevels(1),the_testLevels(end),100);
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
                            paramsValues0 = [mean(the_testLevels) 1/2 0.5 0];
                            options = optimset('fminsearch');   % Type help optimset
                            options.TolFun = 1e-09;             % Increase required precision on LL
                            options.Display = 'off';            % Suppress fminsearch messages
                            lapseLimits = [0 1];                % Limit range for lambda
                            [paramsValues] = PAL_PFML_Fit(...
                                the_testLevels',the_nCorrectResponses',the_nTotalResponses', ...
                                paramsValues0,paramsFree,PF,'searchOptions',options, ...
                                'lapseLimits',lapseLimits);
                            probCorrInterp = PF(paramsValues,testLevelsInterp);
                            thresholdEst = PFI(paramsValues,criterionCorrect);
                        else
                            [alpha,beta] = FitWeibTAFC(the_testLevels,the_nCorrectResponses,the_nTotalResponses-the_nCorrectResponses,[],1/2);
                            thresholdEst = FindThreshWeibTAFC(criterionCorrect,alpha,beta);
                            probCorrInterp = ComputeWeibTAFC(testLevelsInterp,alpha,beta);
                        end
                        
                        % Print threshold
                        if (VERBOSE)
                            fprintf('%d%% correct threshold is %0.1f\n',round(100*criterionCorrect),thresholdEst);
                        end
                        
                        % Finish plot
                        if (DO_PSYCHO_PLOTS)
                            plot(testLevelsInterp,probCorrInterp,'r');
                            plot([thresholdEst thresholdEst],[0.5 criterionCorrect],'g');
                            plot([the_testLevels(1) thresholdEst],[criterionCorrect criterionCorrect],'g');
                            xlim([0 1]);
                            ylim([0.5 1]);
                            drawnow;
                            
                            % Save individual psychometric function figures, mainly for debugging level choice
                            if (DO_TAFC_CLASSIFIER_STATE)
                                outName = sprintf('psychoFig_%s_TAFC_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),round(1000*cdAngle));
                            else
                                outName = sprintf('psychoFig_%s_YN_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),round(1000*cdAngle));
                            end
                            set(gca,'LooseInset',get(gca,'TightInset'));
                            saveas(psychoFig,fullfile(outputDir,psychoPlotDir,outName),'png');
                        end
                        
                        % Store results for this color direction
                        contourThreshResults(cdi).thresholdLevel = thresholdEst;
                        contourThreshResults(cdi).testLMSContrast = the_testLMSContrast;
                        contourThreshResults(cdi).thresholdLMSContrast = thresholdEst*the_testLMSContrast;
                        contourThreshResults(cdi).backgroundLMS = the_backLMS;
                        contourThreshResults(cdi).testLMSGamut = the_testLMSGamut;
                    end
                    
                    % Collect up thresholds for fitting.
                    nColorDirections = length(contourThreshResults);
                    LContourPoints = zeros(nColorDirections,1);
                    MContourPoints = zeros(nColorDirections,1);
                    for cdi = 1:nColorDirections
                        LContourPoints(cdi) = contourThreshResults(cdi).thresholdLMSContrast(1);
                        MContourPoints(cdi) = contourThreshResults(cdi).thresholdLMSContrast(2);
                    end
                    if (dirAngleMax == pi)
                        LContourPoints = [LContourPoints ; -LContourPoints];
                        MContourPoints = [MContourPoints ; -MContourPoints];
                    end
                    
                    %% Make a plot of the threshold contour, and fit it.
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
                    
                    % Save the plot
                    if (DO_TAFC_CLASSIFIER_STATE)
                        outName = sprintf('colorContour_%s_TAFC_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust));
                    else
                        outName = sprintf('colorContour_%s_YN_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust));
                    end
                    drawnow;
                    set(gca,'LooseInset',get(gca,'TightInset'));
                    saveas(contourFig,fullfile(outputDir,outName),'png');
                    
                    %% Close windows
                    %
                    % This can prevent Java heap overflow for big jobs.
                    close all
                end
            end
        end
        
 
    end
    
    %% Close Matlab pool
    if (NEEDTOCLOSEPOOL)
        matlabpool('close');
    end
    
%% Close Matlab pool on error
catch theErr
    if (NEEDTOCLOSEPOOL)
        matlabpool('close');
    end
    
    rethrow(theErr);
end




