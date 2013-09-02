function colorContours(parameterPreset)
% colorContours([parameterPreset])
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
%  * Can do both YN and TAFC svm based observers.
%  * Has a simple surround model, along with simple second site noise.
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
%  * Think about how to add a real surround. The current implementation is a hack
%    that works (sort of, see below) for uniform fields but will not be right
%    for stimuli with spatial structure.
%  * Shouldn't fit ellipses to dichromatic data -- want pairs of lines.  And 
%    I don't think the ellipse fitting enforces a center of zero.  It should.
%    At present I just turned the ellipse plotting off, because I found the
%    ellipses fit to the dichromats distracting.
%
% Known bugs and mysteries:
%  * Classification performance is above chance for some dichromatic conditions
%    where it shouldn't be.  This does require out of gamut contrast, and could
%    be a numerical issue, with nominal isolation not quite being isolation.  But it
%    is bugging me, a little.  I set up the contour plots so they don't plot out
%    of gamut data, but you can see this effect in the saved psychometric functions.
%  * I would expect a bit more dichromat advantage than I see, when SVM is at the cones,
%    simply because the dichromat has more cones for the cone class it has.  We see
%    an appropriate effect for deuteranopes but not for protanopes in the 4:2:1 case.
%    See presentations/EarlyResultsSeptember2013.pptx.  This might just be a law of
%    small numbers or visualization thing.  Could be explored by varying cone numbers used.
%  * The dichromatic sensitivity for cone-specific and cone random dichromats looked like
%    it comes out a bit differently (again see presentations/EarlyResultsSeptember2013.pptx.)
%    I think this is for a known reason in the way dichromats are implemented: we replace
%    the spectral sensitivity of the missing cone by that of the present cone class, but keep
%    the nominal cone labels.  Thus the surrounds for cone random and cone specific surrounds
%    are drawn from pools of different sizes, and in the present implementation those pools
%    are rather small.  So there could be differences in noise averaging properties for the
%    two types of surround.  I think this issue will go away when we implement real surrounds
%    that follow the spatial structure of the mosaic, rather than the current special case
%    that only makes sense for uniform fields.
%
%  Some specific and minor things to patch up are indicated with comments
%  starting with [**] below, where they apply.
%
% Requires:
%   isetbio
%     - Available on gitHub as https://github.com/wandell/isetbio.git
%   PsychophysicsToolbox
%     - Available on gitHub as https://github.com/Psychtoolbox-3/Psychtoolbox-3.git
%
% 8/2/13  dhb/bw Our excellent adventure commences.
% 8/3/13  dhb/bw Tune up and add classifier.
% 8/4/13  dhb/bw Check irradiance calcs between ISET and PTB.
% 8/6/13  dhb/bw Got contour going.  End of this DHB visit west.
% 8/11/13 dhb    Added ellipse fitting, macular pigment adjustment
%         dhb    Added TAFC option.  Not fully tested.
% 8/12/13 dhb    Added dichromats, not tested at all.
% 8/13/13 dhb    Sometimes Matlab's svmtrain lives in the bioinfo toolbox.  Remove that too.
%         dhb    A little work on memory management.  Tweak params to leave running overnight.
% 8/16/13 dhb    Working on parallization.  In a broken state right now but gotta run.
% 8/18/13 dhb    I think cluster stuff is working now.
%         dhb    Surround, second site noise support.
% 8/24/13 dhb    Don't fit/plot thresholds beyond max test level measured.
% 9/2/13  dhb    This is now in pretty good shape.  Good enough that we might learn things from it.

%% Can compute only, analyze only, or do both
%
% Generally do both unless analysis changes without need
% to do the long recompute.
COMPUTE = true;                                % Compute?
ANALYZE = true;                                % Analyze

%% Control diagnostics
%
% These are useful when debugging but messy
% when things are working.  These are set
% in their own structure so they are controlled
% by how they are set here, rather than overridden
% at analysis time by loading in the other saved
% parameter structures.
runtimeParams.DO_SIM_PLOTS = false;
runtimeParams.SIM_QUIET = true;
runtimeParams.DO_PSYCHO_PLOTS = true;
runtimeParams.psychoPlotDir = 'psychometricFcnPlots';
runtimeParams.theContourPlotLim = 0.2;
runtimeParams.plotEllipses = false;

%% Set up parameters
if (nargin < 1 || isempty(parameterPreset))
    parameterPreset = 'QuickTest';
end
[theParams,staticParams] = setParameters(parameterPreset);

%% Set up output directory name.
%
% Try to make name tell us a lot about conditions,
% so that we can keep separate what happened in different runs.
%
% This needs to be set here so that we can load in computed parameters
% from the right file, once we've run a big calculation and just want
% to analyze the saved data.
if (staticParams.dirAngleMax == 2*pi)
    runtimeParams.outputDir = sprintf('%s%s_fullCircle_%d_%d_%d_%s_%d_%d_%d_%d',...
        staticParams.outputRoot,staticParams.parameterPreset,staticParams.nColorDirections,staticParams.nTestLevels,staticParams.nDrawsPerTestStimulus,...
        theParams.surroundType,round(100*theParams.surroundSize),round(100*theParams.surroundWeight),...
        round(100*theParams.integrationArea),round(100*theParams.secondSiteFanoFactor));
else
    runtimeParams.outputDir = sprintf('%s%s_halfCircle_%d_%d_%d__%s_%d_%d_%d_%d',...
        staticParams.outputRoot,staticParams.parameterPreset,staticParams.nColorDirections,staticParams.nTestLevels,staticParams.nDrawsPerTestStimulus,...
        theParams.surroundType,round(100*theParams.surroundSize),round(100*theParams.surroundWeight),...
        round(100*theParams.integrationArea),round(100*theParams.secondSiteFanoFactor));
end

%% Initialization for running n the cluster
%
% The whole program is in a big try/catch conditional, so we
% can close up the matlabpool on error as necessary.
try
    %% Make sure we are in our directory.  This does not
    % happen automatically when launched on the cluster.
    cd(fileparts(mfilename('fullpath'))); %#ok<MCCD>
    
    %% Open Matlab pool if it hasn't been opened for us.
    %
    % If we didn't open it, we don't close it.  You need
    % to have a function called IsCluster on your path
    % when you're running on the cluster and have it return
    % true to let this program know it is running on a DCE
    % cluster.
    if (COMPUTE & exist('IsCluster','file') && IsCluster && exist('matlabpool','file') && matlabpool('size') == 0)
        desiredPoolSize = 15;
        matlabpool(desiredPoolSize);
        NEEDTOCLOSEPOOL = 1;
    else
        NEEDTOCLOSEPOOL = 0;
    end
    
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
        
        %% Make sure random number generator seed is different each run.
        staticComputedValues.ranSeed = ClockRandSeed;
        
        %%  Get display spectra
        d = displayCreate(staticParams.stimulus.monitorName);
        staticComputedValues.displaySpd = displayGet(d,'spd');
        staticComputedValues.wavelengthsNm = displayGet(d,'wave');
        % vcNewGraphWin; plot(w,displaySpd)
        
        %% Set background spectrum
        staticComputedValues.backRGB = [staticParams.stimulus.backRGBValue staticParams.stimulus.backRGBValue staticParams.stimulus.backRGBValue]';
        staticComputedValues.backSpd = staticComputedValues.displaySpd*staticComputedValues.backRGB;

        %% Create standard human polychromatic PSF
        % using wavefront tools, and make it an
        % iset OI thingy.
        wvf = wvfCreate('wave',staticComputedValues.wavelengthsNm);
        sample_mean = wvfLoadThibosVirtualEyes(staticParams.pupilDiameterMm);
        wvf = wvfSet(wvf,'zcoeffs',sample_mean);
        wvf = wvfComputePSF(wvf);
        oiD = wvf2oi(wvf,'shift invariant');
        optics = oiGet(oiD,'optics');
        staticComputedValues.focalLengthMm = opticsGet(optics,'focal length','mm');
        vcAddAndSelectObject(oiD);
        clear wvf
        % oiWindow;
        % vcNewGraphWin; plotOI(oiD,'psf')
        
        %% Store just what we need
        %
        % We only keep one copy of these, because they are big
        staticComputedValues.oiD = oiD;
        clear optics oiD
        
        simParams = constructSimulationParameters(theParams,staticParams);
        nParams = length(simParams);  
        
        %% Make/clear output directory
        if (~exist(runtimeParams.outputDir,'dir'))
            mkdir(runtimeParams.outputDir);
        else
            unix(['rm -rf ' fullfile(runtimeParams.outputDir,'*') ';']);
        end
        
        %% Loop over all the simulations in one big parfor loop.
        %
        % This is the long slow part.
        if (exist('IsCluster','file') && IsCluster)
            mkdir(fullfile(runtimeParams.outputDir,'clusterLogFiles',''));
            parfor p = 1:nParams
                simResults(p) = doOneSimulation(simParams(p),staticParams,runtimeParams,staticComputedValues);
                
                % Write a little log file so we can track what's happening from afar
                fid = fopen(fullfile(runtimeParams.outputDir,'clusterLogFiles',['done.' num2str(p) '_' num2str(nParams)]),'wt');
                fprintf(fid,'\tSimulation %d of %d\n',p,nParams);
                fprintf(fid,'\tCalculations for observer state %s\n',simParams(p).OBSERVER_STATE);
                fprintf(fid,'\tTAFC state %d\n',simParams(p).DO_TAFC_CLASSIFIER);
                fprintf(fid,'\tMacular pigment density adjust %0.2f\n',simParams(p).macularPigmentDensityAdjust);
                fprintf(fid,'\tColor direction %0.3f\n',simParams(p).cdAngle);
                for k = 1:length(simResults(p).fractionCorrect)
                    fprintf(fid,'\t\tTest level %0.3f, fraction correct %0.2f\n',simResults(p).testLevels(k),simResults(p).fractionCorrect(k));
                end
                fclose(fid);
            end
        else
            for p = 1:nParams
                fprintf('\n\tSimulation %d of %d\n',p,nParams);
                fprintf('\tCalculations for observer state %s\n',simParams(p).OBSERVER_STATE);
                fprintf('\tTAFC state %d\n',simParams(p).DO_TAFC_CLASSIFIER);
                fprintf('\tMacular pigment density adjust %0.2f\n',simParams(p).macularPigmentDensityAdjust);
                fprintf('\tColor direction %0.3f\n',simParams(p).cdAngle);
                
                simResults(p) = doOneSimulation(simParams(p),staticParams,runtimeParams,staticComputedValues);
                
                for k = 1:length(simResults(p).fractionCorrect)
                    fprintf('\t\tTest level %0.3f, fraction correct %0.2f\n',simResults(p).testLevels(k),simResults(p).fractionCorrect(k));
                end
            end
        end
        
        %% Save the results
        %
        % This lets us reload
        % and analyze/plot away from the cluster.
        save(fullfile(runtimeParams.outputDir,'simResults.mat'),'theParams','staticParams','staticComputedValues','simParams','simResults');
    end
    
    %% **************
    % Analyze Section
    %% **************
    if (ANALYZE)
        %% Load
        theData = load(fullfile(runtimeParams.outputDir,'simResults'),'theParams','staticParams','staticComputedValues','simParams','simResults');
        
        %% Compare what we load with what was set at the top, to make sure we are analyzing what we think we are.
        %
        % [**] The checks go here, when we get around to writing them.
        
        % Then clear to make sure we are analyzing the data we loaded  
        clear theParams staticParams staticComputedValues simParams simResults

  
        %% Figure out what was run
        %
        % If all is working right, we already have this from the parameters
        % at the top, but it seems wise to recreate from the data.
        %
        % [**] Could implement a check that these match what we think they should
        % be, given what is in the loaded structure theParams.
        OBSERVER_STATES_LIST = {theData.simParams.OBSERVER_STATE};
        THE_OBSERVER_STATES = unique(OBSERVER_STATES_LIST);
        DO_TAFC_CLASSIFIER_STATES_LIST = [theData.simParams.DO_TAFC_CLASSIFIER];
        THE_DO_TAFC_CLASSIFIER_STATES = unique(DO_TAFC_CLASSIFIER_STATES_LIST);
        macularPigmentDensityAdjustments_List = [theData.simParams.macularPigmentDensityAdjust];
        the_macularPigmentDensityAdjustments = unique(macularPigmentDensityAdjustments_List);
        
        %% Make psychometric function output dir, if necessary
        if (runtimeParams.DO_PSYCHO_PLOTS)
            if (~exist(fullfile(runtimeParams.outputDir,runtimeParams.psychoPlotDir,''),'file'))
                mkdir(fullfile(runtimeParams.outputDir,runtimeParams.psychoPlotDir,''));
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
                    useParams0 = theData.simParams(index0);
                    useResults0 = theData.simResults(index0);
                    
                    cdAngles_List = [useParams0.cdAngle];
                    the_cdAngles = unique(cdAngles_List);
                    for cdi = 1:length(the_cdAngles);
                        cdAngle = the_cdAngles(cdi);
                        index1 = find(cdAngle == cdAngles_List);
                        if (length(index1) ~= 1)
                            error('We only expect result for each angle for each condition');
                        end
                        useParams1 = useParams0(index1);
                        useResults1 = useResults0(index1);
                        
                        for t = 1:length(useResults1.testLevels)
                            the_testLevels(t) = useResults1.testLevels(t);
                            the_fractionCorrects(t) = useResults1.fractionCorrect(t);
                            the_nCorrectResponses(t) = useResults1.nCorrectResponses(t);
                            the_nTotalResponses(t) = useResults1.nTotalResponses(t);
                            the_testLMSContrast(:,t) = useResults1.testLMSContrast;
                            the_backLMS(:,t) = useResults1.backgroundLMS;
                            the_testLMSGamut(:,t) = useResults1.testLMSGamut;
                        end
                        [nil,sindex] = sort(the_testLevels);
                        the_testLevels = the_testLevels(sindex);
                        the_fractionCorrects = the_fractionCorrects(sindex);
                        the_nCorrectResponses = the_nCorrectResponses(sindex);
                        the_nTotalResponses = the_nTotalResponses(sindex);
                        the_testLMSContrast = the_testLMSContrast(:,sindex);
                        the_backLMS = the_backLMS(:,sindex);
                        the_testLMSGamut = the_testLMSGamut(:,sindex);
                        
                        % Plot data
                        if (runtimeParams.DO_PSYCHO_PLOTS)
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
                        testLevelsInterp = linspace(0,the_testLevels(end),100);
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
                            thresholdEst = PFI(paramsValues,theData.staticParams.criterionCorrect);
                        else
                            [alpha,beta] = FitWeibTAFC(the_testLevels,the_nCorrectResponses,the_nTotalResponses-the_nCorrectResponses,[],1/2);
                            thresholdEst = FindThreshWeibTAFC(theData.staticParams.criterionCorrect,alpha,beta);
                            probCorrInterp = ComputeWeibTAFC(testLevelsInterp,alpha,beta);
                        end
                        
                        % Print threshold
                        if (~runtimeParams.SIM_QUIET)
                            fprintf('%d%% correct threshold is %0.1f\n',round(100*theData.staticParams.criterionCorrect),thresholdEst);
                        end
                        
                        % Finish plot
                        if (runtimeParams.DO_PSYCHO_PLOTS)
                            plot(testLevelsInterp,probCorrInterp,'r');
                            plot([thresholdEst thresholdEst],[0.5 theData.staticParams.criterionCorrect],'g');
                            plot([0 thresholdEst],[theData.staticParams.criterionCorrect theData.staticParams.criterionCorrect],'g');
                            xlim([0 2]);
                            ylim([0.5 1]);
                            drawnow;
                            
                            % Save individual psychometric function figures, mainly for debugging level choice
                            if (DO_TAFC_CLASSIFIER_STATE)
                                outName = sprintf('psychoFig_%s_TAFC_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),round(1000*cdAngle));
                            else
                                outName = sprintf('psychoFig_%s_YN_%d_%d',OBSERVER_STATE,round(100*macularPigmentDensityAdjust),round(1000*cdAngle));
                            end
                            set(gca,'LooseInset',get(gca,'TightInset'));
                            saveas(psychoFig,fullfile(runtimeParams.outputDir,runtimeParams.psychoPlotDir,outName),'png');
                        end
                        
                        % Store results for this color direction
                        contourThreshResults(cdi).maxTestLevel = max(the_testLevels);
                        contourThreshResults(cdi).thresholdLevel = thresholdEst;
                        contourThreshResults(cdi).testLMSContrast = the_testLMSContrast;
                        contourThreshResults(cdi).thresholdLMSContrast = thresholdEst*the_testLMSContrast;
                        contourThreshResults(cdi).backgroundLMS = the_backLMS;
                        contourThreshResults(cdi).testLMSGamut = the_testLMSGamut;
                    end
                    
                    % Collect up thresholds for fitting.
                    %
                    % Only take those where threshold was inside of gamut
                    % This is always defined as a level of 1.
                    LContourPoints = [];
                    MContourPoints = [];
                    for cdi = 1:theData.staticParams.nColorDirections
                        if (contourThreshResults(cdi).thresholdLevel < 1)
                            LContourPoints = [LContourPoints ; contourThreshResults(cdi).thresholdLMSContrast(1)];
                            MContourPoints = [MContourPoints ; contourThreshResults(cdi).thresholdLMSContrast(2)];
                        end
                    end
                    
                    % Reflect around origin if directions only sampled around hemicircle
                    if (theData.staticParams.dirAngleMax == pi)
                        LContourPoints = [LContourPoints ; -LContourPoints];
                        MContourPoints = [MContourPoints ; -MContourPoints];
                    end
                    
                    %% Make a plot of the threshold contour, and fit it.
                    contourFig = figure; clf; hold on
                    plot(LContourPoints,MContourPoints,'ro','MarkerFaceColor','r','MarkerSize',8);
                    if (length(LContourPoints) > 6)
                        % Sometimes the ellipse fitting routine throws an error if the data aren't close enough
                        % to an ellipse.  We could diagnose this more, but it generally happens only for cases
                        % where the number of parameters is set small to test something, and the try/catch keeps
                        % the program from crashing out.
                        try
                            [ellipseZ, ellipseA, ellipseB, ellipseAlpha] = fitellipse([LContourPoints' ; MContourPoints']);
                            if (runtimeParams.plotEllipses)
                                plotellipse(ellipseZ,ellipseA,ellipseB,ellipseAlpha,'r');
                            end
                        catch
                            fprintf('Ellipse fit failed, skipping and moving on\n');
                        end
                    end
                  
                    plot([-runtimeParams.theContourPlotLim runtimeParams.theContourPlotLim],[0 0],'k:');
                    plot([0 0],[-runtimeParams.theContourPlotLim runtimeParams.theContourPlotLim],'k:');
                    xlim([-runtimeParams.theContourPlotLim runtimeParams.theContourPlotLim]);
                    ylim([-runtimeParams.theContourPlotLim runtimeParams.theContourPlotLim]);
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
                    saveas(contourFig,fullfile(runtimeParams.outputDir,outName),'png');
                    
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




