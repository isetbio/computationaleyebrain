%% s_ccQuickTest(parameterPreset)
%
% Illustrate the steps in producing a color threshold contour using
% ISETBIO.
%
% Requires:
%   ISETBIO
%     - Available on gitHub as https://github.com/isetbio/isetbio.git
%   PsychophysicsToolbox-3
%     - Available on gitHub as https://github.com/Psychtoolbox-3/Psychtoolbox-3.git
%
% We set up a background (gray).
% We create a four color test directions
% We use ISETBIO and PTB to create isomerizations to create
%                background and background + s*test
% for some scalar, s.
% We classify as a function of s for each test direction.
% We plot the psychometric functions and the detection contours
%
% DHB/BW/HJ (c) ISETBIO Team, 2013

%%
s_initISET

%% Main parameters for background
monitorName = 'LCD-Apple.mat';
wave   = 380:4:780;
bLevel = 0.5;             % Linear value of display primaries

%% Make the background file
d       = displayCreate(monitorName);
gTable  = displayGet(d,'gamma table');  % plot(gTable)
igTable = ieLUTInvert(gTable);          % Maps linear values to DAC
bDAC    = ieLUTLinear(repmat(bLevel,1,3),igTable)/size(gTable,1);
bImage = ones(128,128,3);
for ii=1:3
    bImage(:,:,ii) = bImage(:,:,ii)*bDAC(ii);
end
bFile = fullfile(isetbioRootPath,'tmp','bFile.png');
imwrite(bImage,bFile);

%% Initiate the human optics
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);
wvf    = wvfComputePSF(wvf);
oiB    = wvf2oi(wvf,'shift invariant');

%% Build the background scene and oi from the image file.
bScene = sceneFromFile(bFile,'rgb',[],monitorName,wave);
oiB    = oiCompute(oiB,bScene);

% vcAddAndSelectObject(oiB); oiWindow;
% vcNewGraphWin; plotOI(oiB,'psf')

%% Run for one test light

simResults(p) = doOneSimulation(simParams(p),staticParams,runtimeParams,staticComputedValues);




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
                %   nSD = 1.5;
                %   [eVec,h,ptsAndCrv] = covEllipsoid([LContourPoints(:), MContourPoints(:)],nSD, vcNewGraphWin)
                %
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

    
