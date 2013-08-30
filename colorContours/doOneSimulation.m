function results = doOneSimulation(theParams,staticParams,runtimeParams,staticComputedValues)
% results = doOneSimulation(theParams,staticParams,runtimeParams,staticComputedValues)
%
% Do one inner level of the simulations for colorContour.  Written so that
% it can be called in parallel.  See colorContours for comments about
% what this does.
%
% Returns a results structure with everything we want to know and save.
%
% TODO:
%  * Have not tested recently that something inside of various conditionals that control diagnostic plots/prints
%    is not broken.
%
% 8/16/13  dhb  Working on this.
% 8/18/13  dhb  Opponency and second site noise.
% 8/25/13  dhb  Parameter rationalization.

%% Use PTB to compute cone quantal sensitivities.
%
% We will use these in the isetbio scene structures and calculations.
% We compute two versions, one nominal and one with adjusted macular
% pigment densities.  The latter is used to compute the actual
% cone responses, the former to compute spectra for each nominal color
% direction.  This means that our plot in the end can simulate an observer
% with an adjusted macular pigment density.
%
% We also convert to sensivities in energy units.  The conversion
% to energy reads counter-intuitively -- the routine EnergyToQuanta
% is named for what it does to spectra; it does the inverse to
% sensitivities.

% First nominal sensitivities
[ptbNominalBackLMSIsomerizations,staticParams.pupilDiameterMm,ptbNominalPhotorceptorsStruct,ptbNominalIrradianceWattsPerM2] = ptbConeIsomerizationsFromSpectra(staticComputedValues.backSpd,staticComputedValues.wavelengthsNm,...
    staticParams.pupilDiameterMm,staticComputedValues.focalLengthMm,staticParams.integrationTimeSecs,0);
ptbNominalBackLMSIsomerizations = round(ptbNominalBackLMSIsomerizations);
ptbNominalLMSQuantalEfficiency = ptbNominalPhotorceptorsStruct.isomerizationAbsorbtance;
ptbNominalLMSEnergySensitivities = ptbNominalPhotorceptorsStruct.energyFundamentals;

[ptbAdjustedBackLMSIsomerizations,staticParams.pupilDiameterMm,ptbAdjustedPhotorceptorsStruct,ptbAdjustedIrradianceWattsPerM2] = ptbConeIsomerizationsFromSpectra(staticComputedValues.backSpd,staticComputedValues.wavelengthsNm,...
    staticParams.pupilDiameterMm,staticComputedValues.focalLengthMm,staticParams.integrationTimeSecs,theParams.macularPigmentDensityAdjust);
ptbAdjustedBackLMSIsomerizations = round(ptbAdjustedBackLMSIsomerizations);
ptbAdjustedLMSQuantalEfficiency = ptbAdjustedPhotorceptorsStruct.isomerizationAbsorbtance;
ptbAdjustedLMSEnergySensitivities = ptbAdjustedPhotorceptorsStruct.energyFundamentals;

%% Can simulate different types of color observers.  This is a little bit of a kluge, and
% may break comparisons between PTB and isetbio that are currently turned off anyway, for
% the dichromatic cases.
switch (theParams.OBSERVER_STATE)
    case 'LMandS'
        % Don't need to do anything here
    case 'LSonly'
        % Make the 'M' cone sensitivity the same as the L cone sensitivity.  This produces an
        % observer who has only L cones, although they remained labelled LMS.  This seems like
        % what we want if we accept the replacement hypothesis -- we'll use the same number of
        % cones in our decisions, but their spectral information will be limited to the L and S
        % cones.
        ptbAdjustedLMSQuantalEfficiency(2,:) = ptbAdjustedLMSQuantalEfficiency(1,:);
    case 'MSonly'
        ptbAdjustedLMSQuantalEfficiency(1,:) = ptbAdjustedLMSQuantalEfficiency(2,:);
    otherwise
        error('Unknown dichromat/trichromat type specified');
end

%%  Create a human cone mosaic sensor
%
% [**] BAW will fix up the field of view of this thing
% at some point.
scParams.sz = [128,192];
scParams.rgbDensities = staticParams.coneProportions;
scParams.coneAperture = staticParams.coneApertureMeters;
pixel = [];
cSensor = sensorCreate('human',pixel,scParams);
cSensor = sensorSet(cSensor,'exp time',staticParams.integrationTimeSecs);
cSensor = sensorSet(cSensor,'wave',staticComputedValues.wavelengthsNm);

% Put in PTB quantal sensitivities
tempWavelengths = sensorGet(cSensor,'wave');
if (any(find(tempWavelengths ~= staticComputedValues.wavelengthsNm)))
    error('Wavelength sampling not consistent throughout.');
end
cSensor = sensorSet(cSensor,'filter spectra',[zeros(size(ptbAdjustedLMSQuantalEfficiency',1),1) ptbAdjustedLMSQuantalEfficiency']);
sensorSetSizeToFOV(cSensor,0.9*staticParams.fieldOfViewDegrees);
%sensorConePlot(cSensor);

%% Set up cone conversions for computing stimuli.  These
% are done in a nominal LMS energy space, which may differ
% from the space of the simulated observer.
rgb2cones = ptbNominalLMSEnergySensitivities*staticComputedValues.displaySpd;
backLMS = rgb2cones*staticComputedValues.backRGB;

%% Create a test vector in a specified
% direction in cone space.
%
% To do this, we pick a test color direction
% and find the RGB values corresponding to it.
% Then we scale and add these to the background.
%
% Define test direction in cone excitation space
Lval = cos(theParams.cdAngle); Mval = sin(theParams.cdAngle);
testLMSOnUnitCircle = [Lval Mval 0]';

% Compute the RGB direction and scale so that it
% reaches to the edge of the gamut.
testRGBUnitCircle = inv(rgb2cones)*testLMSOnUnitCircle;
gamutScaleFactor = MaximizeGamutContrast(testRGBUnitCircle,staticComputedValues.backRGB);
testRGBGamut = gamutScaleFactor*testRGBUnitCircle;
testLMSGamut = rgb2cones*testRGBGamut;
testLMSContrast = testLMSGamut./backLMS;

% Scale again to put a bound on max contrast length.  This
% makes our sampling of the psychometric functions more
% efficient.
testContrastLength = norm(testLMSContrast);
if (testContrastLength > staticParams.testContrastLengthMax)
    testRGBGamut = (staticParams.testContrastLengthMax/testContrastLength)*testRGBGamut;
    testLMSGamut = (staticParams.testContrastLengthMax/testContrastLength)*testLMSGamut;
    testLMSContrast = testLMSGamut./backLMS;
    testContrastLength = norm(testLMSContrast);
end
results.testLMSContrast = testLMSContrast;
results.backgroundLMS = backLMS;
results.testLMSGamut = testLMSGamut;
        
% Plot comparison of iset and ptb irradiance, optionally
%
% PTB, conversion is pupilArea/(eyeLength^2).
% pi /(1 + 4*fN^2*(1+abs(m))^2)
%
% [**] This plot is currently broken because backOiD is no longer
% computed at this level.  
%if (runtimeParams.DO_SIM_PLOTS)
if (0)
    % Get background irradiance out of the optical image.
    %
    % [**] This currently works be using an ROI that was selected
    % by hand an then stored in a .mat file.  May want to
    % make this more programmatic.  But, we get this just
    temp = load('roiLocs');
    backUdata = plotOI(backOiD,'irradiance energy roi',temp.roiLocs);
    isetIrradianceWattsPerM2 = backUdata.y';
    
    % Make a new plot of PTB and iset irradiance.  Not quite
    % sure why we don't just add this to the window that
    % comes up in the call to PlotOI above.
    figure; hold on
    plot(staticComputedValues.wavelengthsNm,isetIrradianceWattsPerM2,'r');
    plot(staticComputedValues.wavelengthsNm,ptbAdjustedIrradianceWattsPerM2,'k');
    theRatio = isetIrradianceWattsPerM2 ./ ptbAdjustedIrradianceWattsPerM2;
end

%% Get noisy sensor values for both blank (background) and test intervals.
switch (staticParams.stimulus.type)
    case 'rgb_uniform'
        backStimulus = staticParams.stimulus;
        backStimulus.rgbVector = staticComputedValues.backRGB;
        blankConeResponses = getConeResponsesToStimRGB(cSensor,backStimulus,theParams,staticParams,runtimeParams,staticComputedValues);
        
        testStimulus = staticParams.stimulus;
        testStimulus.rgbVector = (staticComputedValues.backRGB + theParams.testLevel*testRGBGamut);
        testConeResponses = getConeResponsesToStimRGB(cSensor,testStimulus,theParams,staticParams,runtimeParams,staticComputedValues);
        
        % Plot out PTB and isetbio cone quantal spectral sensitivities, optionally
        if (runtimeParams.DO_SIM_PLOTS)
            figure; clf; hold on
            plot(SToWls(ptbAdjustedPhotorceptorsStruct.nomogram.S),ptbAdjustedPhotorceptorsStruct.isomerizationAbsorbtance(end:-1:1,:)');
            plot(blankConeResponses.isetLMSQuantalEfficiencyWavelengths,blankConeResponses.isetLMSQuantalEfficiencies(end:-1:1,:)',':');
            xlabel('Wavelength (nm)');
            ylabel('Isomerization Quantal Efficiency');
        end
        
        % Print out comparison of isomerization rates between PTB and iset, as well as PTB photoreceptor parameters.
        if (~runtimeParams.SIM_QUIET)
            fprintf('\tISET computes LMS isomerizations as: %d, %d, %d\n',blankConeResponses.isetLMSIsomerizations(1),blankConeResponses.isetLMSIsomerizations(2),blankConeResponses.isetLMSIsomerizations(3));
            fprintf('\tPTB computes LMS isomerizations as: %d, %d, %d\n',ptbAdjustedBackLMSIsomerizations(1),ptbAdjustedBackLMSIsomerizations(2),ptbAdjustedBackLMSIsomerizations(3));
            PrintPhotoreceptors(ptbAdjustedPhotorceptorsStruct);
        end

    otherwise
        error('Unknown stimulus type specified');
end

%% Compute second site responses
blankSecondSiteResponses = getSecondSiteResponsesFromConeResponses(blankConeResponses,theParams,staticParams,runtimeParams,staticComputedValues);
testSecondSiteResponses = getSecondSiteResponsesFromConeResponses(testConeResponses,theParams,staticParams,runtimeParams,staticComputedValues);

%% Get the data in a form to use with the classifier
classificationData = getTrainingAndValidationData(blankSecondSiteResponses,testSecondSiteResponses);
       
%% Plot the training and test data.  We'll plot the distribution of responses for one cone
if (runtimeParams.DO_SIM_PLOTS)
    % in each sensor class, with the distribution taken over our resampling of each
    % scene by the mosaic.
    %
    % This plot gives an idea of the separation between the classes, but doesn't
    % represent all the info that may be available to the classifier (for example,
    % if we add more receptors or do the TAFC calculation.
    if (~exist('f1','var'))
        f1 = vcNewGraphWin; hold on;
    else
        figure(f1); clf; hold on;
    end
    sym = {'b.','g.','r.','c.','k.'};
    az = 65.5; el = 30;
    index = find(classificationData.trainingLabels == classificationData.blankLabel);
    plot3(classificationData.trainingData(index,1), ...
        classificationData.trainingData(index,2), ...
        classificationData.trainingData(index,3),'ko','MarkerFaceColor','k');
    index = find(classificationData.trainingLabels == classificationData.testLabel);
    plot3(classificationData.trainingData(index,1), ...
        classificationData.trainingData(index,2), ...
        classificationData.trainingData(index,3),'ro','MarkerFaceColor','r');
    
    index = find(classificationData.validateLabels == classificationData.blankLabel);
    plot3(classificationData.validateData(index,1), ...
        classificationData.validateData(index,2), ...
        classificationData.validateData(index,3),'ko');
    index = find(classificationData.validateLabels == classificationData.testLabel);
    plot3(classificationData.validateData(index,1), ...
        classificationData.validateData(index,2), ...
        classificationData.validateData(index,3),'ro');
    xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on
    drawnow;
end

%% Build and test the SVM model
%
% Performance is very sensitive to the options.
%   The default option of a radial basis kernal '-t 2')
%   leads to good performance on the training set but does
%   not generalize well to the test set, at least
%   with the default choice of fgamma.
%
%   Using a linear kernal works quite well on test and training,
%   in initial tests.
%
% Convert training and test data to TAFC form, if we want.
if (theParams.DO_TAFC_CLASSIFIER)
    classificationData = oneIntervalToTwoIntervalClassificationData(classificationData);
end
       
%% Train and evaluate SVM
%
% Train on one part of the data, evaluate on the other.  
svmOpts = '-s 0 -t 0';
predictOpts = '';
if (runtimeParams.SIM_QUIET)
    svmOpts =  [svmOpts ' -q'];
    predictOpts = [predictOpts ' -q'];
end
svmModel = svmtrain(classificationData.trainingLabels, classificationData.trainingData, svmOpts);
[svmTrainingPredictedLabels] = svmpredict(classificationData.trainingLabels, classificationData.trainingData, svmModel, predictOpts);
[svmValidatePredictedLabels] = svmpredict(classificationData.validateLabels, classificationData.validateData, svmModel, predictOpts);
trainingFractionCorrect = length(find(svmTrainingPredictedLabels == classificationData.trainingLabels))/length(classificationData.trainingLabels);
validateFractionCorrect = length(find(svmValidatePredictedLabels == classificationData.validateLabels))/length(classificationData.validateLabels);

%% Store the results for return
results.nCorrectResponses = length(find(svmValidatePredictedLabels == classificationData.validateLabels));
results.nTotalResponses = length(classificationData.validateLabels);
results.fractionCorrect = validateFractionCorrect;
if (~runtimeParams.SIM_QUIET)
    fprintf('\tClassifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
end

% Plot the training and test data.  We'll plot the distribution of responses for one cone
% in each sensor class, with the distribution taken over our resampling of each
% scene by the mosaic.
%
% Some thought is required about how to make a useful plot for the TAFC case, skipping
% it for now.
if (runtimeParams.DO_SIM_PLOTS)
    if (~theParams.DO_TAFC_CLASSIFIER)
        % Indices for coloring the plot
        indexTG = find(svmTrainingPredictedLabels == classificationData.trainingLabels);
        indexTR = find(svmTrainingPredictedLabels ~= classificationData.trainingLabels);
        indexVG = find(svmValidatePredictedLabels == classificationData.validateLabels);
        indexVR = find(svmValidatePredictedLabels ~= classificationData.validateLabels);

        if (~exist('f2','var'))
            f2 = vcNewGraphWin; hold on;
        else
            figure(f2); clf; hold on;
        end
        sym = {'b.','g.','r.','c.','k.'};
        az = 65.5; el = 30;
        plot3(classificationData.trainingData(indexTG,1), ...
            classificationData.trainingData(indexTG,2), ...
            classificationData.trainingData(indexTG,3),'go','MarkerFaceColor','g');
        plot3(classificationData.trainingData(indexTR,1), ...
            classificationData.trainingData(indexTR,2), ...
            classificationData.trainingData(indexTR,3),'ro','MarkerFaceColor','r');
        
        plot3(classificationData.validateData(indexVG,1), ...
            classificationData.validateData(indexVG,2), ...
            classificationData.validateData(indexVG,3),'go');
        plot3(classificationData.validateData(indexVR,1), ...
            classificationData.validateData(indexVR,2), ...
            classificationData.validateData(indexVR,3),'ro');
        xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on
        drawnow;
    end
end