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
if (runtimeParams.DO_SIM_PLOTS)
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


%% Get cell arrays of noisy sensor values for both blank (background) and test intervals.
[blankVectors,blankOneConeEachClassStartIndices,isetBlankLMSIsomerizations,isetLMSQuantalEfficiencies,isetLMSQuantalEfficiencyWavelengths] = ...
    getConeResponsesToStimRGB(cSensor,staticComputedValues.backRGB,theParams,staticParams,runtimeParams,staticComputedValues);
[testVectors,testOneConeEachClassStartIndices,isetTestLMSIsomerizations] = getConeResponsesToStimRGB(cSensor,(staticComputedValues.backRGB + theParams.testLevel*testRGBGamut),theParams,staticParams,runtimeParams,staticComputedValues);

%% Plot out PTB and isetbio cone quantal spectral sensitivities, optionally
if (runtimeParams.DO_SIM_PLOTS)
    figure; clf; hold on
    plot(SToWls(ptbAdjustedPhotorceptorsStruct.nomogram.S),ptbAdjustedPhotorceptorsStruct.isomerizationAbsorbtance(end:-1:1,:)');
    plot(isetLMSQuantalEfficiencyWavelengths,isetLMSQuantalEfficiencies(end:-1:1,:)',':');
    xlabel('Wavelength (nm)');
    ylabel('Isomerization Quantal Efficiency');
end

% Compute baseline Poisson sd from background
%backPoissonSd = sqrt(mean(isetBackLMSIsomerizations));

%% Print out comparison of isomerization rates between PTB and iset, as well as PTB photoreceptor parameters.
if (~runtimeParams.SIM_QUIET)
    fprintf('\tISET computes LMS isomerizations as: %d, %d, %d\n',isetBlankLMSIsomerizations(1),isetBlankLMSIsomerizations(2),isetBlankLMSIsomerizations(3));
    fprintf('\tPTB computes LMS isomerizations as: %d, %d, %d\n',ptbAdjustedBackLMSIsomerizations(1),ptbAdjustedBackLMSIsomerizations(2),ptbAdjustedBackLMSIsomerizations(3));
    PrintPhotoreceptors(ptbAdjustedPhotorceptorsStruct);
end

%% Build a clasifier on the training set
%
% Set up training/test data and label them.
%
% The matrix backVectors has one column for each sample,
% with the entries going down the rows representing
% different cones.
%
% Same for the matrix testVectors.  
%
% Currently just using the responses of just one
% cone from each class, as a place to start, so
% for each column we pull out one entry for cones
% of each class.
%
% There is also a transpose, so that each data
% instance is in a single row, rather than in a
% single column.  This matches what the classification
% routines want.
blankLabel = -1;
testLabel = 1;
blankLMSFullSet = [blankVectors(blankOneConeEachClassStartIndices(1),:) ; ...
    blankVectors(blankOneConeEachClassStartIndices(2),:) ; ...
    blankVectors(blankOneConeEachClassStartIndices(3),:)]';
testLMSFullSet = [testVectors(testOneConeEachClassStartIndices(1),:) ; ...
    testVectors(testOneConeEachClassStartIndices(2),:) ; ...
    testVectors(testOneConeEachClassStartIndices(3),:)]';

% Implement a quick and dirty surround for L and M
% cone. 
switch theParams.surroundType
    case 'none'
        % No surround. Do nothing here.
    case 'rdraw'
        % This subtracts a random mix of L and M cone
        % responses from the single L and M cone response.
        % This means that more cones are being used in the
        % classification than in the case without a surround,
        % which may not be a good idea.
        if (theParams.surroundSize > 0 & theParams.surroundWeight > 0)
            for k = 1:staticParams.nDrawsPerTestStimulus
                % Blank L cone
                LCenterVal = blankLMSFullSet(k,1);
                otherConesIndex = Shuffle(1:blankOneConeEachClassStartIndices(3)-1);
                MSurroundVals = blankVectors(otherConesIndex(1:theParams.surroundSize),k);
                LOpponentVal = LCenterVal - theParams.surroundWeight*sum(MSurroundVals)/theParams.surroundSize;
                blankLMSFullset(1,k) = LOpponentVal;
                
                % Blank M cone
                MCenterVal = blankLMSFullSet(k,2);
                otherConesIndex = Shuffle(1:blankOneConeEachClassStartIndices(3)-1);
                LSurroundVals = blankVectors(otherConesIndex(1:theParams.surroundSize),k);
                MOpponentVal = MCenterVal - theParams.surroundWeight*sum(LSurroundVals)/theParams.surroundSize;
                blankLMSFullset(2,k) = MOpponentVal;
                
                % Test L cone
                LCenterVal = testLMSFullSet(k,1);
                otherConesIndex = Shuffle(1:testOneConeEachClassStartIndices(3)-1);
                MSurroundVals = testVectors(otherConesIndex(1:theParams.surroundSize),k);
                LOpponentVal = LCenterVal - theParams.surroundWeight*sum(MSurroundVals)/theParams.surroundSize;
                testLMSFullset(1,k) = LOpponentVal;
                
                % Test M cone
                MCenterVal = testLMSFullSet(k,2);
                otherConesIndex = Shuffle(1:testOneConeEachClassStartIndices(3)-1);
                LSurroundVals = testVectors(otherConesIndex(1:theParams.surroundSize),k);
                MOpponentVal = MCenterVal - theParams.surroundWeight*sum(LSurroundVals)/theParams.surroundSize;
                testLMSFullset(2,k) = MOpponentVal;
            end      
        end
    case 'determ'
        % This subtracts the M from the L and the L from the M.
        %
        % This version I think should not have an effect unless we add some noise
        % after the recombination.
        if (theParams.surroundWeight > 0)
            for k = 1:staticParams.nDrawsPerTestStimulus
                % Blank
                LCenterVal = blankLMSFullSet(k,1);
                MCenterVal = blankLMSFullSet(k,2);
                blankLMSFullSet(k,1) = LCenterVal - theParams.surroundWeight*MCenterVal;
                blankLMSFullSet(k,2) = MCenterVal - theParams.surroundWeight*LCenterVal;
                
                % Test 
                LCenterVal = testLMSFullSet(k,1);
                MCenterVal = testLMSFullSet(k,2);
                testLMSFullSet(k,1) = LCenterVal - theParams.surroundWeight*MCenterVal;
                testLMSFullSet(k,2) = MCenterVal - theParams.surroundWeight*LCenterVal;
            end
        end
        
    otherwise
        error('Unknown surround type specified');     
end

%% Add Poisson noise at a second site
%
% You can scale the Poisson variance by the parameter
% params.opponetLevelNoiseSd.  Setting this to 1 gives
% Poisson.
if (theParams.opponentLevelNoiseSd > 0)
    blankVarMat = abs(blankLMSFullSet);
    blankSdMat = sqrt(blankVarMat);
    blankLMSFullSet = blankLMSFullSet + theParams.opponentLevelNoiseSd*blankSdMat.*randn(size(blankLMSFullSet));
    
    testVarMat = abs(testLMSFullSet);
    testSdMat = sqrt(testVarMat);
    testLMSFullSet = testLMSFullSet + theParams.opponentLevelNoiseSd*testSdMat.*randn(size(testLMSFullSet));
end

%% Pull apart into training and validation datasets.
fullData = [blankLMSFullSet ; testLMSFullSet];
fullLabels = [blankLabel*ones(size(blankLMSFullSet,1),1) ; testLabel*ones(size(testLMSFullSet,1),1)];
fullDataN = size(fullData,1);
trainingDataN = round(fullDataN/2);
testDataN = fullDataN-trainingDataN;
indices = Shuffle(1:size(fullData,1));
trainingData = fullData(indices(1:trainingDataN),:);
trainingLabels = fullLabels(indices(1:trainingDataN));
validateData = fullData(indices(trainingDataN+1:end),:);
validateLabels = fullLabels(indices(trainingDataN+1:end));

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
    index = find(trainingLabels == blankLabel);
    plot3(trainingData(index,1), ...
        trainingData(index,2), ...
        trainingData(index,3),'ko','MarkerFaceColor','k');
    index = find(trainingLabels == testLabel);
    plot3(trainingData(index,1), ...
        trainingData(index,2), ...
        trainingData(index,3),'ro','MarkerFaceColor','r');
    
    index = find(validateLabels == blankLabel);
    plot3(validateData(index,1), ...
        validateData(index,2), ...
        validateData(index,3),'ko');
    index = find(validateLabels == testLabel);
    plot3(validateData(index,1), ...
        validateData(index,2), ...
        validateData(index,3),'ro');
    xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on
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
    % Build up a TAFC training and test set from the one interval data
    oneIntervalDataDimension = size(trainingData,2);
    nTrainingData = 2*size(trainingData,1);
    trainingBlankIndices = find(trainingLabels == blankLabel);
    trainingTestIndices = find(trainingLabels == testLabel);
    validateBlankIndices = find(validateLabels == blankLabel);
    validateTestIndices = find(validateLabels == testLabel);
    tafcTrainingData = zeros(nTrainingData,2*oneIntervalDataDimension);
    tafcTrainingLabels = zeros(nTrainingData,1);
    tafcValidateData = zeros(nTrainingData,2*oneIntervalDataDimension);
    tafcValidateLabels = zeros(nTrainingData,1);
    for tt = 1:nTrainingData
        % Training set
        %
        % Flip a coin to decide whether test is in first or second interval
        if (CoinFlip(1,0.5))
            temp1 = Shuffle(trainingTestIndices);
            temp2 = Shuffle(trainingBlankIndices);
            tafcTrainingLabels(tt) = blankLabel;
        else
            temp1 = Shuffle(trainingBlankIndices);
            temp2 = Shuffle(trainingTestIndices);
            tafcTrainingLabels(tt) = testLabel;
        end
        tafcTrainingData(tt,1:oneIntervalDataDimension) = trainingData(temp1(1),:);
        tafcTrainingData(tt,oneIntervalDataDimension+1:2*oneIntervalDataDimension) = trainingData(temp2(1),:);
        
        % Validation set
        %
        % Same logic
        if (CoinFlip(1,0.5))
            temp1 = Shuffle(validateTestIndices);
            temp2 = Shuffle(validateBlankIndices);
            tafcValidateLabels(tt) = blankLabel;
        else
            temp1 = Shuffle(validateBlankIndices);
            temp2 = Shuffle(validateTestIndices);
            tafcValidateLabels(tt) = testLabel;
        end
        tafcValidateData(tt,1:oneIntervalDataDimension) = validateData(temp1(1),:);
        tafcValidateData(tt,oneIntervalDataDimension+1:2*oneIntervalDataDimension) = validateData(temp2(1),:);
    end
    
    svmOpts = '-s 0 -t 0';
    predictOpts = '';
    if (runtimeParams.SIM_QUIET)
        svmOpts =  [svmOpts ' -q'];
        predictOpts = [predictOpts ' -q'];
    end
    svmModel = svmtrain(tafcTrainingLabels, tafcTrainingData, svmOpts);
    [svmTrainingPredictedLabels] = svmpredict(tafcTrainingLabels, tafcTrainingData, svmModel, predictOpts);
    [svmValidatePredictedLabels] = svmpredict(tafcValidateLabels, tafcValidateData, svmModel, predictOpts);
    trainingFractionCorrect = length(find(svmTrainingPredictedLabels == tafcTrainingLabels))/length(tafcTrainingLabels);
    validateFractionCorrect = length(find(svmValidatePredictedLabels == tafcValidateLabels))/length(tafcValidateLabels);
    results.nCorrectResponses = length(find(svmValidatePredictedLabels == tafcValidateLabels));
    results.nTotalResponses = length(tafcValidateLabels);
    results.fractionCorrect = validateFractionCorrect;
    if (~runtimeParams.SIM_QUIET)
        fprintf('\tClassifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
    end
    
    % Indices for plots below
    indexTG = find(svmTrainingPredictedLabels == tafcTrainingLabels);
    indexTR = find(svmTrainingPredictedLabels ~= tafcTrainingLabels);
    indexVG = find(svmValidatePredictedLabels == tafcValidateLabels);
    indexVR = find(svmValidatePredictedLabels ~= tafcValidateLabels);
else
    % Just do the one interval analysis
    svmOpts = '-s 0 -t 0';
    predictOpts = '';
    if (runtimeParams.SIM_QUIET)
        svmOpts =  [svmOpts ' -q'];
        predictOpts = [predictOpts ' -q'];
    end
    svmModel = svmtrain(trainingLabels, trainingData, svmOpts);
    [svmTrainingPredictedLabels] = svmpredict(trainingLabels, trainingData, svmModel, predictOpts);
    [svmValidatePredictedLabels] = svmpredict(validateLabels, validateData, svmModel, predictOpts);
    trainingFractionCorrect = length(find(svmTrainingPredictedLabels == trainingLabels))/length(trainingLabels);
    validateFractionCorrect = length(find(svmValidatePredictedLabels == validateLabels))/length(validateLabels);
    results.nCorrectResponses = length(find(svmValidatePredictedLabels == validateLabels));
    results.nTotalResponses = length(validateLabels);
    results.fractionCorrect = validateFractionCorrect;
    if (~runtimeParams.SIM_QUIET)
        fprintf('\tClassifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
    end
    
    % Indices for plots below
    indexTG = find(svmTrainingPredictedLabels == trainingLabels);
    indexTR = find(svmTrainingPredictedLabels ~= trainingLabels);
    indexVG = find(svmValidatePredictedLabels == validateLabels);
    indexVR = find(svmValidatePredictedLabels ~= validateLabels);
end

% Plot the training and test data.  We'll plot the distribution of responses for one cone
% in each sensor class, with the distribution taken over our resampling of each
% scene by the mosaic.
%
% Some thought is required about how to make a useful plot for the TAFC case, skipping
% it for now.
if (runtimeParams.DO_SIM_PLOTS)
    if (~theParams.theParams.DO_TAFC_CLASSIFIER)
        if (~exist('f2','var'))
            f2 = vcNewGraphWin; hold on;
        else
            figure(f2); clf; hold on;
        end
        sym = {'b.','g.','r.','c.','k.'};
        az = 65.5; el = 30;
        plot3(trainingData(indexTG,1), ...
            trainingData(indexTG,2), ...
            trainingData(indexTG,3),'go','MarkerFaceColor','g');
        plot3(trainingData(indexTR,1), ...
            trainingData(indexTR,2), ...
            trainingData(indexTR,3),'ro','MarkerFaceColor','r');
        
        plot3(validateData(indexVG,1), ...
            validateData(indexVG,2), ...
            validateData(indexVG,3),'go');
        plot3(validateData(indexVR,1), ...
            validateData(indexVR,2), ...
            validateData(indexVR,3),'ro');
        xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on
    end
end