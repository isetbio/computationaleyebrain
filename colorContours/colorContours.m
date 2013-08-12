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
%  * Loop over different macular pigment density adjustments.
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
%  * Add option for dichromats.
%  * Add more sensible code to control spatial integration area.
%  * Add eye movements.
%  * Add simple foveal midget ganglion cell model.
%  * Break this big long script into sensible subfunctions.
%  * Smarter choice of test levels to improve efficiency.
%
% Known bugs:
%  * 8/12/13  For some reason, the contour plot doesn't get written with
%             big simulations.  
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
% 8/2/13  DHB/BW Our excellent adventure commences.
% 8/3/13  DHB/BW Tune up and add classifier.
% 8/4/13  DHB/BW Check irradiance calcs between ISET and PTB.
% 8/6/13  DHB/BW Got contour going.  End of this DHB visit west.
% 8/11/13 DHB    Added ellipse fitting, macular pigment adjustment
%         DHB    Added TAFC option.  Not fully tested.

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

nColorDirections = 20;                          % Number of color directions for contour.
dirAngleMax = 2*pi;                             % Use pi for sampling directions from hemicircle, 2*pi for whole circle

nTestLevels = 15;                               % Number of test levels to simulate for each test direction psychometric function.
nDrawsPerTestStimulus = 500;                    % Number of noise draws used in the simulations, per test stimulus
noiseType = 1;                                  % Noise type passed to isetbio routines.  1 -> Poisson.

macularPigmentDensityAdjustments = [-0.25 -0.10 0 0.1 0.25]; % Amount to adjust macular pigment density for cone fundamentals of simulated observer.
                                                % Note that stimuli are computed for a nominal (no adjustment) observer.

criterionCorrect = 0.82;                        % Fraction correct for definition of threshold in TAFC simulations.

DO_TAFC_CLASSIFIER = true;                      % Build the classifier assuming a TAFC design.

QUICK_TEST_PARAMS = true;                       % Set to true to override parameters with a small number of trials for debugging.

%% Process quick test option
if (QUICK_TEST_PARAMS)
    nColorDirections = 4;
    dirAngleMax = pi;                           % Use pi for sampling directions from hemicircle, 2*pi for whole circle
    nTestLevels = 4;
    nDrawsPerTestStimulus = 100;
    macularPigmentDensityAdjustments = [-0.3 0 0.3];
    DO_TAFC_CLASSIFIER = false;
end
    
%% Make sure random number generator seed is different each run.
ClockRandSeed;

%% Get stats path off of the path, temporarily.
%
% This avoids a name space conflict with the libsvn
% routines that we use here.
if (exist('RemoveMatchingPaths','file'))
    path(RemoveMatchingPaths(path,'stats'));
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
% oiWindow;
% vcNewGraphWin; plotOI(oiD,'psf')

%% Loop over macular pigment adjustments
%
% Output saved at end of each loop iteration
nMacularPigmentDensitiesAdjustments = length(macularPigmentDensityAdjustments);
for m = 1:nMacularPigmentDensitiesAdjustments
    macularPigmentDensityAdjust = macularPigmentDensityAdjustments(m);
    
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
    [ptbNominalBackLMSIsomerizations,pupilDiameterMm,ptbNominalPhotorceptorsStruct,ptbNominalIrradianceWattsPerM2] = ptbConeIsomerizationsFromSpectra(backSpd,wavelengthsNm,...
        pupilDiameterMm,focalLengthMm,integrationTimeSecs,0);
    ptbNominalBackLMSIsomerizations = round(ptbNominalBackLMSIsomerizations);
    ptbNominalLMSQuantalEfficiency = ptbNominalPhotorceptorsStruct.isomerizationAbsorbtance;
    ptbNominalLMSEnergySensitivities = ptbNominalPhotorceptorsStruct.energyFundamentals;
    
    [ptbAdjustedBackLMSIsomerizations,pupilDiameterMm,ptbAdjustedPhotorceptorsStruct,ptbAdjustedIrradianceWattsPerM2] = ptbConeIsomerizationsFromSpectra(backSpd,wavelengthsNm,...
        pupilDiameterMm,focalLengthMm,integrationTimeSecs,macularPigmentDensityAdjust);
    ptbAdjustedBackLMSIsomerizations = round(ptbAdjustedBackLMSIsomerizations);
    ptbAdjustedLMSQuantalEfficiency = ptbAdjustedPhotorceptorsStruct.isomerizationAbsorbtance;
    ptbAdjustedLMSEnergySensitivities = ptbAdjustedPhotorceptorsStruct.energyFundamentals;
    
    %%  Create a human cone mosaic sensor
    %
    % [**] BAW will fix up the field of view of this thing
    % at some point.
    params.sz = [128,192];
    params.rgbDensities = coneProportions;
    params.coneAperture = coneApertureMeters;
    pixel = [];
    cSensor = sensorCreate('human',pixel,params);
    cSensor = sensorSet(cSensor,'exp time',integrationTimeSecs);
    cSensor = sensorSet(cSensor,'wave',wavelengthsNm);
    
    % Put in PTB quantal sensitivities
    isetLMSQuantalEfficiencyWavelengths = sensorGet(cSensor,'wave');
    if (any(find(isetLMSQuantalEfficiencyWavelengths ~= wavelengthsNm)))
        error('Wavelength sampling not consistent throughout.');
    end
    cSensor = sensorSet(cSensor,'filter spectra',[zeros(size(ptbAdjustedLMSQuantalEfficiency',1),1) ptbAdjustedLMSQuantalEfficiency']);
    sensorSetSizeToFOV(cSensor,0.9*fieldOfViewDegrees);
    sensorFieldOfView = sensorGet(cSensor,'fov',sceneB,oiD);
    %sensorConePlot(cSensor);
    
    %% Set up cone conversions
    rgb2cones = ptbNominalLMSEnergySensitivities*displaySpd;
    backLMS = rgb2cones*backRGB;
    
    %% Create a test vector in a specified
    % direction in cone space.
    %
    % To do this, we pick a test color direction
    % and find the RGB values corresponding to it.
    % Then we scale and add these to the background.
    %
    % Define test direction in cone excitation space
    cdAngles = linspace(0,dirAngleMax,nColorDirections+1);
    cdAngles = cdAngles(1:end-1);
    for cd = 1:nColorDirections
        Lval = cos(cdAngles(cd)); Mval = sin(cdAngles(cd));
        testLMSUnitCircle = [Lval Mval 0]';
        
        % Compute the RGB direction and scale so that it
        % reaches to the edge of the gamut.
        testRGBUnitCircle = inv(rgb2cones)*testLMSUnitCircle;
        gamutScaleFactor = MaximizeGamutContrast(testRGBUnitCircle,backRGB);
        testRGBGamut = gamutScaleFactor*testRGBUnitCircle;
        testLMSGamut = rgb2cones*testRGBGamut;
        testLMSContrast = testLMSGamut./backLMS;
        
        %% Pass the background through the optics
        backOiD = oiCompute(oiD,sceneB);
        vcAddAndSelectObject(backOiD);
        %oiWindow;
        
        % Plot comparison of iset and ptb irradiance, optionally
        %
        % PTB, conversion is pupilArea/(eyeLength^2).
        % pi /(1 + 4*fN^2*(1+abs(m))^2)
        PLOT_COMPARE_IRRADIANCE = 0;
        if (PLOT_COMPARE_IRRADIANCE)
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
            plot(wavelengthsNm,isetIrradianceWattsPerM2,'r');
            plot(wavelengthsNm,ptbAdjustedIrradianceWattsPerM2,'k');
            theRatio = isetIrradianceWattsPerM2 ./ ptbAdjustedIrradianceWattsPerM2;
        end
        
        %% Compute noise free background sensor image
        backCSensorNF = sensorComputeNoiseFree(cSensor,backOiD);
        vcAddAndSelectObject(backCSensorNF);
        % [**] Something broke this call to sensorWindow
        % Might be related to our fov mucking.  Or maybe
        % it would work again if we just uncommented it.
        % sensorWindow;
        
        %% Get noise free cone isomerizations for background
        %
        % Compare with what PTB routines compute for the same stimuli.
        %
        % We pick the max to avoid the edge effects in the sensor image.
        for ii = 1:nSensorClasses
            backSensorValsNF{ii} = sensorGet(backCSensorNF,'electrons',isetSensorConeSlots(ii));
            isetBackLMSIsomerizations(ii) = round(max(backSensorValsNF{ii}));
        end
        
        % Print out the comparison as well as PTB parameters.
        fprintf('ISET computes LMS isomerizations as: %d, %d, %d\n',isetBackLMSIsomerizations(1),isetBackLMSIsomerizations(2),isetBackLMSIsomerizations(3));
        fprintf('PTB computes LMS isomerizations as: %d, %d, %d\n',ptbAdjustedBackLMSIsomerizations(1),ptbAdjustedBackLMSIsomerizations(2),ptbAdjustedBackLMSIsomerizations(3));
        PrintPhotoreceptors(ptbAdjustedPhotorceptorsStruct);
        
        % Get iset LMS quantal efficiences
        temp = sensorGet(backCSensorNF,'spectralqe')';
        isetLMSQuantalEfficiencyWavelengths = sensorGet(backCSensorNF,'wave');
        isetLMSQuantalEfficiences = temp(isetSensorConeSlots,:);
        
        % Plot out PTB and isetbio cone quantal spectral sensitivities, optionally
        PLOT_COMPARE_CONEQE = 0;
        if (PLOT_COMPARE_CONEQE)
            figure; clf; hold on
            plot(SToWls(ptbAdjustedPhotorceptorsStruct.nomogram.S),ptbAdjustedPhotorceptorsStruct.isomerizationAbsorbtance(end:-1:1,:)');
            plot(isetLMSQuantalEfficiencyWavelengths,isetLMSQuantalEfficiences(end:-1:1,:)',':');
            xlabel('Wavelength (nm)');
            ylabel('Isomerization Quantal Efficiency');
        end
        
        %% Loop over a set of test levels and get classifier
        % performance for each.
        testLevels = linspace(0,1,nTestLevels);
        for t = 1:length(testLevels)
            % Set test level
            fprintf('Calculations for test level %d of %d\n',t,length(testLevels));
            testLevel = testLevels(t);
            testRGBForThisLevel = (backRGB + testLevel*testRGBGamut).^(1/gammaValue);
            
            % Make test scene, in same fashion as we made the
            % background scene.
            tmp = ones(scenePixels,scenePixels,3);
            [tmp,row,col] = RGB2XWFormat(tmp);
            tmp = tmp*diag(testRGBForThisLevel(:));
            testImg = XW2RGBFormat(tmp,row,col);
            imwrite(testImg,'testFile.png','png');
            sceneT = sceneFromFile('testFile.png','rgb',[],'LCD-Apple.mat',wavelengthsNm);
            sceneT = sceneSet(sceneT,'name','test');
            sceneT = sceneSet(sceneT,'fov',2);
            vcAddAndSelectObject(sceneT);
            %sceneWindow;
            
            %% Pass test image through the optics
            testOiD = oiCompute(oiD,sceneT);
            
            %% Get multivariate sample distribution of LMS
            % responses out of the sensor objects.
            %
            % Each time through the loop we do a new instantiation
            % of the poisson noise.
            nSensorClasses = length(isetSensorConeSlots);
            testCSensorNF = sensorComputeNoiseFree(cSensor,testOiD);
            backVoltImage = sensorComputeSamples(backCSensorNF,nDrawsPerTestStimulus,noiseType);
            testVoltImage = sensorComputeSamples(testCSensorNF,nDrawsPerTestStimulus,noiseType);
            for k = 1:nDrawsPerTestStimulus
                for ii = 1:nSensorClasses
                    backCSensorTemp = sensorSet(backCSensorNF,'volts',backVoltImage(:,:,k));
                    backSensorVals{k,ii} = sensorGet(backCSensorTemp,'electrons',isetSensorConeSlots(ii));
                    testCSensorTemp = sensorSet(testCSensorNF,'volts',testVoltImage(:,:,k));
                    testSensorVals{k,ii} = sensorGet(testCSensorTemp,'electrons',isetSensorConeSlots(ii));
                end
            end
            
            %% Pop last instantions into windows for viewing etc.
            backCSensor = sensorSet(backCSensorNF,'name','Background');
            %vcAddAndSelectObject(backCSensorNF); sensorWindow;
            
            testCSensor = sensorSet(testCSensorNF,'name','Test');
            %vcAddAndSelectObject(testCSensorNF); sensorWindow;
            
            %% We want to control the integration area.
            %
            % Here we'll do this by specifying the fraction of
            % the total mosaic to use.
            %
            % First figure out length of sample vector
            fractionUse = 0.001;
            for ii = 1:nSensorClasses
                nUse(ii) = round(fractionUse*length(backSensorVals{1,ii}));
            end
            nUseAll = sum(nUse);
            
            % Now draw samples
            backVec = zeros(nUseAll,nDrawsPerTestStimulus);
            testVec = zeros(nUseAll,nDrawsPerTestStimulus);
            typeVec = zeros(nUseAll,1);
            oneConeEachClassIndices = zeros(nSensorClasses,1);
            for k = 1:nDrawsPerTestStimulus
                if (rem(k,10) == 0)
                    fprintf('\tGetting cone catches for draw %d of %d\n',k,nDrawsPerTestStimulus);
                end
                
                startIndex = 1;
                for ii = 1:nSensorClasses
                    % Pull out a set of randomly chosen responses for this sensor class
                    % and tuck it into the response vector, for both background and test.
                    % Also store the sensor class index for each stored response.
                    endIndex = startIndex + nUse(ii) - 1;
                    
                    % Currently, this shufles the responses of each class before
                    % pulling out the requisite number of responses, for some reason.
                    % Better will be to fix this up to be more realistic.
                    temp = Shuffle(backSensorVals{k,ii});
                    backVec(startIndex:endIndex,k) = temp(1:nUse(ii));
                    temp = Shuffle(testSensorVals{k,ii});
                    testVec(startIndex:endIndex,k) = temp(1:nUse(ii));
                    if (k == 1)
                        typeVec(startIndex:endIndex) = ii;
                        oneConeEachClassIndices(ii) = startIndex;
                    end
                    
                    startIndex = endIndex+1;
                end
            end
            
            %% Build a clasifier on the training set
            %
            % Set up training/test data and label them.  Here
            % we use the reponses of just one cone from each class,
            % just as a place to start.
            blankLabel = -1;
            testLabel = 1;
            backgroundLMSTraining = [backVec(oneConeEachClassIndices(1),:) ; ...
                backVec(oneConeEachClassIndices(2),:) ; ...
                backVec(oneConeEachClassIndices(3),:)]';
            testLMSTraining = [testVec(oneConeEachClassIndices(1),:) ; ...
                testVec(oneConeEachClassIndices(2),:) ; ...
                testVec(oneConeEachClassIndices(3),:)]';
            fullData = [backgroundLMSTraining ; testLMSTraining];
            fullLabels = [blankLabel*ones(size(backgroundLMSTraining,1),1) ; testLabel*ones(size(testLMSTraining,1),1)];
            fullDataN = size(fullData,1);
            trainingDataN = round(fullDataN/2);
            testDataN = fullDataN-trainingDataN;
            indices = Shuffle(1:size(fullData,1));
            trainingData = fullData(indices(1:trainingDataN),:);
            trainingLabels = fullLabels(indices(1:trainingDataN));
            validateData = fullData(indices(trainingDataN+1:end),:);
            validateLabels = fullLabels(indices(trainingDataN+1:end));
            
            %% Plot the training and test data.  We'll plot the distribution of responses for one cone
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
            if (DO_TAFC_CLASSIFIER)
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
                svmModel = svmtrain(tafcTrainingLabels, tafcTrainingData, svmOpts);
                [svmTrainingPredictedLabels] = svmpredict(tafcTrainingLabels, tafcTrainingData, svmModel);
                [svmValidatePredictedLabels] = svmpredict(tafcValidateLabels, tafcValidateData, svmModel);
                trainingFractionCorrect = length(find(svmTrainingPredictedLabels == tafcTrainingLabels))/length(tafcTrainingLabels);
                validateFractionCorrect = length(find(svmValidatePredictedLabels == tafcValidateLabels))/length(tafcValidateLabels);
                nCorrectResponses(t) = length(find(svmValidatePredictedLabels == tafcValidateLabels));
                nTotalResponses(t) = length(tafcValidateLabels);
                fractionCorrect(t) = validateFractionCorrect;
                fprintf('Classifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
                
                % Indices for plots below
                indexTG = find(svmTrainingPredictedLabels == tafcTrainingLabels);
                indexTR = find(svmTrainingPredictedLabels ~= tafcTrainingLabels);
                indexVG = find(svmValidatePredictedLabels == tafcValidateLabels);
                indexVR = find(svmValidatePredictedLabels ~= tafcValidateLabels);
            else
                % Just do the one interval analysis
                svmOpts = '-s 0 -t 0';
                svmModel = svmtrain(trainingLabels, trainingData, svmOpts);
                [svmTrainingPredictedLabels] = svmpredict(trainingLabels, trainingData, svmModel);
                [svmValidatePredictedLabels] = svmpredict(validateLabels, validateData, svmModel);
                trainingFractionCorrect = length(find(svmTrainingPredictedLabels == trainingLabels))/length(trainingLabels);
                validateFractionCorrect = length(find(svmValidatePredictedLabels == validateLabels))/length(validateLabels);
                nCorrectResponses(t) = length(find(svmValidatePredictedLabels == validateLabels));
                nTotalResponses(t) = length(validateLabels);
                fractionCorrect(t) = validateFractionCorrect;
                fprintf('Classifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
                
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
            if (~DO_TAFC_CLASSIFIER)
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
        ylim([0.5 1]);
        drawnow;
        
        % Print threshold
        fprintf('%d%% correct threshold is %0.1f\n',round(100*criterionCorrect),thresholdEst);
        
        % Store results for this color direction
        results(cd).thresholdLevel = thresholdEst;
        results(cd).testLMSContrast = testLMSContrast;
        results(cd).thresholdLMSContrast = thresholdEst*testLMSContrast;
        results(cd).backgroundLMS = backLMS;
        results(cd).testLMSGamut = testLMSGamut;
    end
    
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
        [ellipseZ, ellipseA, ellipseB, ellipseAlpha] = fitellipse([LContourPoints' ; MContourPoints']);
        plotellipse(ellipseZ,ellipseA,ellipseB,ellipseAlpha,'r');
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
        outName = sprintf('colorContour_TAFC_%d',round(100*macularPigmentDensityAdjust));
    else
        outName = sprintf('colorContour_YN_%d',round(100*macularPigmentDensityAdjust));
    end
    drawnow;
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(contourFig,outName,'png');
    
    % These files are sort of big, so don't always save.
    %
    % Probably could figure out what not to save and get 
    % a useful piece of the data if we wanted to.
    %save(outName);
    
    % Close windows
    close all
    
end

