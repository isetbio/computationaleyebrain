% colorContours
%
%  * Read in an RGB background value from a display
%  * Read in a color direction (RGB? LMS?)
%  * Find the threshold in that direction based on linear discriminant or
%       svm or something
%  * Plot it.
%
% Notes: We are currently simulating a one-interval experiment and
% getting the fraction correct averaged over a mixture of blank and
% test trials.  Since the data we are interested in is from TAFC, we
% should probably handle this differently in the long run.
%
% 8/2/13  DHB/BW Our excellent adventure commences.
% 8/3/12  DHB/BW Tune up and add classifier

%% Clear out the junk.  Remember where you are.
s_initISET
talkD = pwd;
saveFlag = 0;           % Don't save results

%% Get stats path off of the path, temporarily.
%
% This avoids a name space conflict with the libsvn
% routines that we use here.
if (exist('RemoveMatchingPaths','file'))
    path(RemoveMatchingPaths(path,'stats'));
end

%%  Get display spectra
d = displayCreate('LCD-Apple');
displaySpd = displayGet(d,'spd');
w = displayGet(d,'wave');
% vcNewGraphWin; plot(w,displaySpd)

%% Create a scene with the background spectrum
%
% To patch into scene creation routines, we create an image
% on disk with spatially uniform RGB values.  These are then treated
% as frame buffer values, and will be raised to the 2.2 (gamma uncorrected)
% by the scene creation routine.
scenePixels = 64;
backRGBValue = 0.5;
gammaValue = 2.2;
backImg = ones(scenePixels,scenePixels,3)*(backRGBValue)^(1/gammaValue);
imwrite(backImg,'backFile.png','png');
sceneB = sceneFromFile('backFile.png','rgb',[],'LCD-Apple.mat',w);
sceneB = sceneSet(sceneB,'name','background');
sceneB = sceneSet(sceneB,'fov',2);
vcAddAndSelectObject(sceneB);
%sceneWindow;

%% Create standard human polychromatic PSF
% using wavefront tools, and make it an
% iset OI thingy.
wvf = wvfCreate('wave',w);
pDiameter = 3;
sample_mean = wvfLoadThibosVirtualEyes(pDiameter);
wvf = wvfSet(wvf,'zcoeffs',sample_mean);
wvf = wvfComputePSF(wvf);
oiD = wvf2oi(wvf,'shift invariant');
% vcNewGraphWin; plotOI(oiD,'psf')

%%  Create a sample cone mosaic
% and take a look at it.
params.sz = [128,192];
params.rgbDensities = [0.1 .6 .2 .1]; % Empty, L,M,S
params.coneAperture = [3 3]*1e-6;     % In meters
pixel = [];
cSensor = sensorCreate('human',pixel,params);
%sensorConePlot(cSensor);
cSensor = sensorSet(cSensor,'exp time',0.050);
    
%% Get spectrum of background directly,
% mostly for debugging.
backRGB = [backRGBValue backRGBValue backRGBValue]';
backSpd = displaySpd*backRGB;
% vcNewGraphWin; plot(w,stockman);

%% Create a test vector in a specified
% direction in cone space.
%
% To do this, we pick a test color direction
% and find the RGB values corresponding to it.
% Then we scale and add these to the background.
%
% Define test direction in cone excitation space
testLMS = [1 0 0]';

% Get the matrix from RGB to cone space
T_stockman = vcReadSpectra('stockman',w);  % Energy
rgb2cones = T_stockman'*displaySpd;

% Compute the RGB direction and scale so that it
% reaches to the edge of the gamut.  Note that
% this isn't a robust caclulation because it works
% only when the background value is 0.5.  We're 
% probably going to convert this into contrast units
% at some point anyway, so we don't bother to do it
% right.  See PTB routine MaximizeGamutContrast for
% a more general computation.
if (backRGBValue ~= 0.5)
    error('Stimulus scaling only does what we wanted if background RGB values are all 0.5');
end
testRGB = inv(rgb2cones)*testLMS;
testRGB = testRGB/max(testRGB(:));
 
%% Pass the background through the optics
backOiD = oiCompute(oiD,sceneB);

%% Compute noise free background sensor image
backCSensorNF = sensorComputeNoiseFree(cSensor,backOiD);

%% Loop over a set of test levels and get classifier
% performance for each.
nLevels = 6;
testLevels = linspace(0,1,nLevels);
for t = 1:length(testLevels)
    % Set test level
    fprintf('Calculations for test level %d of %d\n',t,length(testLevels));
    testLevel = testLevels(t);
    testRGBForThisLevel = (backRGB + testLevel*0.5*testRGB).^(1/gammaValue);
    
    % Make test scene, in same fashion as we made the
    % background scene.
    tmp = ones(scenePixels,scenePixels,3);
    [tmp,row,col] = RGB2XWFormat(tmp);
    tmp = tmp*diag(testRGBForThisLevel(:));
    testImg = XW2RGBFormat(tmp,row,col);
    imwrite(testImg,'testFile.png','png');
    sceneT = sceneFromFile('testFile.png','rgb',[],'LCD-Apple.mat',w);
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
    % of the poisson/sensor noise.
    %
    % Still need to control the sources of sensor noise.
    nDraws = 100;
    noiseType = 1;
    slot = [2 3 4];   % Typical human 1621 case.  1 empty, 6 L, 2 M, 1 S
    nSensorClasses = length(slot);
    testCSensorNF = sensorComputeNoiseFree(cSensor,testOiD);
    backVoltImage = sensorComputeSamples(backCSensorNF,nDraws,noiseType);
    testVoltImage = sensorComputeSamples(testCSensorNF,nDraws,noiseType);  
    for k = 1:nDraws
        for ii = 1:nSensorClasses
            backCSensorTemp = sensorSet(backCSensorNF,'volts',backVoltImage(:,:,k));
            %backCSensorTemp = backCSensorNF;
            backSensorVals{k,ii} = sensorGet(backCSensorTemp,'electrons',slot(ii));
            testCSensorTemp = sensorSet(testCSensorNF,'volts',testVoltImage(:,:,k));
            %testCSensorTemp = testCSensorNF;
            testSensorVals{k,ii} = sensorGet(testCSensorTemp,'electrons',slot(ii));
        end
    end
    
    %% Pop last instantions into windows for viewing etc.
    backCSensor = sensorSet(backCSensorNF,'name','Background');
    vcAddAndSelectObject(backCSensorNF); sensorWindow;
    
    testCSensor = sensorSet(testCSensorNF,'name','Test');
    vcAddAndSelectObject(testCSensorNF); sensorWindow;
    
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
    backVec = zeros(nUseAll,nDraws);
    testVec = zeros(nUseAll,nDraws);
    typeVec = zeros(nUseAll,1);
    oneConeEachClassIndices = zeros(nSensorClasses,1);
    for k = 1:nDraws
        if (rem(k,10) == 0)
            fprintf('\tGetting cone catches for draw %d of %d\n',k,nDraws);
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
    f = vcNewGraphWin; hold on;
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
    %   leads to good performance on the training set, at least
    %   with the default choice of fgamma.
    %
    %   Using a linear kernal works quite well in my initial test.
    svmOpts = '-s 0 -t 0';
    svmModel = svmtrain(trainingLabels, trainingData, svmOpts);
    [svmTrainingPredictedLabels] = svmpredict(trainingLabels, trainingData, svmModel);
    [svmValidatePredictedLabels] = svmpredict(validateLabels, validateData, svmModel);
    trainingFractionCorrect = length(find(svmTrainingPredictedLabels == trainingLabels))/length(trainingLabels);
    validateFractionCorrect = length(find(svmValidatePredictedLabels == validateLabels))/length(validateLabels);
    fprintf('Classifier percent correct: %d (training data), %d (validation data)\n',round(100*trainingFractionCorrect),round(100*validateFractionCorrect));
    
    % Plot the training and test data.  We'll plot the distribution of responses for one cone
    % in each sensor class, with the distribution taken over our resampling of each
    % scene by the mosaic.
    f = vcNewGraphWin; hold on;
    sym = {'b.','g.','r.','c.','k.'};
    az = 65.5; el = 30;
    index = find(svmTrainingPredictedLabels == trainingLabels);
    plot3(trainingData(index,1), ...
        trainingData(index,2), ...
        trainingData(index,3),'go','MarkerFaceColor','g');
    index = find(svmTrainingPredictedLabels ~= trainingLabels);
    plot3(trainingData(index,1), ...
        trainingData(index,2), ...
        trainingData(index,3),'ro','MarkerFaceColor','r');
    
    index = find(svmValidatePredictedLabels == validateLabels);
    plot3(validateData(index,1), ...
        validateData(index,2), ...
        validateData(index,3),'go');
    index = find(svmValidatePredictedLabels ~= validateLabels);
    plot3(validateData(index,1), ...
        validateData(index,2), ...
        validateData(index,3),'ro');
    xlabel('L-absorptions'); ylabel('M-Absorptions'); zlabel('S-absorptions'); axis square; grid on
    
    %% Store out performance measure
    nCorrectResponses(t) = length(find(svmValidatePredictedLabels == validateLabels));
    nTotalResponses(t) = length(validateLabels);
    fractionCorrect(t) = validateFractionCorrect;
end

%% Plot and fit psychometric function, extract threshold
%
% Requires Palamedes toolbox
%
% Plot data
figure; clf; hold on
plot(testLevels,fractionCorrect,'ro','MarkerSize',10,'MarkerFaceColor','r');

% Fit and find thresholdﬂ
criterionCorr = 0.82;
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
testLevelsInterp = linspace(testLevels(1),testLevels(end),100);
probCorrInterp = PF(paramsValues,testLevelsInterp);
threshPal = PFI(paramsValues,criterionCorr);
plot(testLevelsInterp,probCorrInterp);
plot([threshPal threshPal],[0.5 criterionCorr],'g');
plot([testLevels(1) threshPal],[criterionCorr criterionCorr],'g');


% Print threshold 
fprintf('%d%% correct threshold is %0.1f\n',round(100*criterionCorr),threshPal);


