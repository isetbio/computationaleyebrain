%% s_ccColorPatch
%
%  Illustrate the steps in producing a color threshold contour using
%  ISETBIO.
%
% Requires:
%   ISETBIO
%     - Available on gitHub:
%       https://github.com/isetbio/isetbio.git
%   PsychophysicsToolbox-3
%     - Available on gitHub:
%       https://github.com/Psychtoolbox-3/Psychtoolbox-3.git
%
%  Process:
%    1. Set up a reference color (0.5 gray)
%    2. Set color test directions
%    3. Use ISETBIO and PTB to create isomerizations to create
%       reference color and reference + s*test for some scalar, s.
%    4. Classify as a function of s for each test direction.
%    5. Plot the psychometric functions and the detection contours
%
%  See also:
%    ccAccuracy
%
%  DHB/BW/HJ (c) ISETBIO Team, 2013

%% Init ISET
s_initISET

%% Main parameters for background
monitorName = 'LCD-Apple.mat';
wave        = 380:4:780;
refColor    = [0.5 0.5 0.5];
bgColor     = [0.5 0.5 0.5];

% Init static parameter structure
staticValues.refColor = refColor;
staticValues.doSecondSiteNoise = true;
staticValues.nFrames = 1500;
staticValues.scenePixels = [64 64];

%% Create reference color patch
%  Assume that display has been linearized before experiment
staticValues.display = displayCreate(monitorName);
refImage = ones([staticValues.scenePixels 3]);
for ii = 1 : 3
    refImage(:,:,ii) = refImage(:,:,ii) * refColor(ii);
end
refFile = fullfile(isetbioRootPath,'tmp','refFile.png');
imwrite(refImage,refFile);

%% Initiate the human optics
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);
wvf    = wvfComputePSF(wvf);
staticValues.refOI    = wvf2oi(wvf,'shift invariant');

%% Build the background scene and oi from the image file
refScene = sceneFromFile(refImage,'rgb',[],monitorName,wave);
refScene = sceneSet(refScene,'h fov', 0.5);
staticValues.refScene = refScene;
staticValues.refOI = oiCompute(staticValues.refOI, refScene);
% vcAddAndSelectObject(staticValues.refOI); oiWindow;

%% Create human sensor
coneDensity = [.1 .6 .2 .1];
sensor = sensorCreate('human');
sensor = sensorSet(sensor,'exp time',0.05);
sensor = sensorSetSizeToFOV(sensor,sceneGet(refScene,'hfov'), ...
    refScene, staticValues.refOI);
[sensor, ~, coneType] = sensorCreateConeMosaic(sensor, ...
                          sensorGet(sensor, 'size'), coneDensity);
sensor = sensorCompute(sensor,staticValues.refOI);

staticValues.sensor   = sensor;
staticValues.coneType = coneType;
% vcAddAndSelectObject(staticValues.sensor); sensorWindow('scale',1);

%% Create simulation parameters
[theParams, staticParams] = setParameters('QuickTest');
simParams = constructSimulationParameters(theParams, staticParams);

%% Simulate under each conditions
% Try open matlabpool
% matlabpool open 4
% Loop over and compute classification accuracy
refLMS   = RGB2ConeContrast(staticValues.display, refColor, bgColor);
for curSim = 1 : length(simParams)
    % Compute match value
    params       = simParams(curSim);
    staticParams = staticValues; % Make a local copy, parfor only
    
    % Show some debug information here, should be deleted when computed
    % parallel
    fprintf('Simulation %d: Angle - %d, Level - %f...', curSim, ...
        round(params.cdAngle*180/pi), params.nTestLevels);
    
    dir = [cos(params.cdAngle) sin(params.cdAngle) 0]';
    matchLMS = refLMS + params.nTestLevels * dir;
    
    params.matchRGB = coneContrast2RGB(staticParams.display,...
                                       matchLMS, bgColor);
    % Do simulation
    [simResults(curSim),~,staicValues] = ccAccuracy(params, staticParams);
    
    % Show debug information here, should be deleted in parallel
    % computation
    fprintf('Acc - %.2f...Done!\n', simResults(curSim));
end

save deleteMe.mat simResults staticValues simParams
%  Close matlabpool
%  matlabpool close

%%  Re-organize data
results.angle   = unique([simParams.cdAngle]);
results.level  = unique([simParams.nTestLevels]);
results.predAcc = zeros(length(results.angle), length(results.level));

for i = 1 : length(simParams)
    angleIndx = find(results.angle == simParams(i).cdAngle, 1);
    levelIndx = find(results.level == simParams(i).nTestLevels, 1);
    results.predAcc(angleIndx, levelIndx) = simResults(i);
end


%% Weibull fit and plot
%  init plot figure
figure; hold on;

%  fit for every direction and plot weibull curve
pCorrect = 0.51:0.01:0.99;
results.ccThresh = zeros(length(results.angle), 1);
results.threshColor = zeros(2, length(results.angle));

for curDir = 1 : length(results.angle)
    dir = [cos(results.angle(curDir)) sin(results.angle(curDir))];
    
    [alpha,beta,~] = FitWeibAlphTAFC(results.level, ...
        results.predAcc(curDir,:) * staticValues.nFrames, ...
        (1-results.predAcc(curDir, :)) * staticValues.nFrames,[],2.2);
    results.ccThresh(curDir) = FindThreshWeibTAFC(0.75,alpha,beta);
    results.threshColor(:, curDir) = refLMS(1 : 2)' + ...
        results.ccThresh(curDir)*dir;
    
    % Plot
    threshX = alpha*(-log(2*(1-pCorrect))).^(1/beta);
    plot(threshX,pCorrect);
    plot(results.level,results.predAcc(curDir,:),'or');
end

%  plot threshold data and fitted ellipse
figure; hold on;

plot(results.threshColor(1,:), results.threshColor(2,:), 'ro');
[zg, ag, bg, alphag] = fitellipse(results.threshColor);
plotellipse(zg, ag, bg, alphag, 'b--')