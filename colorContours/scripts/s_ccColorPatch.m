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
%    doOneSimulation
%
%  DHB/BW/HJ (c) ISETBIO Team, 2013

%% Init ISET
s_initISET

%% Main parameters for background
monitorName = 'LCD-Apple.mat';
wave        = 380:4:780;
refColor    = [0.5 0.5 0.5];
staticValues.refColor = refColor;

%% Create reference color patch
%  Assume that display has been linearized before experiment
staticValues.display = displayCreate(monitorName);
refImage = ones(64,64,3);
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

%% Build the background scene and oi from the image file.
refScene = sceneFromFile(refFile,'rgb',[],monitorName,wave);
refScene = sceneSet(refScene,'h fov',2);
staticValues.refScene = refScene;
staticValues.refOI = oiCompute(staticValues.refOI, refScene);
vcAddAndSelectObject(staticValues.refOI); oiWindow;

%% Create human sensor
coneDensity = [.1 .6 .2 .1];
sensor = sensorCreate('human');
sensor = sensorSet(sensor,'exp time',0.05);
[sensor,xy,coneType] = sensorCreateConeMosaic(sensor, [], coneDensity);
sensor = sensorSetSizeToFOV(sensor,sceneGet(refScene,'hfov'), ...
    refScene, staticValues.refOI);
sensor = sensorCompute(sensor,staticValues.refOI);
staticValues.sensor = sensor;
vcAddAndSelectObject(staticValues.sensor); sensorWindow('scale',1);

%% Create simulation parameters
[theParams, staticParams] = setParameters('QuickTest');
simParams = constructSimulationParameters(theParams, staticParams);

%% Simulate under each conditions
% Try open matlabpool
% matlabpool open 4
%  Loop over and compute classification accuracy
for curSim = 2 % 1 : length(simParams)
    % Compute match value
    params     = simParams(curSim);
    
    % This isn't right.  Need to move into LMS space.
    params.matchRGB = refColor + params.DO_TAFC_CLASSIFIER .* ...
                                            0.0058*[cos(params.cdAngle),sin(params.cdAngle),0];
    % Do simulation
    simResults(curSim) = ccAccuracy(params, staticValues);
end

%  Close matlabpool
% matlabpool close
%% Plot Result