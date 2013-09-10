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
refColor    = 0.5;             % Linear value of display primaries

%% Create reference color patch
%  Assume that display has been linearized before experiment
d      = displayCreate(monitorName);
refImage = ones(128,128,3);
for ii=1:3
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
oiB    = wvf2oi(wvf,'shift invariant');

%% Build the background scene and oi from the image file.
bScene = sceneFromFile(refFile,'rgb',[],monitorName,wave);
oiB    = oiCompute(oiB,bScene);
% vcAddAndSelectObject(oiB); oiWindow;

%% Create simulation parameters
[theParams,staticParams] = setParameters('QuickTest');
simParams = constructSimulationParameters(theParams,staticParams);

%% Simulate under each conditions
%  Try open matlabpool
for curSim = 1 : length(simParams)
    % Compute match value
    staticComputedValues;
    % Do simulation
    simResults(curSim) = doOneSimulation(simParams(curSim),staticParams, ...
        runtimeParams,staticComputedValues);
end

%% Plot Result