%% s_ccQuickTest(parameterPreset)
%
% Illustrate the steps in producing a color threshold contour using
% ISETBIO.
%
% Requires:
%   ISETBIO
%     - Available on gitHub as 
%         https://github.com/isetbio/isetbio.git
%   PsychophysicsToolbox-3
%     - Available on gitHub as 
%         https://github.com/Psychtoolbox-3/Psychtoolbox-3.git
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

%%
s_initISET

%% Main parameters for background
monitorName = 'LCD-Apple.mat';
wave        = 380:4:780;
refColor    = [0.5 0.5 0.5];
bgColor     = [0.5 0.5 0.5];

% Init static parameter structure
nFrames = 1500;
scenePixels = [64 64];
display = displayCreate(monitorName);
refLMS = RGB2ConeContrast(display, refColor, bgColor);

%% Create reference color patch
%  Assume that display has been linearized before experiment
refImage = ones([scenePixels 3]);
for ii = 1 : 3
    refImage(:,:,ii) = refImage(:,:,ii) * refColor(ii);
end

%% Initiate the human optics
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf = wvfSet(wvf,'zcoeffs',sample_mean);
wvf = wvfComputePSF(wvf);
oi  = wvf2oi(wvf,'shift invariant');
refOI = oi; matchOI = oi;

%% Build the background scene and oi from the image file
refScene = sceneFromFile(refImage,'rgb',[], monitorName, wave);
refScene = sceneSet(refScene,'h fov', 0.5);
refOI    = oiCompute(refOI, refScene);

%% Create human cone
coneDensity = [.1 .6 .2 .1];

% create cone structure
cone = coneCreate;
cone = coneSet(cone, 'exp time', 0.05);
cone = coneSet(cone, 'h fov', sceneGet(refScene,'hfov'), refScene, refOI);
cone = coneSet(cone, 'density', coneDensity);
cone = coneSet(cone, 'frames per position', nFrames);

%% Compute reference absorptions
cone = coneCompute(cone, refOI);
refPhotons = double(coneGet(cone, 'photons'));
refPhotons = refPhotons(50:60, 50:60, :); % Should change to a roiLoc
refPhotons   = RGB2XWFormat(refPhotons)';

%% Simulate under each conditions
% Loop over and compute classification accuracy
% Currently, the classification accuracy is a little too high. This is
% because we didn't add any second site noise here. We hope to add the
% second site noise in rgc level

angle = (0:45:179)*pi/180;
for i = 1 : length(angle)
    dir = [cos(angle(i)) sin(angle(i)) 0]';
    for nTestLevels = 0.002 + 0.002 * (1:5);
        % Compute match color
        matchLMS = refLMS + nTestLevels * dir;
        matchRGB = coneContrast2RGB(display, matchLMS, bgColor);
        
        % Create match scene
        matchImg   = ones([scenePixels 3]);
        for k = 1 : 3
            matchImg(:,:,k) = matchRGB(k);
        end
        matchScene = sceneFromFile(matchImg,'rgb', [], monitorName, wave);
        
        % Compute match absorptions
        matchOI = oiCompute(matchScene, matchOI);
        cone = coneCompute(cone, matchOI);
        matchPhotons = double(coneGet(cone, 'photons'));
        matchPhotons = matchPhotons(50:60, 50:60, :);
        
        % Classification
        svmOpts = '-s 0 -q';
        
        nFolds = 10;
        labels = [ones(nFrames,1); -1*ones(nFrames,1)];
        matchPhotons = RGB2XWFormat(matchPhotons)';
        acc = svmClassifyAcc(cat(1,refPhotons, matchPhotons), ...
            labels, nFolds, 'svm', svmOpts);
        acc = acc(1);
        
        result(i,round(nTestLevels/0.002)) = acc;
        % Show debug information
        fprintf('Simulation %d: Angle - %d, Level - %f...', i, ...
            round(angle(i)*180/pi), nTestLevels);
        fprintf('Acc - %.2f...Done!\n', acc);
    end
end

    
