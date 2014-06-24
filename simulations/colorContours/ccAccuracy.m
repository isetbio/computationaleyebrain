function [acc, err, sParams] = ccAccuracy(simParams, sParams)
%% function ccAccuracy(simParams)
%    find classification accuracy for color discrimination experiment
%    defined in simParams
%
%  Inputs:
%    simParams    - structure, contains all information about experiment,
%                   for more information about the structure, please refer 
%                   to setParameters (Name might be changed afterwards) and
%                   constructSimulationParameters
%    sParams      - contains precomputed values
%                   Required fields include
%                     .refRGB,  reference color RGB value in 0~1
%                     .display, ISET compatible dispaly structure
%                     .sensor,  human eye sensor structure
%                   Optional fields include
%                     .cbType, colorblind type
%                     .scenePixels, size of patches in pixels
%                     .refScene, ISET scene structure for reference color
%                     .refOI,    optical image for reference color
%                     .isRefVoltsImgComputed, bool, indicating whether
%                                volts image for reference color exists in
%                                sensor structure
%
%
%  Outputs:
%    acc          - classification accuracy
%    err          - standard deviation of accuracies in cross validation
%    staticValues - static value structure with optional fields computed,
%                   this is useful only when program is not run in
%                   parellel, otherwise, do not use this output
%
%  Example:
%
%
%  See also:
%    setParameters, 
%
%  DHB/BW/HJ ISETBIO Team, 2013

%% Check inputs
if nargin < 1, error('simulation parameter structure required'); end
if nargin < 2, error('static values structure required'); end

%  set number of frames to be used in train and test
try nFrames = sParams.nFrames; catch, nFrames = 3000; end
try display = sParams.display; catch, display = 'OLED-Sony'; end

if ~isfield(sParams, 'doSecondSiteNoise')
    sParams.doSecondSiteNoise = false;
end

%% Compute photon images for reference image
%  if scene for reference color is not set, create it
if ~isfield(sParams, 'refScene')
    refImg   = ones([sParams.scenePixels 3]);
    for i = 1 : 3
        refImg(:,:,i) = simParams.matchRGB(i);
    end
    sParams.refScene = sceneFromFile(refImgName, 'rgb', [], display);
end

%  if refOI is not computed, create and compute it
if ~isfield(sParams, 'refOI')
    sParams.refOI = oiCreate('wvf human');
end

%  if photon images is not computed, compute and store it
if ~isfield(sParams, 'isRefVoltsImgComputed') || ...
        ~sParams.isRefVoltsImgComputed
    sensor = sensorSet(sParams.sensor, 'sensor positions', ...
                        zeros(nFrames, 2));
    sParams.sensor = coneAbsorptions(sensor, sParams.refOI);
    sParams.isRefVoltsImgComputed = 1;
end

% get photon absorptions from each cone in the sensor array
% if doSecondSiteNoise is true, the sensor should contain photon
% absorptions with second site noise here. This is not reasonable to store
% it there, should change it to a cone structure later
refPhotons = sensorGet(sParams.sensor, 'photons');

if sParams.doSecondSiteNoise
    coneType = sensorGet(sParams.sensor, 'cone type');
    refPhotons = coneComputeSSNoise(refPhotons, coneType);
end
refPhotons = refPhotons(50:60, 50:60, :);

%% Compute photon images for match image
%  Compute match image color
if ~isfield(simParams, 'matchRGB')
    % Compute match value
    dir = [cos(simParams.cdAngle) sin(simParams.cdAngle) 0]';
    matchLMS = sParams.refLMS + simParams.nTestLevels * dir;
    
    simParams.matchRGB = coneContrast2RGB(sParams.display,...
        matchLMS, sParams.bgColor);
end
%  Create image for match color
if ~isfield(sParams, 'scenePixels')
    sParams.scenePixels = [64 64]; 
end
if isscalar(sParams.scenePixels)
    val = sParams.scenePixels;
    sParams.scenePixels = [val val];
end

% This creates the data that will be written out into the scene file
matchImg   = ones([sParams.scenePixels 3]);
for i = 1 : 3
    matchImg(:,:,i) = simParams.matchRGB(i);
end

% Create scene for match patch 
matchScene = sceneFromFile(matchImg, 'rgb', [], display);
matchScene = sceneSet(matchScene, 'fov', sceneGet(sParams.refScene,'fov'));

%  Compute cone samples
sensor = coneSamples(matchScene, nFrames, sParams.sensor, sParams.refOI);
matchPhotons = sensorGet(sensor, 'photons');
matchPhotons = double(matchPhotons);

%  Compute second site noise if needed
if sParams.doSecondSiteNoise
    matchPhotons = coneComputeSSNoise(matchPhotons, sParams.coneType);
end

%  Sample cone response in ROI
matchPhotons = matchPhotons(50:60, 50:60, :);


%% Classification
%  Set svm options
svmOpts = '-s 0 -q';

% We have matchPhotons - these are the test stimuli
% We have refPhotons   - these are the absorptions to the background
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
refPhotons   = RGB2XWFormat(refPhotons)';
matchPhotons = RGB2XWFormat(matchPhotons)';
acc = svmClassifyAcc(cat(1,refPhotons, matchPhotons), ...
                     labels, nFolds, 'svm', svmOpts);
err = acc(2);
acc = acc(1);


%% Tell the calling routine that we computed the refPhotons
sParams.isRefVoltsImgComputed = true;

end