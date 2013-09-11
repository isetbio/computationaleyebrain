function [acc, err, staticValues] = ccAccuracy(simParams, staticValues)
%% function ccAccuracy(simParams)
%    find classification accuracy for color discrimination experiment
%    defined in simParams
%
%  Inputs:
%    simParams    - structure, contains all information about experiment,
%                   for more information about the structure, please refer 
%                   to setParameters (Name might be changed afterwards) and
%                   constructSimulationParameters
%    staticValues - contains precomputed values
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

%% Compute photon images for reference image
%  if scene for reference color is not set, create it
if ~isfield(staticValues, 'refScene')
end

%  if refOI is not computed, create and compute it
if ~isfield(staticValues, 'refOI')
end

%  if photon images is not computed, compute and store it
if ~isfield(staticValues, 'isRefVoltsImgComputed') || ...
        ~staticValues.isRefVoltsImgComputed
end

% get photon images from sensor
refSensor  = staticValues.sensor; % store current state for return value
refPhotons = sensorGet(staticValues.sensor,'photons');

%% Compute photon images for match image
%  Compute match image color
if ~isfield(simParams, 'matchRGB')
end
%  Create image for match color
if ~isfield(staticValues, 'scenePixels')
    staticValues.scenePixels = [64 64]; 
end
if isscalar(staticValues.scenePixels)
    val = staticValues.scenePixels;
    staticValues.scenePixels = [val val];
end

matchImg   = ones([staticValues.scenePixels 3]);
for i = 1 : 3
    matchImg(:,:,i) = simParams.matchRGB(i);
end

% Create scene for match patch
% should change the random number to some simulation 
matchImgName = fullfile(isetbioRootPath,'tmp',['matchImg_' ...
                        num2str(randi(1e6)) '.png']);
imwrite(matchImg, matchImgName);

matchScene = sceneFromFile(matchImgName, 'rgb', [], staticValues.display);

% clean up
delete(matchImgName);

%  Compute cone samples
matchPhotons = coneSamples(matchScene, 10, staticValues.sensor, ...
    staticValues.refOI);
%% Classification
[acc, err] = getClassificationAccuray(refPhotons, matchPhotons);

%% Set values for output
staticValues.isRefVoltsImgComputed = true;
staticValues.sensor = refSensor;

end