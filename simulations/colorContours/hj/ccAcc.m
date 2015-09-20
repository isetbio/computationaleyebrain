function [acc, err, params] = ccAcc(rColor, mColor, params)
%% Compute classification accuracy
%    [acc, err, params] = ccAcc(rColor, mColor, params)
%
%  Inputs:
%    rColor   - 3x1 vector, reference color in RGB
%    mColor   - 3x1 vector, test/match color in RGB
%    params   - parameter structure, could include
%               .nSamples, number of samples per color
%               .sensor, human eye sensor structure
%               .oi, human optics structure
%               .d, display structure
%               .svmOpts, params for svm classification
%               .cropSz, croping size for classification
%
%  Outputs:
%    acc     - classification accuracy
%    err     - standard deviation of classification accuracy
%
%  Example:
%    [acc, err] = ccAcc([0.5 0.5 0.5], [0.49 0.49 0.5]);
%
%  See also:
%    ccThreshold, svmClassifyAcc
%
%  (HJ), ISETBIO TEAM, 2014

%% Check inputs
if notDefined('rColor'), error('reference color required'); end
if notDefined('mColor'), error('match color required'); end
if notDefined('params'), params = []; end

%% Set up scene
if isfield(params, 'd'), d = params.d; 
else
    d = displayCreate('OLED-Sony');
    d = displaySet(d, 'gamma', 'linear');
end
if isfield(params, 'sceneSz'), sz = params.sceneSz; else sz = 64; end

% Create reference scene
rColor = reshape(rColor, [1 1 3]);
rScene = sceneFromFile(rColor(ones(sz, 1), ones(sz, 1), :), 'rgb', [], d);
rScene = sceneSet(rScene, 'fov', 5);

% Create test/match scene
mColor = reshape(mColor, [1 1 3]);
mScene = sceneFromFile(mColor(ones(sz, 1), ones(sz, 1), :), 'rgb', [], d);
mScene = sceneSet(mScene, 'fov', 5);

%% Compute optical image
if isfield(params, 'oi')
    oi = params.oi;
else
    oi = oiCreate('wvf human');
end

rOI = oiCompute(rScene, oi);
mOI = oiCompute(mScene, oi);

%% Compute cone absorptions
if isfield(params, 'sensor')
    sensor = params.sensor;
else
    if isfield(params, 'sensorSz'), sz = params.sensorSz;
    else sz = [45, 45]; end
    
    sensor = sensorCreate('human', [], params.cone);
    % fov = sceneGet(rScene, 'fov');
    % sensor = sensorSetSizeToFOV(sensor, fov, rScene, oi);
    sensor = sensorSet(sensor, 'size', sz);
end

% set up number of samples
if isfield(params, 'nSamples')
    nSamples = params.nSamples;
else
    nSamples = 5000;
    params.nSamples = nSamples;
end

% disable eye movement
sensor = sensorSet(sensor, 'sensor positions', zeros(nSamples, 2));
sensor = coneAbsorptions(sensor, rOI);
rPhotons = sensorGet(sensor, 'photons');

sensor = coneAbsorptions(sensor, mOI);
mPhotons = sensorGet(sensor, 'photons');

%% Add second site noise (cone opponency)
% coneType = sensorGet(sensor, 'cone type');
% cg = sensorGet(sensor, 'conversion gain');
% rVolts = coneComputeSSNoise(rVolts / cg, coneType) * cg;
% mVolts = coneComputeSSNoise(mVolts / cg, coneType) * cg;
% rVolts = coneComputeCenterSurround(rVolts);
% mVolts = coneComputeCenterSurround(mVolts);

%% SVM classification
nFolds = 5;
if isfield(params, 'svmOpts')
    svmOpts = params.svmOpts;
else
    svmOpts = [];
    params.svmOpts = [];
end

labels = [ones(nSamples,1); -1*ones(nSamples,1)];
rPhotons = RGB2XWFormat(rPhotons)';
mPhotons = RGB2XWFormat(mPhotons)';

acc = svmAcc([rPhotons; mPhotons], labels, nFolds, 'linear', svmOpts);
err = nan;

end