function [acc, err] = ccAcc(rColor, mColor, params)
%% Compute classification accuracy
%    ccAcc(rColor, mColor, params)
%
%  Inputs:
%    rColor   - 3x1 vector, reference color in RGB
%    mColor   - 3x1 vector, test/match color in RGB
%    params   - parameter structure, could include
%               - nSamples, number of samples per color
%               - sensor, human eye sensor structure
%               - oi, human optics structure
%               - d, display structure
%               - svmOpts, params for svm classification
%               - cropSz, croping size for classification
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
if isfield(params, 'd'), d = params.d; else d = 'OLED-Sony'; end
if isfield(params, 'sceneSz'), sz = params.sceneSz; else sz = 64; end

% Create reference scene
rColor = reshape(rColor, [1 1 3]);
rScene = sceneFromFile(rColor(ones(sz, 1), ones(sz, 1), :), 'rgb', [], d);
rScene = sceneSet(rScene, 'fov', 0.5);

% Create test/match scene
mColor = reshape(mColor, [1 1 3]);
mScene = sceneFromFile(mColor(ones(sz, 1), ones(sz, 1), :), 'rgb', [], d);
mScene = sceneSet(mScene, 'fov', 0.5);

%% Compute optical image
if isfield(params, 'oi'), oi = params.oi; else oi = oiCreate('human'); end
rOI = oiCompute(rScene, oi);
mOI = oiCompute(mScene, oi);

%% Compute cone absorptions
if isfield(params, 'sensor')
    sensor = params.sensor;
else
    sensor = sensorCreate('human');
    fov = sceneGet(rScene, 'fov');
    sensor = sensorSetSizeToFOV(sensor, fov, rScene, oi);
end

% set up number of samples
if isfield(params, 'nSamples')
    nSamples = params.nSamples;
else
    nSamples = 3000;
end
sensor = sensorSet(sensor, 'sensor positions', zeros(nSamples, 2));

% compute adapted cone samples
sensor = coneAbsorptions(sensor, rOI); % reference cone absorptions
[~, rVolts] = coneAdapt(sensor);       % reference adapted data

sensor = coneAbsorptions(sensor, mOI); % test cone absorptions
[~, mVolts] = coneAdapt(sensor);       % test adapted data

%% Add second site noise (cone opponency)
coneType = sensorGet(sensor, 'cone type');
cg = sensorGet(sensor, 'conversion gain');
rVolts = coneComputeSSNoise(rVolts / cg, coneType) * cg;
mVolts = coneComputeSSNoise(mVolts / cg, coneType) * cg;

%% Crop from center
if isfield(params, 'cropSz'), cropSz = params.cropSz;
else cropSz = 24; end
rVolts = getMiddleMatrix(rVolts, cropSz); % Get patch from center
mVolts = getMiddleMatrix(mVolts, cropSz); % Get patch from center

%% SVM classification
nFolds = 10;
if isfield(params, 'svmOpts')
    svmOpts = params.svmOpts;
else
    svmOpts = [];
end

labels = [ones(nSamples,1); -1*ones(nSamples,1)];
rVolts = RGB2XWFormat(rVolts)';
mVolts = RGB2XWFormat(mVolts)';

acc = svmClassifyAcc(cat(1, rVolts, mVolts), ...
                        labels, nFolds, 'ranksvm', svmOpts);
err = acc(2);
acc = acc(1);

end