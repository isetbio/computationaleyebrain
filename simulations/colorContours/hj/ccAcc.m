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
    
    sensor = sensorCreate('human', params.cone);
    % fov = sceneGet(rScene, 'fov');
    % sensor = sensorSetSizeToFOV(sensor, fov, rScene, oi);
    sensor = sensorSet(sensor, 'size', sz);
    % sensor.human.macular.density = 0;
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

%% Compute outer-segment and retina ganglion cell response
% osL = osCreate('linear');
% osL = osSet(osL, 'patch size', sensorGet(sensor, 'width', 'um'));
% osL = osSet(osL, 'time step', sensorGet(sensor, 'time interval', 'sec'));
% 
% sensor = sensorSet(sensor, 'photons', rPhotons(:, :, 1));
% osL = osCompute(osL, sensor);
% 
% irParams.name      = 'Macaque inner retina 1'; % This instance
% irParams.eyeSide   = 'left';   % Which eye
% irParams.eyeRadius = 0;        % Radius in mm
% irParams.eyeAngle  = 90;       % Polar angle in degrees
% 
% innerRetina = irCreate(osL, irParams);
% innerRetina.mosaicCreate('model','glm','type','on midget');
% 
% for ii = 1 : size(rPhotons, 3)
%     photons = repmat(rPhotons(:, :, ii), [1 1 50]);
%     sensor = sensorSet(sensor, 'photons', photons);
%     osL = osCompute(osL, sensor);
%     
%     innerRetina = irCompute(innerRetina, osL);
% end

%% Add second site noise (cone opponency)
coneType = sensorGet(sensor, 'cone type');
rPhotons = coneComputeSSNoise(rPhotons, coneType);
mPhotons = coneComputeSSNoise(mPhotons, coneType);
rPhotons = coneComputeCenterSurround(rPhotons);
mPhotons = coneComputeCenterSurround(mPhotons);

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