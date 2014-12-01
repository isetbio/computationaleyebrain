function [acc, err] = coPixelVisibilityAcc(d, tDist, params)
%% function coPixelVisibilityAcc(params)
%    Compute classification accuracy of white display pixels vs uniform
%    scene on certain display and viewing distance
%
%  Inputs:
%    d     : isetbio display structure, see displayCreate for more detail
%    tDist : scalar, viewing distance to be tested
%    params: parameter structure, could include:
%      coneDensity   - human cone spatial density, default 0.6 0.3 0.1
%      expTime       - eye integration time
%      nFrames       - number of samples to be used by the classifier
%      cropSz        - number of cones to be used in classification
%
%  Outputs:
%    acc:     classification accuray
%    err:     classification standard deviation
%
%  See also:
%    coPixelVisibilityThreshold, s_coPixelVisibility_Proclus
%
% (HJ) ISETBIO TEAM, 2014

%% Check inputs & init parameters
%  check inputs
if notDefined('d'), error('display not defined'); end
if notDefined('params'), params = []; end

%  init simulation parameters
try density = params.coneDensity; catch, density = [0 .6 .3 .1]; end
try expTime = params.expTime;     catch, expTime = 0.05; end
try nFrames = params.nFrames;     catch, nFrames = 3000; end
try cropSz  = params.cropSz;      catch, cropSz  = [12 12]; end

% init parameters for classifiers
svmOpts = '-s 0 -q';
nFolds = 10;         % folds of cross-validation
labels = [ones(nFrames,1); -ones(nFrames,1)];

%% Generate scene
%  Set display to white and generate scene
imgSz  = 20; doSub = true;
d = displaySet(d, 'viewing distance', tDist);
sceneD = sceneFromFile(ones(imgSz), 'rgb', [], d, [], doSub);

%  Generate control scene by spatially averaging radiance
p  = mean(mean(sceneGet(sceneD, 'photons')));
sz = sceneGet(sceneD, 'size');
sceneU = sceneSet(sceneD, 'photons', repmat(p, [sz 1]));

%% Compute optical image
%  Create human optics
oi = oiCreate('wvf human');

%  Compute irradiance map for scenes
oiD = oiCompute(sceneD, oi);
oiU = oiCompute(sceneU, oi);

%% Compute cone absorptions
%  Create human sensor
params.humanConeDensity = density;
params.sampTime = 0.001; % 1 ms
sensor = sensorCreate('human', [], params);
sensor = sensorSet(sensor, 'exp time', expTime);
sensor = sensorSet(sensor, 'sample time interval', params.sampTime); % 1 ms

% Generate eye movement sequence
% params.emType = [1 0 0]; % include tremor only
% params.nSamples = nFrames * expTime / sampTime;
% sensor = eyemoveInit(sensor, params);
sensor = sensorSet(sensor, 'sensorpositions', zeros(nFrames, 2));

% Compute cone absorptions
sensor = coneAbsorptions(sensor, oiD);
refPhotons = getMiddleMatrix(sensorGet(sensor, 'photons'), cropSz);
refPhotons = RGB2XWFormat(refPhotons);

sensor = coneAbsorptions(sensor, oiU);
matchPhotons = getMiddleMatrix(sensorGet(sensor, 'photons'), cropSz);
matchPhotons = RGB2XWFormat(matchPhotons);

%% Classification
accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
    labels, nFolds, 'svm', svmOpts);

err = accuracy(2);
acc = accuracy(1);

end