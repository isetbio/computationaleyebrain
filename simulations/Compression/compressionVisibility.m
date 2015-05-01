function [acc, err] = compressionVisibility(refImg, testImg, d, params)
%% Computing classification accuracy for compressed images
%    function compressionVisibility(refImg, testImg, d, [params]);
%
%  Inputs:
%    refImg  - reference RGB image
%    testImg - compressed RGB image
%    d       - display structure, see displayCreate
%    params  - parameter structure, can include
%            .vd      - viewing distance
%            .cone    - human cone structure, see coneCreate
%            .nFrames - number of frames
%            .nFolds  - cross-validation folds
%            .svmOpts - parameters for libsvm
%
%  Output:
%    acc - classification accuracy
%    err - standard deviation of acc in cross validation
%
%  See also:
%    s_compressionVisibility
%
%  (HJ) ISETBIO TEAM, 2015

%% Check inputs and init
if notDefined('refImg'), error('reference image required'); end
if notDefined('testImg'), error('test image required'); end
if notDefined('d'), error('display structure required'); end
if notDefined('params'), params = []; end

try vd = params.vd; catch, vd = 1; end
try nFrames = params.nFrames; catch, nFrames = 3000; end
try cone = params.cone; catch, cone = coneCreate; end
try nFolds = params.nFolds; catch, nFolds = 5; end
try svmOpts = params.svmOpts; catch, svmOpts = []; end

d = displaySet(d, 'viewing distance', vd);

%% Compute radiance, irradiance and cone absorptions
%  Create scene and compute radiance
refScene = sceneFromFile(refImg, 'rgb', [], d);
testScene = sceneFromFile(testImg, 'rgb', [], d);
fov = sceneGet(refScene, 'h fov');

%  Human optics model and compute irradiance
oi = oiCreate('human');
refOI  = oiCompute(oi, refScene);
testOI = oiCompute(oi, testScene);

%  Human cone mosaic and compute cone absorptions
sensor = sensorCreate('human', cone);
sensor = sensorSetSizeToFOV(sensor, fov, refScene, refOI);
sensor = sensorSet(sensor, 'sensor positions', zeros(nFrames, 2));

sensor = coneAbsorptions(sensor, refOI);
refP   = sensorGet(sensor, 'photons');

sensor = coneAbsorptions(sensor, testOI);
testP  = sensorGet(sensor, 'photons');

%% Classify
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
data   = cat(1, RGB2XWFormat(refP)', RGB2XWFormat(testP)');
acc = svmClassifyAcc(data, labels, nFolds, 'linear', svmOpts);

err = acc(2);
acc = acc(1);

%%