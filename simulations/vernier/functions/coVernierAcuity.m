function [jndDist, acc, err, tDist] = coVernierAcuity(params)
%% function coVernierAcuity(params)
%    Compute just noticeable difference of vernier acuity under certain
%    conditions
%
%  Inputs:
%    params: parameter structure, could include:
%      sceneFov      - scene field of view, in degree
%      barWidth      - width of the bar
%      tDist         - test misalignment vector, in degree
%      coneDensity   - human cone spatial density, default 0.6 0.3 0.1
%      expTime       - eye integration time
%      emDuration    - eye saccade duration
%      nFrames       - number of samples per scene
%      threshold     - threshold used to compute jndIntensity, default 80%
%      meanLum       - mean luminance of the scene
%      blurStd       - std of Gaussian blur of the bar
%
%  Outputs:
%    jndDist: just noticeable difference, in degree
%    acc:     classification accuray
%    err:     classification standard deviation
%    tDist:   misalignment vector tested
%
%  Notes:
%    The structure of this code is designed for proclus computation. Thus,
%    most of the intermediate output is turned off.
%
%  See also:
%    s_coVernierAcuity_Proclus
%
% (HJ) March, 2014

%% Check inputs & Init Parameters
if notDefined('params'), params = []; end

try sceneFov = params.sceneFov; catch, sceneFov = .2; end
try nFrames  = params.nFrames; catch, nFrames = 3000; end
try barWidth = params.barWidth; catch, barWidth = 5; end 
try tDist    = params.tDist; catch, tDist = ([2 6 12:4:40])/3600; end
try density  = params.coneDensity; catch, density = [0 .6 .3 .1]; end
try expTime  = params.expTime; catch, expTime = 0.05; end % integration time
try emDuration = params.emDuration; catch, emDuration = 0.01; end
try meanLum = params.meanLuminance; catch, meanLum = 100; end
try threshold = params.threshold; catch, threshold = 0.8; end

% This method is not the best trial, but I think it should be fine (HJ)
sceneSz = ceil(sceneFov/min(tDist));
if sceneSz > 10000, warning('scene resolution too high required'); end
offset = round(tDist / sceneFov * sceneSz);
tDist  = offset / sceneSz * sceneFov;

if barWidth < min(offset), warning('barwidth is too small'); end

% Init classification parameters
svmOpts = '-s 0 -q';
nFolds = 10;
labels = [ones(nFrames,1); -ones(nFrames,1)];
acc = zeros(length(tDist), 1);
err = zeros(length(tDist), 1);

%% Create Scene
%  scene{1} - scene with no misalignment
%  scene{2} - two bars with offset
params.offset = 0;
params.barWidth = barWidth;
params.sceneSz = sceneSz;

scene{1} = sceneCreate('vernier', 'object', params);
scene{1} = sceneSet(scene{1}, 'h fov', sceneFov);
scene{1} = sceneAdjustLuminance(scene{1}, meanLum);

% Show radiance image (scene)
% vcAddAndSelectObject('scene', scene{1}); sceneWindow;

%% Create Human Lens
%  Create a typical human lens
oi = oiCreate('wvf human');

% Compute optical image
% Actually, we could wait to compute it in coneSamples
% But, we compute it here to do sanity check
% vcAddAndSelectObject('scene', scene{1});
OI{1} = oiCompute(scene{1}, oi);

% Show irradiance (optical image) 
%vcAddAndSelectObject('oi', OIs{1}); oiWindow;

%% Create Sensor
% generate human photoreceptors structure
params.humanConeDensities = density;
sensor = sensorCreate('human', [], params);
sensor = sensorSetSizeToFOV(sensor, sceneFov, scene{1}, OI{1});
sensor = sensorSet(sensor, 'exp time', emDuration);

params.center   = [0,0];
%params.Sigma    = 1e-4 *[0.3280 0.0035; 0.0035 0.4873]*emDuration*1000;
params.Sigma    = 1e-4 *[0.3280 0.0; 0.0 0.0]*emDuration*1000;
%params.Sigma = zeros(2);
emPerExposure = round(expTime / emDuration);
params.nSamples = nFrames * emPerExposure;
params.fov      = sensorGet(sensor,'fov',scene{1}, OI{1});

% Set up the eye movement properties
sensor = emInit('fixation gaussian', sensor, params);

%% Generate cone samples & classify
vcAddAndSelectObject('oi', OI{1});
sensor = coneAbsorptions(sensor, OI{1}, 2);
refPhotons = RGB2XWFormat(double(sensorGet(sensor, 'photons')));

sz = size(refPhotons);
refPhotons = sum(reshape(refPhotons, ...
                [sz(1) sz(2)/emPerExposure emPerExposure]), 3);

for ii = 1 : length(tDist)
    vcDeleteSelectedObject('scene');
    vcDeleteSelectedObject('oi');
    
    params.offset = offset(ii);
    scene{2} = sceneCreate('vernier', 'object', params);
    scene{2} = sceneSet(scene{2},'h fov',sceneFov);
    scene{2} = sceneAdjustLuminance(scene{2}, meanLum);

    vcAddAndSelectObject('scene', scene{2});
    OI{2} = oiCompute(scene{2}, oi);
    vcAddAndSelectObject('oi', OI{2});
    
    sensor = coneAbsorptions(sensor, OI{2}, 2);
    matchPhotons = RGB2XWFormat(double(sensorGet(sensor, 'photons')));
    matchPhotons = sum(reshape(matchPhotons, ...
                    [sz(1) sz(2)/emPerExposure emPerExposure]), 3);
    
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'svm', svmOpts);
    
    err(ii) = accuracy(2);
    acc(ii) = accuracy(1);
end

% Find JND
try
    tRange = min(tDist):(min(tDist)/100):max(tDist);
    interpolatedAcc = interp1(tDist, acc, tRange, 'linear');
    [~, ind] = min(abs(interpolatedAcc - threshold));
    jndDist = tRange(ind);
catch
    jndDist = 0;
end


%% END