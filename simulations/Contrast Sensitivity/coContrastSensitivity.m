function [jndContrast, acc, err, tContrast] = coContrastSensitivity(frequency, params)
%% function coContrastSensitivity(contrast, frequency, [params])
%    Calculate the contrast sensitivity functions of the computational
%    observer based on the cone absorptions.
%
%  Inputs:
%    frequency  - spatial frequency of the harmonics (cpd)
%    params     - parameter structure, can include
%      sceneSz       - scene field of view (degree)
%      coneDensity   - human cone spatial density, default 0.6 0.3 0.1
%      expTime       - eye integration time
%      emDuration    - eye saccade duration
%      nFrames       - number of samples per scene
%      threshold     - threshold used to compute jndContrast, default 80%
%      testContrast  - array of contrast to be tested
%      meanLuminance - mean luminance of the scenes in cd/m2
%
%  Outputs:
%    jndContrast - JND for given threshold
%    acc         - classification accuracy
%    err         - std of classificaiton error
%    tContrast   - array of contrast tested
%
%  (HJ) March, 2014

%% Check inputs and set parameters
%  check inputs
if notDefined('frequency'), error('freqency required'); end
if notDefined('params'), params = []; end

%  set parameters
try sceneSz = params.sceneSz; catch, sceneSz = 0.5; end % fov (degree)
try density = params.coneDensity; catch, density = [0 .6 .3 .1]; end
try expTime = params.expTime; catch, expTime = 0.05; end % integration time
try emDuration = params.emDuration; catch, emDuration = 0.01; end
try nFrames = params.nFrames; catch, nFrames = 3000; end
try meanLum = params.meanLuminance; catch, meanLum = 100; end
try sceneWave = params.sceneWave; catch, sceneWave = 380:780; end
try threshold = params.threshold; catch, threshold = 0.8; end

if isfield(params, 'testContrast')
    tContrast = params.testContrast;
else
    tContrast = [.1 .05 .04 .03 .02 .01 .005 0.004 0.003 .002];
end

svmOpts = '-s 0 -q';
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
acc = zeros(length(tContrast), 1);
err = zeros(length(tContrast), 1);

%% Create scene
%  scene{1} - uniform patch
%  scene{2} - patch with spatially varying patterns

pparams.freq = frequency;    % spatial frequency
pparams.ph = pi; % phase, 0 - black in center, pi - white in center
pparams.ang = 0; % direction, 0 - vertical
pparams.row = 192; pparams.col = 192; % sample size 
pparams.GaborFlag = .2; % Gabor blur

pparams.contrast = 0;
scene{1} = sceneCreate('harmonic',pparams, sceneWave);
scene{1} = sceneSet(scene{1},'fov',1);
scene{1} = sceneAdjustLuminance(scene{1}, meanLum);

%% Create Human Optics
%  create lens for standard human observer
%  Optics for both scenes are the same. But we use different oi structure
%  to store the computed optical image
%

% generate human optics structure
wave   = 380 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

% compute optical image
vcAddAndSelectObject('scene', scene{1});
OI{1} = oiCompute(scene{1}, oi);

%% Create Human Photoreceptors (cones)
%  create standard human photoreceptors  
%  the absorption of the photoreceptors is adjusted by macular pigment and
%  lens density
%

% generate human photoreceptors structure
params.humanConeDensities = density;
sensor = sensorCreate('human', [], params);
sensor = sensorSetSizeToFOV(sensor, sceneSz, scene{1}, OI{1});
sensor = sensorSet(sensor, 'exp time', emDuration);

params.center   = [0,0];
params.Sigma    = 1e-4 *[0.3280 0.0035; 0.0035 0.4873]*emDuration*1000;
%params.Sigma = zeros(2);
emPerExposure = round(expTime / emDuration);
params.nSamples = nFrames * emPerExposure;
params.fov      = sensorGet(sensor,'fov',scene{1}, OI{1});

% Set up the eye movement properties
sensor = emInit('fixation gaussian', sensor, params);


%% Compute Human cone absorptions
%  Compute samples of cone absorptions
%  The absorptions for the two groups are normalized to have same mean
%
sensorSz = sensorGet(sensor, 'size');
sensorC = round(sensorSz/2);
vcAddAndSelectObject('oi', OI{1});

sensor = coneAbsorptions(sensor, OI{1}, 2);
refPhotons = double(sensorGet(sensor, 'photons'));

% take a small portion out
refPhotons = refPhotons(round(sensorC(1)/2):round(sensorC(1)*3/2), ...
    round(sensorC(2)/2):round(sensorC(2)*3/2), :);
refPhotons   = RGB2XWFormat(refPhotons);
sz = size(refPhotons);
refPhotons = sum(reshape(refPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);


for ii = 1 : length(tContrast)
    vcDeleteSelectedObject('scene');
    pparams.contrast = tContrast(ii); % contrast
    scene{2} = sceneCreate('harmonic',pparams, sceneWave);
    scene{2} = sceneSet(scene{2},'fov',1);
    scene{2} = sceneAdjustLuminance(scene{2}, meanLum);

    vcAddAndSelectObject('scene', scene{2});
    OI{2} = oiCompute(scene{2}, oi);
    
    sensor = coneAbsorptions(sensor, OI{2}, 2);
    matchPhotons = double(sensorGet(sensor, 'photons'));
    matchPhotons = matchPhotons(round(sensorC(1)/2):round(sensorC(1)*3/2), ...
        round(sensorC(2)/2):round(sensorC(2)*3/2), :);
    
    matchPhotons = RGB2XWFormat(matchPhotons);
    matchPhotons = sum(reshape(matchPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);
    
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'svm', svmOpts);
    
    err(ii) = accuracy(2);
    acc(ii) = accuracy(1);
end

% Find JND
try
    [~, ind] = sort(acc);
    jndContrast = interp1(acc(ind), tContrast(ind), threshold, 'linear');
catch
    tRange = min(tContrast):(min(tContrast)/100):max(tContrast);
    interpolatedAcc = interp1(tContrast, acc, tRange, 'linear');
    [~, ind] = min(abs(interpolatedAcc - threshold));
    jndContrast = tRange(ind);
end
%% END