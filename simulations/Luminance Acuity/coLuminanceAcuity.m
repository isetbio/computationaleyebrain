function [jndIntensity, acc, err, tIntensity] = coLuminanceAcuity(intensity, params)
%% function coLuminanceAcuity
%    This function computes the just noiticable different for certain
%    intensity level by using the computational observer
%
%  Inputs:
%    intensity - reference intensity in cd/m2
%    params    - parameter structure, could include
%      sceneSz       - scene field of view (degree)
%      coneDensity   - human cone spatial density, default 0.6 0.3 0.1
%      expTime       - eye integration time
%      emDuration    - eye saccade duration
%      nFrames       - number of samples per scene
%      threshold     - threshold used to compute jndIntensity, default 80%
%      testIntensity - array of test intensity to be tested
%
%  Outputs:
%      jndIntensity  - JND of intensity
%      acc           - classification accuracy
%      err           - std of classification accuracy
%      tIntensity    - intensities tested
%
%  Notes:
%    The structure of this code is designed for proclus computation. Thus,
%    most of the intermediate output is turned off.
%
%  See also:
%    s_luminanceAcuity_Proclus
%
%  (HJ) March 2014

%% Check inputs and set parameters
%  check inputs
if notDefined('intensity'), error('reference intensity required'); end
if notDefined('params'), params = []; end

%  set parameters
try sceneSz = params.sceneSz; catch, sceneSz = 0.1; end % fov (degree)
try density = params.coneDensity; catch, density = [0 .6 .3 .1]; end
try expTime = params.expTime; catch, expTime = 0.05; end % integration time
try emDuration = params.emDuration; catch, emDuration = 0.01; end
try nFrames = params.nFrames; catch, nFrames = 3000; end

if isfield(params, 'tIntensity')
    tIntensity = params.tIntensity;
else
    tIntensity = intensity + [.01 .02 .04 .08 .15 .3 .45 0.6 0.8 1 2 4];
end

svmOpts = '-s 0 -q';
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
acc = zeros(length(tIntensity), 1);
err = zeros(length(tIntensity), 1);

%% Create scene
scene = cell(2,1);
scene{1} = sceneCreate('uniform');
scene{1} = sceneSet(scene{1}, 'fov', sceneSz);
scene{1} = sceneAdjustLuminance(scene{1}, intensity);

%% Create OI
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

emPerExposure = round(expTime / emDuration);
params.nSamples = nFrames * emPerExposure;
params.fov      = sensorGet(sensor,'fov',scene{1}, OI{1});

% Set up the eye movement properties
sensor = emInit('fixation gaussian', sensor, params);

%% Compute cone absorptions and do classification
sensorSz = sensorGet(sensor, 'size');
sensorC = round(sensorSz/2);
vcAddAndSelectObject('oi', OI{1});

sensor = coneAbsorptions(sensor, OI{1}, 2);
refPhotons = double(sensorGet(sensor, 'photons'));

% take a small portion out
refPhotons = refPhotons(sensorC(1)-3 : sensorC(1)+3, ...
    sensorC(2)-3 : sensorC(2)+3, :);
refPhotons   = RGB2XWFormat(refPhotons);
sz = size(refPhotons);
refPhotons = sum(reshape(refPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);


for ii = 1 : length(tIntensity)
    vcDeleteSelectedObject('scene');
    scene{2} = sceneAdjustLuminance(scene{1}, tIntensity(ii));

    vcAddAndSelectObject('scene', scene{2});
    OI{2} = oiCompute(scene{2}, oi);
    
    sensor = coneAbsorptions(sensor, OI{2}, 2);
    matchPhotons = double(sensorGet(sensor, 'photons'));
    matchPhotons = matchPhotons(sensorC(1)-3 : sensorC(1)+3, ...
        sensorC(2)-3 : sensorC(2)+3, :);
    
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
    jndIntensity = interp1(acc(ind), tIntensity(ind), threshold, 'linear');
catch
    tRange = min(tIntensity):(min(tIntensity-intensity)/100):max(tIntensity);
    interpolatedAcc = interp1(tIntensity, acc, tRange, 'linear');
    [~, ind] = min(abs(interpolatedAcc - threshold));
    jndIntensity = tRange(ind);
end
%% END