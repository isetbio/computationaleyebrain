function [jndWave, acc, err, wave] = wavelengthAcuity(refWave, params)
%% function wavelengthAcuity(wave, [params])
%    This function calculates the acuity for given wavelength. The acuity
%    is given by wavelength distance of two just noiticable difference
%    lights
%
%  Inputs:
%    refWave  - reference wavelength
%    params   - parameter structure, could contain
%      .threshold  - percentage of correctness of JND, default 0.8
%      .sceneSz    - scene fov in degrees, default 0.1
%      .direction  - searching direction, could be 'up', 'down'
%                    for 'up', it finds jndWave that is larger than wave
%      .density    - human cone spatial density, default [0 .6 .3 .1]
%      .expTime    - human cone integration time, default 0.05 seconds
%      .nFrames    - number of samples to be used in classification
%      .wavelength - wavelength of harmonic scene
%
%  Outputs:
%    jndWave  - wavelength of light that is barely noticable to reference
%               light
%
%  Notes:
%    The computation is done without cone opponency or any rgc features
%
%  (HJ) March, 2014

%% Set parameters
if notDefined('params'), params = []; end
assert(isscalar(refWave), 'reference wave should be scalar');

try sceneSz = params.sceneSz; catch, sceneSz = 0.05;end % scene fov
try density = params.density; catch, density = [0 .6 .3 .1]; end  % KLMS
try expTime = params.expTime; catch, expTime = 0.05; end % integration time
try nFrames = params.nFrames; catch, nFrames = 3000; end
try direction = params.direction; catch, direction = 'up'; end
try threshold = params.threshold; catch, threshold = 0.8; end

wave   = 380 : 780;
pupilDiameterMm = 3;
svmOpts = '-s 0 -q';
nFolds = 10;

%% Create reference scene
%  scene{1} - uniform patch with reference wavelength
scene{1} = sceneCreate('uniform monochromatic', refWave, 128);
scene{1} = sceneSet(scene{1}, 'fov', sceneSz);

% XYZ = ieReadSpectra('XYZ', 380:780);
% adjLum = XYZ(wave==refWave,2) / XYZ(wave==550, 2) * 100;
% scene{1} = sceneAdjustLuminance(scene{1}, adjLum);
scene{1} = sceneAdjustLuminance(scene{1}, 1);
vcAddAndSelectObject('scene', scene{1});
% sceneWindow;

%% Create Human Optics
%  create lens for standard human observer
%  Optics for both scenes are the same. But we use different oi structure
%  to store the computed optical image

% generate human optics structure
wvf    = wvfCreate('wave',wave);
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

% compute optical image
OI{1} = oiCompute(scene{1}, oi);
vcAddAndSelectObject('oi', OI{1});

%% Create Human Photoreceptors (cones)
%  create standard human photoreceptors  
%  the absorption of the photoreceptors is adjusted by macular pigment and
%  lens density
%
% generate human photoreceptors structure
pparams.humanConeDensities = density;
sensor = sensorCreate('human', [], pparams);
sensor = sensorSet(sensor, 'exp time', expTime);
sensor = sensorSetSizeToFOV(sensor, sceneSz, scene{1}, OI{1});

sensor = sensorSet(sensor,'movement positions', [0 0]);
sensor = sensorSet(sensor,'frames per position', nFrames);



%% Compute Human cone absorptions & Classification
%  Compute samples of cone absorptions
%  The absorptions for the two groups are normalized to have same mean
%
sampleWave = [1 2 4 8 13 18 24 28 32];
switch direction
    case 'up'
        wave = refWave + sampleWave;
    case 'down'
        wave = refWave - sampleWave;
    otherwise
        error('Unknown direction');
end

sensor = coneAbsorptions(sensor, OI{1}, 0);
refPhotons = RGB2XWFormat(double(sensorGet(sensor, 'photons')));

labels = [ones(nFrames,1); -1*ones(nFrames,1)];
acc = zeros(length(sampleWave), 1);
err = zeros(length(sampleWave), 1);

for ii = 1 : length(sampleWave)
    vcDeleteSelectedObject('scene');
    scene{2} = sceneCreate('uniform monochromatic', wave(ii), 128);
    scene{2} = sceneSet(scene{2}, 'fov', sceneSz);
    vcAddAndSelectObject('scene', scene{2});
    vcDeleteSelectedObject('oi');
    OI{2} = oiCompute(scene{2}, oi);
    vcAddAndSelectObject('oi', OI{2});
    sensor = coneAbsorptions(sensor, OI{2}, 0);
    matchPhotons = RGB2XWFormat(double(sensorGet(sensor, 'photons')));
    
    % accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
    %     labels, nFolds, 'svm', svmOpts);
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'svm', svmOpts);
    
    err(ii) = accuracy(2);
    acc(ii) = accuracy(1);
end

%% Identify JND wavelength
jndWave = 0;
try
    [~, ind] = sort(acc);
    jndWave = interp1(acc(ind), wave(ind), threshold, 'linear');
catch
    warning('JND wave cannot be found');
end

end
%% END