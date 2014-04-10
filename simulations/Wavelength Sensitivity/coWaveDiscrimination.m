function [jndWave, acc, err, wave] = coWaveDiscrimination(refWave, params)
%% function coWaveDiscrimination(refWave, [params])
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
if notDefined('refWave'), error('reference wavelength required'); end
if notDefined('params'), params = []; end
if isscalar(refWave), warning('10nm band assumed'); end

try sceneSz = params.sceneSz; catch, sceneSz = 0.1;end % scene fov
try density = params.density; catch, density = [0 .6 .3 .1]; end  % KLMS
try expTime = params.expTime; catch, expTime = 0.05; end % integration time
try nFrames = params.nFrames; catch, nFrames = 3000; end
try direction = params.direction; catch, direction = 'up'; end
try threshold = params.threshold; catch, threshold = 0.8; end

if isfield(params, 'tWave')
    tWave = params.tWave;
else
    tWave = [];
end

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
scene{1} = sceneAdjustLuminance(scene{1}, 100);
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

sensorSz = sensorGet(sensor, 'size');
sensorC = round(sensorSz/2);


%% Compute Human cone absorptions & Classification
%  Compute samples of cone absorptions
%  The absorptions for the two groups are normalized to have same mean
%

sensor = coneAbsorptions(sensor, OI{1}, 0);
refPhotons = double(sensorGet(sensor, 'photons'));
refPhotons = refPhotons(sensorC(1)-2 : sensorC(1)+2, ...
    sensorC(2)-2 : sensorC(2)+2, :);
refPhotons = RGB2XWFormat(refPhotons);

labels = [ones(nFrames,1); -1*ones(nFrames,1)];
acc = zeros(length(tWave), 1);
err = zeros(length(tWave), 1);

% Now if the test wave values are not given, we should find a proper range
% for it by ourselves
% The general process for doing this is as below:
%   1) find a wave w_min with accuracy between 70~threshold
%   2) find a wave w_max with accuracy between threshold~90
%   3) linear interpolate between w_min and w_max
if isempty(tWave)
    w_min = 0; acc_min = 0.5;
    w_max = 16; acc_max = 1;
    iter = 1;
    while true
        wave = refWave + (w_min + w_max) / 2;
        vcDeleteSelectedObject('scene');
        scene{2} = sceneCreate('uniform monochromatic', wave, 128);
        scene{2} = sceneSet(scene{2}, 'fov', sceneSz);
        vcAddAndSelectObject('scene', scene{2});
        vcDeleteSelectedObject('oi');
        OI{2} = oiCompute(scene{2}, oi);
        vcAddAndSelectObject('oi', OI{2});
        sensor = coneAbsorptions(sensor, OI{2}, 0);
        matchPhotons = double(sensorGet(sensor, 'photons'));
        matchPhotons = matchPhotons(sensorC(1)-2 : sensorC(1)+2, ...
            sensorC(2)-2 : sensorC(2)+2, :);
        matchPhotons = RGB2XWFormat(matchPhotons);
        
        matchPhotons = matchPhotons * mean(refPhotons(:))/ mean(matchPhotons(:));
        % accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        %     labels, nFolds, 'svm', svmOpts);
        accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
            labels, nFolds, 'svm', svmOpts);
        acc = accuracy(1);
        fprintf('wave: %f\t acc:%f\n', (w_min + w_max)/2, acc);
        if acc < threshold
            w_min = (w_min + w_max)/2;
            acc_min = acc;
        else
            w_max = (w_min + w_max)/2;
            acc_max = acc;
        end
        
        % Make sure it will not loop infinitely
        iter = iter + 1;
        if iter > 8
            tWave = [0.1 0.2 0.3 0.4 0.6 0.8 1 1.3 1.5 1.7 2:0.5:10];
            break;
        end
        
        if acc_min > 0.7 && acc_min < threshold && ...
                acc_max > threshold && acc_max < 0.9
            tWave = (w_min-0.5):0.1:(w_max+0.5);
            break;
        end
    end
end

[X,Y] = meshgrid(refWave, tWave);
switch direction
    case 'up'
        wave = X + Y;
    case 'down'
        wave = X - Y;
    otherwise
        error('Unknown direction');
end

for ii = 1 : length(tWave)
    vcDeleteSelectedObject('scene');
    scene{2} = sceneCreate('uniform monochromatic', wave(ii,:), 128);
    scene{2} = sceneSet(scene{2}, 'fov', sceneSz);
    vcAddAndSelectObject('scene', scene{2});
    vcDeleteSelectedObject('oi');
    OI{2} = oiCompute(scene{2}, oi);
    vcAddAndSelectObject('oi', OI{2});
    sensor = coneAbsorptions(sensor, OI{2}, 0);
    matchPhotons = double(sensorGet(sensor, 'photons'));
    matchPhotons = matchPhotons(sensorC(1)-2 : sensorC(1)+2, ...
        sensorC(2)-2 : sensorC(2)+2, :);
    matchPhotons = RGB2XWFormat(matchPhotons);
    
    matchPhotons = matchPhotons * mean(refPhotons(:))/ mean(matchPhotons(:));
    % accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
    %     labels, nFolds, 'svm', svmOpts);
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'svm', svmOpts);
    
    err(ii) = accuracy(2);
    acc(ii) = accuracy(1);
end

%% Identify JND wavelength
wave = wave(:,1);
try
    [~, ind] = sort(acc);
    jndWave = interp1(acc(ind), wave(ind), threshold, 'linear');
catch
    tRange = min(wave):(min(wave)/100):max(wave);
    interpolatedAcc = interp1(wave, acc, tRange, 'linear');
    [~, ind] = min(abs(interpolatedAcc - threshold));
    jndWave = tRange(ind);
end

end
%% END