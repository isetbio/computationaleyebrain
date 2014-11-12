function [jndWave, expData] = coWaveDiscrimination(refWave, params)
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
%      .rgbDensity - human cone spatial density
%      .expTime    - human cone integration time, default 0.05 seconds
%      .nFrames    - number of samples to be used in classification
%
%  Outputs:
%    jndWave  - wavelength of light that is barely noticable to reference
%               light
%    expData  - experiment data structure
%
%
%  (HJ) ISETBIO TEAM, 2014

%% Set parameters
if notDefined('refWave'), error('reference wavelength required'); end
if notDefined('params'), params = []; end
if isscalar(refWave), warning('10nm band assumed'); end

try sceneSz = params.sceneSz; catch, sceneSz = 0.3;  end % scene fov
try expTime = params.expTime; catch, expTime = 0.05; end % integration time
try nFrames = params.nFrames; catch, nFrames = 3000; end % number of frames
try cropSz  = params.cropSz;  catch, cropSz  = 4;    end
try threshold = params.threshold; catch, threshold = 0.8; end

svmOpts = '-s 0 -q';
nFolds = 10;

%% Create reference scene
%  scene{1} - uniform patch with reference wavelength
scene{1} = sceneCreate('uniform monochromatic', refWave, 128);
scene{1} = sceneSet(scene{1}, 'fov', sceneSz);

% by default, scene{1} is equal luminance with mean luminance equals 100
% we can adjust it to have equal energy, say 1
%
% energy = sum(sceneGet(scene{1}, 'energy'), 3);
% p = sceneGet(scene{1}, 'photons');
% scene{1} = sceneSet(scene{1}, 'photons', p / mean(energy(:)));

%% Create Human Optics
%  create lens for standard human observer
%  Optics for both scenes are the same. But we use different oi structure
%  to store the computed optical image

% generate human optics structure
oi = oiCreate('human');

% compute optical image
OI{1} = oiCompute(scene{1}, oi);

%% Create Human Photoreceptors (cones)
%  create standard human photoreceptors
%  the absorption of the photoreceptors is adjusted by macular pigment and
%  lens density
%
% generate human photoreceptors structure
sensor = sensorCreate('human', [], params);
sensor = sensorSet(sensor, 'exposure time', expTime);
sensor = sensorSetSizeToFOV(sensor, sceneSz, scene{1}, OI{1});
sensor = sensorSet(sensor, 'sensorpositions', zeros(nFrames, 2));

%% Compute Human cone absorptions & Classification
%  Compute samples of cone absorptions
%  The absorptions for the two groups are normalized to have same mean
%
sensor = coneAbsorptions(sensor, OI{1}, 0);
refPhotons = sensorGet(sensor, 'photons');

% do L-M wiring and add second site noise
coneType = sensorGet(sensor, 'cone type');
refPhotons = coneComputeSSNoise(refPhotons, coneType);

% crop out a small area and transform to XW format
refPhotons = RGB2XWFormat(getMiddleMatrix(refPhotons, cropSz));

% Init labels and expData
labels = [ones(nFrames,1); -ones(nFrames,1)];
expData.acc = []; expData.err = []; expData.tWave = [];
w_min = 0; acc_min = 0.5; w_max = 16; acc_max = 1;

% find an upper bound for w_max
while acc_max < threshold
    tWave = refWave + w_max;
    scene{2} = sceneCreate('uniform monochromatic', tWave, 128);
    scene{2} = sceneSet(scene{2}, 'fov', sceneSz);
    
    % adjust to equal energy
    % energy = sum(sceneGet(scene{2}, 'energy'), 3);
    % p = sceneGet(scene{2}, 'photons');
    % scene{2} = sceneSet(scene{2}, 'photons', p / mean(energy(:)));
    
    OI{2}  = oiCompute(scene{2}, oi);
    sensor = coneAbsorptions(sensor, OI{2}, 0);
    matchPhotons = sensorGet(sensor, 'photons');
    
    % L-M wiring and second site noise
    matchPhotons = coneComputeSSNoise(matchPhotons, coneType);
    
    % croping and transform to XW format
    matchPhotons = RGB2XWFormat(getMiddleMatrix(matchPhotons, cropSz));
    
    % do classification
    accuracy = svmClassifyAcc(cat(1, refPhotons', matchPhotons'), ...
        labels, nFolds, 'linear', svmOpts);
    
    acc_max = accuracy(1);
    w_max = w_max + 4;
    
    assert(w_max < 48, 'Cannot classify with very large wave difference');
end

% Now if the test wave values are not given, we should find a proper range
% for it by ourselves
% The general process for doing this is as below:
%   1) find a wave w_min with accuracy between 70~threshold
%   2) find a wave w_max with accuracy between threshold~90
%   3) linear interpolate between w_min and w_max
while true
    tWave = refWave + (w_min + w_max) / 2;
    scene{2} = sceneCreate('uniform monochromatic', tWave, 128);
    scene{2} = sceneSet(scene{2}, 'fov', sceneSz);
    
    % adjust to equal energy
    % energy = sum(sceneGet(scene{2}, 'energy'), 3);
    % p = sceneGet(scene{2}, 'photons');
    % scene{2} = sceneSet(scene{2}, 'photons', p / mean(energy(:)));
    
    OI{2}  = oiCompute(scene{2}, oi);
    sensor = coneAbsorptions(sensor, OI{2}, 0);
    matchPhotons = sensorGet(sensor, 'photons');
    
    % L-M wiring and second site noise
    matchPhotons = coneComputeSSNoise(matchPhotons, coneType);
    
    % croping and transform to XW format
    matchPhotons = RGB2XWFormat(getMiddleMatrix(matchPhotons, cropSz));
    
    % do classification
    accuracy = svmClassifyAcc(cat(1, refPhotons', matchPhotons'), ...
        labels, nFolds, 'linear', svmOpts);
    fprintf('wave: %f\t acc:%f\n', (w_min + w_max)/2, accuracy(1));
    
    % log experiment data
    expData.acc   = cat(1, expData.acc, accuracy(1));
    expData.err   = cat(1, expData.err, accuracy(2));
    expData.tWave = cat(1, expData.tWave, tWave(1));
    
    % update test condition
    if accuracy(1) < threshold
        w_min = (w_min + w_max)/2;
        acc_min = accuracy(1);
    else
        w_max = (w_min + w_max)/2;
        acc_max = accuracy(1);
    end
    
    % check whether or not to continue
    if acc_min > threshold - 0.05 && acc_max < threshold + 0.05
        break;
    end
    
    if acc_max < threshold + 0.01 || acc_min > threshold - 0.01
        jndWave = tWave(1) - refWave(1);
        return;
    end
end

%% Identify JND wavelength by interpolation
jndWave = w_min + (w_max-w_min) * (threshold-acc_min) / (acc_max-acc_min);


end
%% END