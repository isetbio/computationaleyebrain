%% s_contrastSensitivity
%
%  Calculate the contrast sensitivity functions of the computational
%  observer based on the cone absorptions.
%
%
%  (HJ) March, 2014

%% Set parameters
if notDefined('dpi'), ppi = 200; end % 100 pixel per inch
if notDefined('ppc'), ppc = 7;   end % pixels per half cycle
if notDefined('sceneSz'), sceneSz = 0.5; end % scene field of view
if notDefined('refColor'), refColor = 0.5; end
if notDefined('testColor'), testColor = 0.48; end
if notDefined('viewingDst'), viewingDst = 1; end % 1 meter
if notDefined('density'), density = [0 .6 .3 .1]; end % cone density
if notDefined('expTime'), expTime = 0.05; end % eye integration time
if notDefined('emDuration'), emDuration = 0.01; end % eye saccade duration
if notDefined('nFrames'), nFrames = 1500; end

%% Create display model
display = displayCreate('LCD-Apple');
display = displaySet(display, 'gTable', 'linear');
display = displaySet(display, 'dpi', ppi);

%% Create scene
%  scene{1} - patch with spatially varying patterns
%  scene{2} - uniform patch
%

% compute patch size in number of pixels
patchSz = 2 * tand(sceneSz/2) * viewingDst * 39.37 * ppi;
patchSz = round(patchSz/2/ppc)*2*ppc;
sceneSz = atand(patchSz/2/viewingDst/39.37/ppi)*2;

% create image to be shown on virtual display
% image = ones(patchSz) * refColor;

% create image with spatially varying patterns
nPeriod = floor(patchSz/2/ppc);
assert(nPeriod > 0, 'Frequency is too low.');
image = repmat(refColor + (refColor-testColor)*cos((1:patchSz)/ppc*pi) ,[patchSz, 1]);
%for ii = 1 : nPeriod
%    image(:, (ii-1)*ppc*2+1: ii*ppc*2 - ppc) = testColor;
%end

%%
params.freq = 40;   % Per FOV
params.contrast = 1;
params.ph = pi;
params.ang = 0;
params.row = 192; parsams.col = 192;
params.GaborFlag = .2;
clear scene

scene{1} = sceneCreate('harmonic',params);
scene{1} = sceneSet(scene{1},'fov',1);
scene{1} = sceneAdjustLuminance(scene{1},0.1);
vcAddObject(scene{1}); sceneWindow;

params.contrast = 0;
scene{2} = sceneCreate('harmonic',params);
scene{2} = sceneSet(scene{2},'fov',1);
scene{2} = sceneAdjustLuminance(scene{2},0.1);

% % create scene for patterned image
% scene{1} = sceneFromFile(image, 'rgb', [], display);
% scene{1} = sceneSet(scene{1}, 'distance', viewingDst);
% scene{1} = sceneSet(scene{1}, 'h fov', sceneSz);
% 
% % create uniform patch scene
% image(:) = mean(image(:));
% scene{2} = sceneFromFile(image, 'rgb', [], display);
% scene{2} = sceneSet(scene{2}, 'distance', viewingDst);
% scene{2} = sceneSet(scene{2}, 'h fov', sceneSz);

%% Create Human Optics
%  create lens for standard human observer
%  Optics for both scenes are the same. But we use different oi structure
%  to store the computed optical image
%

% generate human optics structure
wave   = 380 : 4 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

% compute optical image
vcAddAndSelectObject('scene', scene{1});
OI{1} = oiCompute(scene{1}, oi);
vcAddAndSelectObject('scene', scene{2});
OI{2} = oiCompute(scene{2}, oi);

%% Create Human Photoreceptors (cones)
%  create standard human photoreceptors  
%  the absorption of the photoreceptors is adjusted by macular pigment and
%  lens density
%

% generate human photoreceptors structure
if exist('density', 'var')
    pparams.humanConeDensities = density;
else
    pparams = [];
end

sensor = sensorCreate('human', [], pparams);
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
absorptions = cell(2, 1);
sensorSz = sensorGet(sensor, 'size');
sensorC = round(sensorSz/2);
vcAddAndSelectObject('oi', OI{1});
for ii = 1 : 2
    sensor = coneAbsorptions(sensor, OI{ii}, 2);
    absorptions{ii} = double(sensorGet(sensor, 'photons'));
    % take a small portion out
    % this should not be hard coded, change it (HJ)
    absorptions{ii} = absorptions{ii}(sensorC(1)-2:sensorC(1)+2, ...
            round(sensorC(2)/2):round(sensorC(2)*3/2), :);
end

absorptions{1} = absorptions{1} * sum(absorptions{2}(:)) / sum(absorptions{1}(:));

%% SVM Classification
%  Do classification
%
% svmOpts = '-s 0 -q';
svmOpts = '-q';

nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
refPhotons   = RGB2XWFormat(absorptions{1});
%refPhotons = cumsum(refPhotons, 2);
%refPhotons = refPhotons(:, emPerExposure + 1:end) - ...
%    refPhotons(:, 1:end-emPerExposure);
sz = size(refPhotons);
refPhotons = sum(reshape(refPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);

matchPhotons = RGB2XWFormat(absorptions{2});
%matchPhotons = cumsum(matchPhotons, 2);
%matchPhotons = matchPhotons(:, emPerExposure + 1:end) - ...
%    matchPhotons(:, 1:end-emPerExposure);
matchPhotons = sum(reshape(matchPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);

% accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
%     labels, nFolds, 'svm', svmOpts);
accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
    labels, nFolds, 'linear', svmOpts);

err = accuracy(2);
acc = accuracy(1);

%% END