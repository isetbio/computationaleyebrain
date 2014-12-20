%% s_computationalEquivalence
%    This script demonstrates the information loss in different stages of
%    human vision system.
%    Two cases are shown here:
%    1) Different radiance that has same / similar optical irradiance
%       (high spatial frequency harmonincs)
%    2) Different optical image that has same / similar cone photon
%       absoptions (vernier acuity)
%
%  (HJ) Dec, 2014

%% Init
s_initISET;

%% Create harmonic scene
%  Set parameters
%  the final harmonics are computed as
%       contrast*window.*cos(2*pi*f*([cos(ang)*X + sin(ang)*Y] + ph)) + 1
params.freq = 120;       % spatial frequency
params.contrast = 1;     % full contrast
params.GarborFlag = 0.1; % Garbor filter
params.ph = 0;           % phase
params.ang = 0;          % angle
params.row = 240;        % number of samples in vertical direction
params.col = 240;        % number of samples in horizontal direction

%  create scene
sceneH = sceneCreate('harmonic', params);

%  adjust illuminant to D65
il = illuminantCreate('D65', sceneGet(sceneH, 'wave'));
sceneH = sceneAdjustIlluminant(sceneH, il);

%  visualize
vcAddObject(sceneH); sceneWindow;

%% Create uniform scene
%  set parameters
sceneSz = [params.row params.col];

%  create uniform scene with D65 illuminance
sceneU = sceneCreate('uniform D65', sceneSz);
sceneU = sceneSet(sceneU, 'h fov', sceneGet(sceneH, 'h fov')); % set fov
sceneU = sceneAdjustLuminance(sceneU, sceneGet(sceneH, 'mean luminance'));

%  visualize
vcAddObject(sceneU); sceneWindow;

%% Effect of human optics
%  create optics for standard human observer
oi = oiCreate('human');

%  compute irradiance
oiH = oiCompute(sceneH, oi); % opitcal image for harmonic scene
oiU = oiCompute(sceneU, oi); % optical image for uniform scene

%  visualize
vcAddObject(oiH); vcAddObject(oiU);
oiWindow;

%% Compare harmonic and uniform
%  Plot scene radiance
vcNewGraphWin;

subplot(2,2,1); imshow(sceneGet(sceneH, 'rgb image'));
 title('Scene Radiance (Harmonics)');

subplot(2,2,3); imshow(sceneGet(sceneU, 'rgb image'));
title('Scene Radiance (Uniform)');

% Plot optical irradiance map
subplot(2,2,2); imshow(oiGet(oiH, 'rgb image'));
title('Optical Image (Harmonics)');

subplot(2,2,4); imshow(oiGet(oiU, 'rgb image'));
title('Optical Image (Uniform)');

%% Create scenes for aligned and mis-aligned lines
%  set parameters
params = []; % clean up
params.sceneSz   = [240 240]; % scene size (number of samples)
params.barWidth  = 2;         % bar width in number of samples
params.barReflect = 1;        % bar reflectance
params.bgReflect  = 0;        % background reflectance
params.il         = il;       % illuminance

%  create scene - misaligned lines
params.offset = 1;         % misaligned offset in number of samples
sceneV = sceneCreate('vernier', 'object', params); % misaligned
sceneV = sceneSet(sceneV, 'h fov', 0.5);

%  create scene - aligned lines
params.offset = 0;
sceneA = sceneCreate('vernier', 'object', params); % aligned
sceneA = sceneSet(sceneA, 'h fov', 0.5);

%  visualize
vcAddObject(sceneV); vcAddObject(sceneA);
sceneWindow;

%% Optical irradiance image
%  Compute irradiance map
oiV = oiCompute(sceneV, oi);
oiA = oiCompute(sceneA, oi);

%  Visualize
vcAddObject(oiV); vcAddObject(oiA);
oiWindow;

%% Compute cone absorptions
%  Create cone mosaics for standard human observer
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, oiGet(oi, 'fov'), [], oi);

expTime  = sensorGet(sensor, 'exp time');
sampTime = sensorGet(sensor, 'time interval');
nSamples = round(expTime / sampTime);

%  Set eye-movement params
sensor = sensorSet(sensor, 'exp time', sampTime); 
sensor = sensorSet(sensor, 'em type', [1 0 1]); % tremor and micro-saccade

params = []; % clean up
params.nSamples = nSamples;

sensor = eyemoveInit(sensor, params);

%  Compute photon absorptions
sensorV = sensorCompute(sensor, oiV);
p = sensorGet(sensorV, 'photons'); % photon absorption for each 1 ms
p = sum(p, 3); % photon absorption during the exposure time
sensorV = sensorSet(sensorV, 'exp time', expTime);
sensorV = sensorSet(sensorV, 'photons', p);

sensorA = sensorCompute(sensor, oiA);
p = sensorGet(sensorA, 'photons'); % photon absorption for each 1 ms
p = sum(p, 3); % photon absorption during the exposure time
sensorA = sensorSet(sensorA, 'exp time', expTime);
sensorA = sensorSet(sensorA, 'photons', p);

%  Visualize
vcAddObject(sensorV); vcAddObject(sensorA);
sensorWindow('scale', true);

%% Compare aligned and mis-aligned lines
%  Plot optical irradiance image
vcNewGraphWin;
subplot(2,2,1); imshow(oiGet(oiV, 'rgb image'));
title('Optical Image (Mis-aligned)');

subplot(2,2,3); imshow(oiGet(oiA, 'rgb image'));
title('Optical Image (Aligned)');

%  Plot cone absorptions
gam = 1; % gamma
scale = true; % scale image


subplot(2,2,2); 
imshow(sensorGet(sensorV, 'rgb image', 'volts', gam, scale));
title('Cone Photon Absorptions (Misaligned)');

subplot(2,2,4); 
imshow(sensorGet(sensorA, 'rgb image', 'volts', gam, scale));
title('Cone Photon Absorption (Aligned)');
