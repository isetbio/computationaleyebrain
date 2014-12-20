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
subplot(2,2,2);imshow(oiGet(oiH, 'rgb image'));
title('Optical Image (Harmonics)');

subplot(2,2,4); imshow(oiGet(oiU, 'rgb image'));
title('Optical Image (Uniform)');

%% Create scenes for aligned and mis-aligned lines

%% Optical irradiance image

%% Compute cone absorptions

%% Compare aligned and mis-aligned lines