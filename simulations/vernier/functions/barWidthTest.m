%% s_barWidthTest
%
%  This is just a test scritp. Nothing special.
%
% (HJ) April, 2014

%% Init
wave = 400:10:700;
sceneSz = 1200;
fov = 1;
vDist = 1;
pupilDiameterMm = 3; 

%% Create scenes
param.sceneSz = sceneSz;
param.offset = 0;
param.barWidth = 1;
param.barReflect = 1;
param.bgReflect = 0.1;

scene = sceneCreate('vernier', 'object', param);
scene = sceneSet(scene, 'fov', fov);
scene = sceneSet(scene, 'distance', vDist);

vcAddAndSelectObject('scene', scene);
sceneWindow;

%% Create human optics and compute optical image
wvf    = wvfCreate('wave',wave);
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

% compute oi
oi = oiCompute(scene, oi);

% visualize
vcAddAndSelectObject('oi', oi);
oiWindow;

%% 