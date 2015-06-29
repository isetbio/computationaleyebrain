%% t_Andres
%    This is a starting file for Andres project for course Psych 221
%    In the script, we will use ISETBIO to generate a scene file for the
%    psychophysics image data set and compute the irradiance map after
%    human optics. We will compute cone absorption afterwards and make
%    analysis there based on the photon absorption rate
%
%    HJ knows how to do the computation
%
%  (HJ) ISETBIO TEAM

%% Init
ieInit;

%% Create display
d = displayCreate('OLED-Sony');
I = imread('structured_D62_s.jpg');

%% Create scene
scene = sceneFromFile(I, 'rgb', [], d);
scene = sceneSet(scene, 'fov', 2);

vcAddObject(scene); sceneWindow;

%% Compute irradiance map
oi = oiCreate('wvf human');
oi = oiCompute(scene, oi);

vcAddObject(oi); oiWindow;

%% Compute cone absorptions
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, sceneGet(scene, 'fov'), scene, oi);
sensor = sensorCompute(sensor, oi);

vcAddObject(sensor); sensorWindow('scale', 1);

%% Do Analysis