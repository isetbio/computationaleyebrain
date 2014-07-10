%%  Test for speed
%
%
%  (HJ) ISETBIO TEAM, 2014

%% Init
I = rand(1000);
d = displayCreate('OLED-Sony');

%% Scene
%  Create the scene from image
%  The scene is 10 x 10 deg
tic;
fov = 10;
scene = sceneFromFile(I, 'rgb', [], d);
scene = sceneSet(scene, 'h fov', fov);
toc;

%% Optical Image
tic;
oi = oiCreate('wvf human');
oi = oiCompute(scene, oi);
toc;

%% Sensor
tic;
% create sesnor
sensor = sensorCreate('human');

% adjust fov
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

% set eye movement
em = emCreate;
sensor = sensorSet(sensor, 'eye movement', em);
sensor = sensorSet(sensor, 'positions', zeros(100, 2));
sensor = emGenSequence(sensor);

% Compute cone absorptions
sensor = coneAbsorptions(sensor, oi);
toc;

% Add a little adaptations
tic;
sensor = coneAdapt(sensor);
toc;