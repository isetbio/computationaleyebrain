%% t_scene2coneAbsorption
%
%  This is the tutorial script for Psych 221 project.
%
%  In this script, we illustrate how to build up a scene based on a
%  calibrated display and generate cone absorption samples
%
% (HJ) Feb, 2014

%% Init Parameters
if notDefined('ppi'), ppi = 500; end            % points per inch
if notDefined('imgFov'), imgFov = [60 60]/60; end % visual angle (degree)
if notDefined('nFrames'), nFrames = 2000; end   % Number of samples

vDist  = 1.0;                                   % viewing distance (meter)
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/39.37/vDist);     % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
img = ones(imgSz)*0.5;            % Init to black
img(:, round(imgSz(1)/2)) = .99;  % Draw vertical straight line in middle
img(1:round(imgSz(1)/2), :) = circshift(img(1:round(imgSz(1)/2), :),[0 1]);

% Create scene from file
scene = sceneFromFile(img, 'rgb', [], display);

% set scene fov
scene = sceneSet(scene, 'h fov', imgFov);
scene = sceneSet(scene, 'distance', vDist);

% Visualize scene
vcAddAndSelectObject('scene', scene); sceneWindow;

%% Create Human Lens
%  Create a typical human lens
oi = oiCreate('wvf human');

% Compute optical image
oi = oiCompute(scene, oi);

% Visualize optical image
vcAddAndSelectObject('oi', oi); oiWindow;

%% Create Sensor and Compute Samples

%  Create human sensor
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, imgFov, scene, oi); % set fov
sensor = sensorSet(sensor, 'exp time', 0.05); % integration time: 50 ms

%  Set exposure time to 1 ms
sensor = sensorSet(sensor, 'exp time', 0.001);

% Set up the eye movement properties
sensor = eyemoveInit(sensor);

% Compute the cone absopritons
sensor = coneAbsorptions(sensor, oi);

vcAddObject(sensor); sensorWindow;

% Store the photon samples
pSamples = double(sensorGet(sensor, 'photons'));

%% END