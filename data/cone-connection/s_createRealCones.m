%% s_realHumanSensor
%    Create real human sensor structure
%
%  For real human data, the cone mosaic is not aligned as a grid and thus
%  we could not use sensorCreate('human') to create computational human
%  cone structure.
%  Instead, we create it step by step here.
%
%  (HJ) April, 2014


%% Init
wave = 400:10:700;
expTime = 0.001;

coneAperture = [2 2]*1e-6;

%% Create basic sensor structure
sensor = sensorCreate;
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'exp time', expTime);
sensor = sensorSet(sensor,'pixel',pixelCreate('human',wave));
sensor = sensorSet(sensor,'name', 'sensor-Human-real');

%% Load cone mosaic data
xy = importdata('cones_plantain_subset_locations.txt');
coneType = importdata('cones_plantain_subset_coneTypes.txt');

% The position data is from (120~220, 160~210)
% we transform it to handle the eccentricity
xy = xy - repmat([170 185],[length(xy) 1]);

% we only have L/M in the cone type data
% transform L to 2 and M to 3
coneType = strcmp(coneType, 'M') + 2;
density = sum(coneType-2)/length(coneType);
density = [0 1-density density 0];

% Set cone types
sensor.cols = 1; sensor.rows = length(coneType);
sensor = sensorSet(sensor, 'pattern', coneType);

%% Create cone structure
% Add the default lens structure
lens = lensCreate([], wave);
sensor = sensorSet(sensor, 'human lens', lens);

% Add the default macular structure
macular = macularCreate([], wave);
sensor = sensorSet(sensor, 'human macular', macular);

cone = coneCreate('human');
cone = coneSet(cone, 'spatial density', density);
sensor = sensorSet(sensor, 'human cone', cone);
fN = {'kBlack', 'rLong', 'gMiddle', 'bShort'};
sensor = sensorSet(sensor,'filterNames',fN);

pixel  = sensorGet(sensor,'pixel');
pixel  = pixelSet(pixel,'sizeSameFillFactor',coneAperture);
sensor = sensorSet(sensor,'pixel',pixel);

sensor = sensorSet(sensor,'cone locs',xy);
sensor = sensorSet(sensor,'cone type',coneType);

%% END