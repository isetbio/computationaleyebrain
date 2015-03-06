%% s_deltaEVariabilityForUniform
%    This script computes the variance for cone absorptions in units of
%    deltaEab
%  
%    The intuition for this calculation is to see how large deltaE is
%    'statistically significant'
%
%  (HJ) ISETBIO TEAM, 2015

%% Init & Parameters
ieInit; % initialize a new ISET session

fov = 2;  % field of view in degrees
lum = 10; % mean luminance level for scene in cd/m^2
nSample = 1000; 

%% Create a uniform scene
scene = sceneCreate('uniformEqualEnergy');
scene = sceneSet(scene, 'fov', fov);
scene = sceneAdjustLuminance(scene, lum);

%% Compute cone absorptions
%  comptue irradiance
oi = oiCreate('wvf human');
oi = oiCompute(scene, oi);

%  compute cone absorptions
sensor = sensorCreate('human');
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);
sensor = sensorCompute(sensor, oi);

%% Analysis distribution in deltaEab
%  Get photons
p = sensorGet(sensor, 'photons');

%  Compute mean LMS values (noise free)
coneType = sensorGet(sensor, 'cone type');
meanLMS  = zeros(3,1);
for ii = 1 : 3
    meanLMS(ii) = round(mean(p(coneType == ii+1)));
end

%  Sample and compute XYZ values
for ii = 1 : nSamples
    
end