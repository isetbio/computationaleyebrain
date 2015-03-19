%% s_eyemove
%    Compute cone absorptions with long range eye movements
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit;

%% Generate Eye Movements
I  = imread('eagle.jpg'); % load image
ts = 0.001; % Sample interval
nSamples = 10;

% Eye position at each sample time
% should change this one with more practical eye movement patterns
eyePos = 50 + round(bsxfun(@times, rand(nSamples, 2), ...
    [size(I,1) size(I,2)] - 100));


% If we know the actual size, we don't need to set fov
fov = 1;

%% Comptue cone absorption for each eye position
sensor = sensorCreate('human');
oi = oiCreate('human');

sensor = sensorSetSizeToFOV(sensor, fov);

p = zeros([sensorGet(sensor, 'size') nSamples]);

for ii = 1 : nSamples
    pos = eyePos(ii, :);
    
    % crop a region around eye pos
    img = I(pos(1)-50:pos(1)+50, pos(2)-50:pos(2)+50, :);
    
    % create scene
    scene = sceneFromFile(img, 'rgb', [], 'LCD-Apple');
    scene = sceneSet(scene, 'fov', fov); 
    
    % compute irradiance map
    oi = oiCompute(scene, oi);
    
    % compute cone absorptions
    sensor = sensorCompute(sensor, oi);
    
    % Get photon absorptions
    p(:,:,ii) = sensorGet(sensor, 'photons');
end

%% Visualize
%  we can make a small video based on photon absorptions here