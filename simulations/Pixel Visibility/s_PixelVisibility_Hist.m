%% s_PixelVisibility_Hist
%
%  This script computes the cone absorption for display patch in standard
%  human cone mosaic.
%
%  For analysis, we plot the histogram of the theoretical uniform field and
%  display pixeleted design
%
%  See also:
%    coPixelVisibilityThreshold
%
%  (HJ) ISETBIO TEAM, 2015

%% Init Session
ieInit; % Initialize a new ISET session

%% Set Up Display
%  In this section, we create a display structure with calibrated data.
%  Also, we set up the viewing distance of experiment settings.

vDist = 1.5; % viewing distance

d = displayCreate('LCD-Gechic'); % create display structure
d = displaySet(d, 'viewing distance', vDist); % set viewing distance

% visualize
% vcAddObject(d); displayWindow;

%% Scene Radiance
%  In this section, we create a scene radiance map for certain image shown
%  on that calibrated display

img = zeros(40,40,3); % image to be shown on the display
% img(1:4:end, 1:4:end, 3) = 1;
% img(2:4:end, 2:4:end, 3) = 1;

doSub = true;   % we do care the subpixel structure
oSample = 20;   % upsampling rate for each pixel

% create scene structure
scene = sceneFromFile(img, 'rgb', [], d, [], doSub, [], oSample);

% uniform scene for comparison
p  = mean(mean(sceneGet(scene, 'photons'))); % average photon rate
sz = sceneGet(scene, 'size');
sceneU = sceneSet(scene, 'photons', repmat(p, [sz 1]));

% visualize
% vcAddObject(scene); sceneWindow;

%% Optical Irradiance
%  In this section, we create the optics stucture for standard human
%  observer and compute the corresponding irradiance map

% create optics structure
oi = oiCreate('wvf human');

% compute irradiance map
oiU = oiCompute(oi, sceneU); % irradiance map for uniform field
oi  = oiCompute(oi, scene);  % irradiance map for image on display

% visualize
% vcAddObject(oi); oiWindow;

%% Cone Absorption
%  In this section, we compute the cone absorption for standard human
%  observer

cone = coneCreate;
cone.spatialDensity = [0 0.6 0.3 0.1]; % cone density for K,L,M,S
expTime = 0.05; % cone integration time
sampleTime = 0.001; % sample interval

% create human cone mosaic
sensor = sensorCreate('human', cone);
sensor = sensorSet(sensor, 'exp time', expTime);
sensor = sensorSet(sensor, 'sample time interval', sampleTime);
sensor = sensorSetSizeToFOV(sensor, sceneGet(scene, 'h fov'), scene, oi);

sensor = sensorSet(sensor, 'sensor positions', zeros(1000, 2));

% compute cone absorptions
sensor = coneAbsorptions(sensor, oi);
sensorU = sensorComputeNoiseFree(sensor, oiU);

% visualize
% vcAddObject(sensor); sensorWindow;

%% Analysis
%  In this section, we plot the histogram of each cone type

coneType = sensorGet(sensor, 'cone type');
p = sensorGet(sensor, 'photons');
pU = sensorGet(sensorU, 'photons');

figure;
for ii = 2 : 4
    subplot(1,3,ii-1); hold on;
    
    % plot histogram of image on display
    indx = (coneType == ii);
    [f, x] = hist(p(repmat(indx, [1 1 size(p,3)])), 20);
    bar(x, f/trapz(x,f));
    
    % plot theoretical curve
    meanP = median(pU(indx));
    plot(x, 1/sqrt(2*pi*meanP)*exp(-0.5*(x-meanP).^2/meanP), ...
        '--r','lineWidth', 2);
end