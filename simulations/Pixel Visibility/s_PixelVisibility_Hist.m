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

% 2x2 pixel
% img(1:4:end, 1:4:end, :) = 1;
% img(2:4:end, 1:4:end, :) = 1;
% img(1:4:end, 2:4:end, :) = 1;
% img(2:4:end, 2:4:end, :) = 1;

% 3x3 pixel
% img(1:9:end, 1:9:end, :) = 1;
% img(2:9:end, 1:9:end, :) = 1;
% img(3:9:end, 1:9:end, :) = 1;
% img(1:9:end, 2:9:end, :) = 1;
% img(2:9:end, 2:9:end, :) = 1;
% img(3:9:end, 2:9:end, :) = 1;
% img(1:9:end, 3:9:end, :) = 1;
% img(2:9:end, 3:9:end, :) = 1;
% img(3:9:end, 3:9:end, :) = 1;

% every other pixel
img(1:2:end, 1:2:end, :) = 1;


imgFov = 40 * displayGet(d, 'deg per pixel');

doSub = true;   % we do care the subpixel structure
oSample = 20;   % upsampling rate for each pixel

% create scene structure
scene = sceneFromFile(img, 'rgb', [], d, [], doSub, [], oSample);

% uniform scene for comparison
% p  = mean(mean(sceneGet(scene, 'photons'))); % average photon rate
% sz = sceneGet(scene, 'size');
% sceneU = sceneSet(scene, 'photons', repmat(p, [sz 1]));

% visualize
vcAddObject(scene); sceneWindow;

% The analysis is based on the perfectly uniform field with matched
% luminance and spd. It's not uniform image on that display
%
% We could also compare with the uniform image on display

imgU = ones(40, 40, 3);
sceneU = sceneFromFile(imgU, 'rgb', [], d, [], doSub, [], oSample);
sceneU = sceneAdjustLuminance(sceneU, sceneGet(scene, 'mean luminance'));
vcAddObject(sceneU); sceneWindow;

%% Optical Irradiance
%  In this section, we create the optics stucture for standard human
%  observer and compute the corresponding irradiance map

% create optics structure
oi = oiCreate('wvf human');

% compute irradiance map
oiU = oiCompute(oi, sceneU); % irradiance map for uniform field
oi  = oiCompute(oi, scene);  % irradiance map for image on display

% oiFov = oiGet(oi, 'fov');

% visualize
% vcAddObject(oi); oiWindow;
% vcAddObject(oiU); oiWindow;

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

% a cone is about 2 microns in size 
% a blur circle is about 3-5 cones
% so we want 10 micron sFOV
% sFOV is in degrees of visual angle so convert 
% One degree of visual angle = 330 microns 
% so 10 microns = 1/33 degree or 0.0303 degrees
sFOV = 0.0303
% sensor = sensorSetSizeToFOV(sensor, sceneGet(scene, 'h fov'), scene, oi);
sensor = sensorSetSizeToFOV(sensor, sFOV, scene, oi);
sensor = sensorSet(sensor, 'sensor positions', zeros(1000, 2));

% compute cone absorptions
sensor = coneAbsorptions(sensor, oi);  
sensorU = coneAbsorptions(sensor, oiU);

% visualize
% vcAddObject(sensor); sensorWindow;
% vcAddObject(sensorU); sensorWindow;

%% Analysis
%  In this section, we plot the histogram of each cone type

coneType = sensorGet(sensor, 'cone type');
p = sensorGet(sensor, 'photons');
pU = sensorGet(sensorU, 'photons');

figure;
% for ii = 2 : 4
%     subplot(1,3,ii-1); hold on;
%     
%     % plot histogram of image on display
%     indx = (coneType == ii);
%     [f, x] = hist(p(repmat(indx, [1 1 size(p,3)])), 20);
%     bar(x, f/trapz(x,f), 'g');
%     
% %     % plot theoretical curve
% %     meanP = median(pU(indx));
% %     plot(x, 1/sqrt(2*pi*meanP)*exp(-0.5*(x-meanP).^2/meanP), ...
% %         '--r','lineWidth', 2);
%     
%     % plot the uniform display image case
%     [f, x] = hist(pU(repmat(indx, [1 1 size(p,3)])), 20);
%     bar(x, f/trapz(x,f), 'r');
% end


for ii = 2 : 4
    subplot(1,3,ii-1); hold on;
    
    % plot histogram of image on display
    indx = (coneType == ii);
    [f, x] = hist(p(repmat(indx, [1 1 size(p,3)])), 20);
    h = bar(x, f/trapz(x,f), 'g');
    set(h,'facecolor','g','edgecolor','g')
%     h.FaceColor = [0 0.5 0.0];
    h.EdgeColor = 'g';
    
%     % plot theoretical curve
%     meanP = median(pU(indx));
%     plot(x, 1/sqrt(2*pi*meanP)*exp(-0.5*(x-meanP).^2/meanP), ...
%         '--r','lineWidth', 2);
    
    % plot the uniform display image case
    [f, x] = hist(pU(repmat(indx, [1 1 size(p,3)])), 20);
    h = bar(x, f/trapz(x,f), 'r');
    set(h,'facecolor','r','edgecolor','r')
%     h.FaceColor = [0 0.0 0.5];
    h.EdgeColor = 'r';
end


y = randn(10000,1);
[n,x] = hist(y,30);
h = bar(x,n);
set(h,'facecolor','r','edgecolor','b')
h.FaceColor = [0 0.0 0.5];
h.EdgeColor = 'r';