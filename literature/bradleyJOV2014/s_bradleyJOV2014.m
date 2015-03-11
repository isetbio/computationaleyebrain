%% s_bradleyJOV2014
%
%    This script replicates Bradley et al. paper in JOV 2014
%
%       Retina-V1 model of detectability across the visual field
%
%    This script compares the components, performance and results between
%    their implementation and ISETBIO implementation
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit;

%% Compare human optics
%  create human optics in ISETBIO
oi = oiCreate('human'); % standard human optics

%  get otf function of human optics
otf = oiGet(oi, 'optics otf', 540); % otf at 540 nm

%  get frequency support in cycles / degree
fx = oiGet(oi, 'optics otf fx', 'cycles/deg'); % frequency support
fx = fx(fx >= 0); % use non-negative frequencies only

%  plot mtf function
vcNewGraphWin; grid on; hold on;
plot(fx, otf(1, 1:size(otf, 2)/2), '--r');

% plot MTF in Bradley et al. paper
% MTF function is given by function:
%   MTF(f) = 0.78 exp(-0.172f) + 0.22 exp(-0.037f)
MTF = @(f) 0.78 * exp(-0.172*f) + 0.22 * exp(-0.037*f);
plot(fx, MTF(fx), '-b');

% add label, title and lengend
xlabel('Spatial Frequency (cycles/deg)'); ylabel('MTF');
title('Comparison of Optics'); legend('ISETBIO Optics', 'Bradley Paper');

% The MTF function used in the paper is not wavelength dependent and it
% seems like it's close to MTF in ISETBIO for wavelength around 540 nm

%% Light adaptation
%  In ISETBIO, light adaptation in cone photoreceptors is implemented in
%  coneAdapt. Light adaption in RGC has not been implemented yet
%
%  In Bradley et al. paper, the light adaptation is by local luminance gain
%  controll, which scales the response by its neighborhood mean luminance
%
%  

%% Ganglion cell
%    A simple center surround mechanism is implemented in the paper
%    Here, we will re-write that code in an easy to understand way
%
%    Filter is defined as difference of Guassian: 
%      D(y;x) = wc* Gauss_c(y;x) - (1-wc)*Gauss_s(y;x)
%

% init parameters
fov = 3; % field of view
wc  = 0.53; % weight for center gaussian
ks  = 10.1; % standard deviation for surround filter
kc  = 1;

% create a simple scene and compute cone absorptions
scene = sceneFromFile('eagle.jpg', 'rgb', [], 'LCD-Apple');
scene = sceneSet(scene, 'h fov', fov); % set field of view to smaller value
oi = oiCreate('wvf human'); % human optics
oi = oiCompute(scene, oi);
sensor = sensorCreate('human'); % create human cone mosaic
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi); % adjust sensor size
sensor = sensorCompute(sensor, oi);

% visualize sensor image
vcAddObject(sensor); sensorWindow;

p = sensorGet(sensor, 'photons'); % get photons

% create filter
% actually, we could call coneComputeCenterSurround(p) instead
fc = fspecial('gaussian', 20, kc); % center filter
fs = fspecial('gaussian', 20, ks); % surround filter
p_filtered = wc * imfilter(p, fc) - (1-wc) * imfilter(p, fs);

% visualize
vcNewGraphWin;
pf = p_filtered;
coneType = sensorGet(sensor, 'cone type');
for ii = 2:4
    indx = (coneType == ii);
    pf(indx) = p_filtered(indx) - min(p_filtered(indx));
    pf(indx) = pf(indx) / max(pf(indx));
end
imshow(pf);

