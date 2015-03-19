%% t_VernierAcuity
%    This turtorial script uses biological and computational methods to
%    explain vernier acuity (super acuity) in human vision
%
%    Vernier acuity (or positional acuity) is a measurement of sensitivity
%    of human eye in detecting mis-alignment of simple object (lines, etc.)
%
%    In this script, we compute the irradiance and human cone absorptions
%    for a scene with two mis-aligned lines. Then, we try to discriminate
%    the aligned and mis-aligned cases by using first order statistics and
%    machine learning classifiers
%
%  HJ/BW, ISETBIO TEAM, 2015

%% Init
%  Initialize a new session
ieInit;

%% Create Scene
%  In the section, we create a verneir scene radiance image by specifying a
%  image on some calibrated displays. This method can take some steps, but
%  it's more flexible. Another way to create vernier scenes is to call
%  sceneCreate('vernier') directly

% Create scene with RGB image on display
% set parameters
viewDist = 2; % viewing distance in meters
imgSz    = [200 80]; % image size in pixels
barColor = [1 1 1]; % RGB value for foreground bar
bgColor  = [0.5 0.5 0.5]; % Background color
barWidth = 3; % width of the bar in pixels
offset   = 1; % mis-alignment size in pixels

doSub = false; % rendering scene at pixel level (no subpixel rendering)
wave  = 400:10:700; % wavelength sample points
meanLum = [];  % adjustment to scene mean luminance - don't adjust it

% create display and linearize its gamma table
d = displayCreate('LCD-Apple');
d = displaySet(d, 'gamma', repmat(linspace(0, 1, 256)', [1 3]));
d = displaySet(d, 'viewing distance', viewDist);
vcAddObject(d); displayWindow;

% create image
img = zeros([imgSz 3]);
barIndx = (1:barWidth) - floor((barWidth+1)/2) + imgSz(2)/2;
for ii = 1 : 3 % loop over color primaries (RGB)
    img(:, :, ii) = bgColor(ii) * ones(imgSz); % set background
    img(:, barIndx, ii) = barColor(ii); % set foreground bar
end
imgA = img; imgM = img; % A for aligned, M for mis-aligned
imgM(1:imgSz(1)/2, :,:) = circshift(img(1:imgSz(1)/2, :,:), offset, 2);

vcNewGraphWin([], 'tall'); 
subplot(2,1,1); imshow(imgA); title('Aligned Image');
subplot(2,1,2); imshow(imgM); title('Misaligned Image');

% create a scene with the image and display
sceneA = sceneFromFile(imgA, 'rgb', meanLum, d, wave, doSub);
sceneM = sceneFromFile(imgM, 'rgb', meanLum, d, wave, doSub);

vcAddObject(sceneA); vcAddObject(sceneM); sceneWindow;

%% Compute Irradiance with Optics Wavefront
%    In this section, we compute the optical image (irradiance map) by
%    using human optics wavefront. Another good way to do this computation
%    is using a shift-invariant model with simulated chromatic abbretion
%    (Wandell, 1994).
%
%    If we just want a standard human optics model, we can simplify the
%    code as oiCraete('wvf human');

%  Load Zernike coefficient
pupilSize = 3; % pupil size in mm
zCoefs = wvfLoadThibosVirtualEyes(pupilSize);

%  Create wavefront structure
wvf = wvfCreate('wave', wave, 'zcoeffs', zCoefs, 'name', 'human optics');
wvf = wvfSet(wvf, 'calc pupil size', pupilSize); 

% Adjust for individuals
% Here, we use defocus as an example. For more adjustable entries, see
% wvfOSAIndexToVectorIndex

% ajust zernike coefficient for defocus 
% if we need to use defocus in diopters, use wvfDefocusDioptersToMicrons to
% do the conversion
zDefocus = -0.0104; 
wvf = wvfSet(wvf, 'zcoeffs', zDefocus, {'defocus'});

% compute psf and convert to optical image structure
wvf = wvfComputePSF(wvf);
oi = wvf2oi(wvf, 'human');
oi = oiSet(oi, 'name', sprintf('Human WVF %.1f mm', pupilSize));

% compute irradiance map (optical image)
oiA = oiCompute(sceneA, oi);
oiM = oiCompute(sceneM, oi);

vcAddObject(oiA); vcAddObject(oiM); oiWindow;

%% Compute Photon Absorptions of Human Cone Receptors
%    In this section, we compute the human cone absorption samples with
%    fixational eye movement.

%  set parameters
params.humanConeDensities = [0 0.6 0.3 0.1]; % cone spatial density KLMS
params.wave = wave; % wavelength
sampTime = 0.001; % sample time interval
params.totTime = 0.05; % set total time
params.emFlag = [1 0 0]; % eye movement type flag

%  create human sensor
sensor = sensorCreate('human', [], params);

%  adjust sensor size
sensor = sensorSetSizeToFOV(sensor, sceneGet(sceneA, 'fov'), sceneA, oiA);

%  set up eye movements
sensor = sensorSet(sensor, 'time interval', sampTime);
sensor = eyemoveInit(sensor, params);

%  randomize initial points
%% Analysis