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

%% Create the display

% In this example we impose a linear gamma table, though
% in general it could be the default or anything.
d = displayCreate('LCD-Apple');
% d = displaySet(d, 'gamma', repmat(linspace(0, 1, 256)', [1 3]));

viewDist = 2; % viewing distance in meters
d = displaySet(d, 'viewing distance', viewDist);

%% Create the the RGB image
imgSz    = [200 80]; % image size in pixels
barColor = [1 1 1]; % RGB value for foreground bar
bgColor  = [0.5 0.5 0.5]; % Background color
barWidth = 3; % width of the bar in pixels
offset   = 1; % mis-alignment size in pixels, converted to deg later

doSub = false; % rendering scene at pixel level (no subpixel rendering)
wave  = 400:10:700; % wavelength sample points
meanLum = [];  % adjustment to scene mean luminance - don't adjust it

img = zeros([imgSz 3]);
barIndx = (1:barWidth) - floor((barWidth+1)/2) + imgSz(2)/2;
for ii = 1 : 3 % loop over color primaries (RGB)
    img(:, :, ii) = bgColor(ii) * ones(imgSz); % set background
    img(:, barIndx, ii) = barColor(ii); % set foreground bar
end
imgA = img; imgM = img; % A for aligned, M for mis-aligned
imgM(1:imgSz(1)/2, :,:) = circshift(img(1:imgSz(1)/2, :,:), [0 offset 0]);
        
vcNewGraphWin([], 'tall'); 
subplot(2,1,1); imshow(imgA); title('Aligned Image');
subplot(2,1,2); imshow(imgM); title('Misaligned Image');

%% Create Vernier Scene in its full radiance representation
%
%  In the section, we create a vernier scene radiance image by specifying a
%  image on some calibrated displays. This method makes each of the
%  parameters explicit. This provides flexibility for scripting. 
%
%  Another way to create a vernier scene with default parameters is to call 
%
%     scene = sceneCreate('vernier');
%     vcAddObject(scene); sceneWindow;
%
% directly

% Create a scene with the image using the display parameters
% The scene spectral radiance is created using the RGB image and the
% properties of the display.
sceneA = sceneFromFile(imgA, 'rgb', meanLum, d, wave, doSub); % aligned
sceneM = sceneFromFile(imgM, 'rgb', meanLum, d, wave, doSub); % mis-aligned

vcAddObject(sceneA); vcAddObject(sceneM); sceneWindow;

%% Examine some of the scene properties

% This is the scene luminance at the different sample points on the display
sz = sceneGet(sceneM,'size');
plotScene(sceneM,'luminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
plotScene(sceneM,'radiance hline',[1,round(sz(1)/2)]);


%% Compute Irradiance with Optics Wavefront
%    In this section, we compute the optical image (irradiance map) by
%    using human optics wavefront.
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

%% Examine the point spread function at different wavelengths

% This is the standard model based on the Thibos AO estimates of the human
% wavefront function
plotOI(oi,'psf 550');
nPoints = 50;
plotOI(oi,'ls wavelength',[],nPoints);


%% compute irradiance map (optical image)
oiA = oiCompute(sceneA, oi);
oiM = oiCompute(sceneM, oi);

vcAddObject(oiA); vcAddObject(oiM); oiWindow;

% Another way to do this computation is using the chromatic aberration in
% the Marimont and Wandell model (1994, JOSA).  You can create that oi
% structure simply by calling oi = oiCreate('human');

%% Examine the irradiance at the retina, prior to absorption

% This is the scene luminance at the different sample points on the display
sz = sceneGet(oiA,'size');
plotOI(oiA,'illuminance hline',[1,round(sz(1)/2)]);

% This is the full spectral radiance on the same line
plotOI(oiA,'irradiance hline',[1,round(sz(1)/2)]);
view(89.5,55);
% Notice that the short-wavelength light is spread a lot more at the
% retinal surface than the middle wavelength light.
% The long-wavelength is spread a little bit more.

%% Compute Photon Absorptions of Human Cone Receptors
%    In this section, we compute the human cone absorption samples with
%    fixational eye movement.

%  set retinal cone parameters, but no eye movements
params.humanConeDensities = [0 0.6 0.3 0.1]; % cone spatial density KLMS
params.wave = wave; % wavelength

%  create human sensor
cones = sensorCreate('human', [], params);

%  adjust sensor size
cones = sensorSetSizeToFOV(cones, sceneGet(sceneA, 'fov'), sceneA, oiA);

cones = sensorCompute(cones,oiM);
vcAddObject(cones); sensorWindow('scale',1);

%% Static analysis, without eye movements

% Plot the average of the top and bottom half, looking for displacement
% of peak
e  = sensorGet(cones,'photons');
sz = sensorGet(cones,'size');
midSensor = floor(sz(2)/2) + 1;

midRow    = floor((sz(1)/2));
topSensor = 1:midRow;
botSensor = (midRow+1):sz(1);
topE  = mean(e(topSensor,:),1)/max(e(:));
botE  = mean(e(botSensor,:),1)/max(e(:));

% Graph the top and bottom
vcNewGraphWin;
plot([topE(:),botE(:)])
% line([midSensor,midSensor],[min(e(:)) max(e(:))],'color','r')
grid on; legend('top','bottom')
title(sprintf('Bar offset %0.f (arc sec)',sceneGet(sceneM,'degrees per sample','arcsec')))
ylabel('Normalized cone absorptions');
xlabel('Position (um)');

%% Analysis with eye movements

% Set up eye movement parameters
sampTime = 0.001; % sample time interval
cones = sensorSet(cones, 'time interval', sampTime);
cones = sensorSet(cones,'exp time',sampTime);

% Eye movement type flag
% Each entry is a type of eye movement
% Tremor, drift, microsaccade

params.emFlag = [1 1 1]; 

% This is how long we run the sequence for
params.totTime = 0.500;  % set total time, 200 ms

% Attach the eye movement parameters to the cone array object
cones = eyemoveInit(cones, params);

% For demo, we enlarge the eye movement tremo
tremor = sensorGet(cones,'em tremor');
tremor.amplitude = tremor.amplitude*5;
cones = sensorSet(cones,'em tremor',tremor);

% Generate the eye movement sequence
cones = emGenSequence(cones);

% Show the eye movement positions
p = sensorGet(cones,'sensor positions');

vcNewGraphWin;
plot(p(:,1),p(:,2),'-o')
xlabel('Cone position');
ylabel('Cone Position')
set(gca,'xlim',[-15 15],'ylim',[-15 15]);
grid on

%% Show the cone absorptions each millisecond

% Notice that when we make an eye movement video we call coneAbsorptions,
% not sensorCompute.  The coneAbsorptions function is an interface to
% sensorCompute
cones = coneAbsorptions(cones,oiM);

% Get the time series out from the cone photon data
absorptions = sensorGet(cones,'photons');

% This should work, but it depends on having a later version of Matlab than
% 2013b
%   mplay(absorptions,'intensity',10);

vcNewGraphWin;
vObj = VideoWriter('coneAbsorptions.avi');
open(vObj);
colormap(gray);
nframes = size(absorptions,3);
% Record the movie
mx = max(absorptions(:));
for j = 1:nframes 
    image(absorptions(:,:,j)*(255/mx)); drawnow;
    F = getframe;
    writeVideo(vObj,F);
end
close(vObj);
fprintf('Max cone absorptions %.0f\n',mx);

%% Temporal integration from Rieke applied to the cone absorptions
[cones,adaptedData] = coneAdapt(cones,'rieke');
adaptedData = ieScale(adaptedData,0,1);

vcNewGraphWin;
vObj = VideoWriter('coneVoltage.avi');
open(vObj);
colormap(gray);
nframes = size(adaptedData,3);
% Record the movie
mx = max(adaptedData(:));
for j = 1:nframes 
    image(255*adaptedData(:,:,j)); drawnow;
    F = getframe;
    writeVideo(vObj,F);
end
close(vObj);


%%