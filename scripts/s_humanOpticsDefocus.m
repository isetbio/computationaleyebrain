%% s_humanOpticsDefocus
%    Generate psf and show effect of defocus in human optics
%    This script is based on changing Zernike coefficient
%    For oi created by other method, see s_opticsDefocus
%
% (HJ) ISETBIO TEAM, 2015

%% Init & Set Parameters
ieInit; % initialize a new ISET session

pupilSize = 3; % pupil size diameter in mm, can choose from 7.5, 6, 4.5, 3
defocus   = 2; % defocus in diopters

fov = 1; % field of view of the scene

%% Load Zernike Coefficient
% load standard 
zCoefs = wvfLoadThibosVirtualEyes(pupilSize);

% adjust for defocus
microns = wvfDefocusDioptersToMicrons(defocus, pupilSize);
zCoefs(4) = zCoefs(4) + microns;

%% Create Scene and Human optics
%  create a scene
scene = sceneFromFile('eagle.jpg', 'rgb', [], 'OLED-Sony');
scene = sceneSet(scene, 'fov', fov);

% create human optics structure
oi = oiCreate('wvf human'); % in focus
oiD = oiCreate('wvf human', [], [], [], pupilSize, zCoefs); % defocused

% compute irradiance map
oi = oiCompute(scene, oi);
oiD = oiCompute(scene, oiD);

% visualize
vcAddObject(oi); vcAddObject(oiD); oiWindow;

%% Get and Plot PSF (or OTF)
%  We will get the point spread function from the optics and plot it
%  Actually, we could directly call plotOI(oi, 'psf 550') to get the plot

%  Get point spread from optics and compare
psf  = oiGet(oi,  'optics psf data', 550);
psfD = oiGet(oiD, 'optics psf data', 550);

psf_support  = oiGet(oi,  'optics psf support');
psf_supportD = oiGet(oiD, 'optics psf support');

% visualize
vcNewGraphWin;
subplot(1,2,1); mesh(psf_support{1},  psf_support{2},  psf);
xlabel('Position (um)'); ylabel('Position (um)'); 
zlabel('Irradiance (relative)'); title('Point Spread (in focus)');

subplot(1,2,2); mesh(psf_supportD{1}, psf_supportD{2}, psfD);
xlabel('Position (um)'); ylabel('Position (um)'); 
zlabel('Irradiance (relative)'); title('Point Spread (Defocused)');

%% Scene with depth
%  In this section, we will apply oi to scene with depth informtion. 
%  In the computation, we assume we are focusing on either background or
%  foreground (discretize the depth map to two planes)
%
%  Note that the computation here is just an approximation and can be
%  inaccurate. A better implementation can be found in s_opticsDepthScene
%  in ISET. For the best solution, you will need to some specially designed
%  ray-tracing algorithms (scene3D project)
%
%  There are functions can be called to complete this. But here, we just
%  use plane code to let people know how to do it manually. To see the more
%  professional ISET coding, see s_opticsDepthScene
%
%  see also: s_opticsDepthScene, oiDepthCompute

%  load scene with depth information
load(fullfile(isetRootPath,'data','scenes','piano3d.mat'));
scene = sceneSet(scene,'fov', fov);

% for this scene, we simply assign it as foreground and background (well,
% not that good, but let's try it)
depthMap = sceneGet(scene, 'depth map');
[idx, C] = kmeans(depthMap(:), 2); % classify to fore / background
idx = reshape(idx, size(depthMap)) - 1; % make idx 0/1 instead of 1 and 2

% visualize the depth map
vcNewGraphWin;
subplot(1,2,1); imagesc(depthMap); title('True Depth Map');
subplot(1,2,2); imagesc(idx); title('Simplified Depth Map');

% Create optics for the defocused
zCoefs  = wvfLoadThibosVirtualEyes(pupilSize);
imgDist = 1/(1/oiGet(oi, 'focal length') - 1/min(C));
defocus = opticsDepthDefocus(max(C), oiGet(oi, 'optics'), imgDist);
microns = wvfDefocusDioptersToMicrons(defocus, pupilSize);
zCoefs(4) = zCoefs(4) + microns;
oiD = oiCreate('wvf human', [], [], [], pupilSize, zCoefs);

% Compute in-focus and defocused image
oi  = oiCompute(scene, oi);  % compute in focus
oiD = oiCompute(scene, oiD); % compute defocused

% Combine the optical irradiance map
% Actually, for simplicity, you could directly do combination in rgb images
p  = oiGet(oi,  'photons'); % in focus photon image
p  = p(26:225, 26:225, :);  % get rid of the border
pD = oiGet(oiD, 'photons'); % defocused photon image
pD = pD(26:225, 26:225, :);

pCombined = bsxfun(@times, p, idx) + bsxfun(@times, pD, 1 - idx); % combine
img_fore = imageSPD(pCombined, [], [], [], [], -1);

% Now, focus on background and defocus foreground
zCoefs  = wvfLoadThibosVirtualEyes(pupilSize);
imgDist = 1/(1/oiGet(oi, 'focal length') - 1/max(C));
defocus = opticsDepthDefocus(min(C), oiGet(oi, 'optics'));
microns = wvfDefocusDioptersToMicrons(defocus, pupilSize);
zCoefs(4) = zCoefs(4) + microns;
oiD = oiCreate('wvf human', [], [], [], pupilSize, zCoefs);

% comptue defocused image and combine
oiD = oiCompute(scene, oiD);
pD  = oiGet(oiD, 'photons');
pD = pD(26:225, 26:225, :);

pCombined = bsxfun(@times, p, 1-idx) + bsxfun(@times, pD, idx); % combine
img_back = imageSPD(pCombined, [], [], [], [], -1);

% visualize the combined image
vcNewGraphWin; 
subplot(1,2,1); imshow(img_fore); title('foreground in focus');
subplot(1,2,2); imshow(img_back); title('background in focus');