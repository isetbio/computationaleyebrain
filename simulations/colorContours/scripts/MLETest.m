%% This script computes color discrimination contours with MLE
%
%
%
%  (HJ) Jan, 2014

%%
s_initISET;

%% Init Parameters
monitorName = 'LCD-Apple.mat';
wave = 380:4:780;
%refColor = [127 127 127]/255; 
refColor = [127 127 128]/255;
bgColor = [127 127 127]/255;

M = 64; N = 64; % scene size
scenePixels = [M N];
refImage = ones([scenePixels 3]);
for ii = 1 : 3
    refImage(:,:,ii) = refImage(:,:,ii) * refColor(ii);
end
display = displayCreate(monitorName);

%%
for ii = 1 : 3
    refImage(:,:,ii) = refImage(:,:,ii) * refColor(ii);
end
scene = sceneFromFile(refImage, 'rgb', [], monitorName, wave);
scene = sceneSet(scene, 'h fov', 1); % 1 degree

%% Initiate the human optics
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf = wvfSet(wvf,'zcoeffs',sample_mean);
wvf = wvfComputePSF(wvf);
oi  = wvf2oi(wvf,'shift invariant');
oi  = oiCompute(scene,oi);

%%
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', 0.05); % 50 ms
sensor = sensorSetSizeToFOV(sensor, 1,  scene, oi);
sensor = sensorSet(sensor,'noise flag',0); % noise free
sensor = sensorCompute(sensor, oi);

%%
LMS0 = [72 60 27.73]; %sort(unique(sensorGet(sensor,'photons')), 'descend');
LMS1 = [72 60 29];

%%
mu0 = LMS0 + sqrt(LMS0); sigma0 = diag(sqrt(LMS0)./[60 30 10]);
mu1 = LMS1 + sqrt(LMS1); sigma1 = diag(sqrt(LMS1)./[60 30 10]);
sigma = (sigma0 + sigma1) / 2;

%%
err = (det(sigma)^2/det(sigma0)/det(sigma1))^(-1/4)*exp(-1/8*(mu1-mu0)*inv(sigma)*(mu1-mu0)')/2;
acc = 1 - err;