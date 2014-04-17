%% s_coneTypeLearning
%
%    This is the script that uses unspervised learning method to classify
%    cones into groups
%
%  (HJ) April, 2014

%% Init
imageDir = '~/Pictures/cd02A/';
imageS = dir([imageDir '*.JPG']);

fov = 1; % 1 deg
density = [0 0.6 0.3 0.1];

d = displayCreate('LCD-Apple');

%% Create human optics and sensor
%  Create a typical human lens
wave   = 380 : 4 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);
wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

%  Create a human sensor
params.humanConeDensities = density;
sensor = sensorCreate('human', [], params);
sensor = sensorSet(sensor, 'exp time', 0.05);
sensor = sensorSet(sensor, 'noise flag', 1); %only photon noise

%% Compute cone absorption for all natural images
imFName = fullfile(imageDir, imageS(1).name);
imData = im2double(imread(imFName));
imData = imresize(imData, 0.1);
scene = sceneFromFile(imData, 'rgb', [], d);
scene = sceneSet(scene, 'fov', fov);

oi = oiCompute(scene, oi);
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

sz = sensorGet(sensor, 'size');
cc = round(sz/2);
voltsImage = zeros([21*21 length(imageS)]);

for ii = 1 : length(imageS)
    fprintf('Loading image %d...\n', ii);
    imFName = fullfile(imageDir, imageS(ii).name);
    imData = im2double(imread(imFName));
    imData = imresize(imData, 0.1);
    imData(imData < 0) = 0;
    imData(imData > 1) = 1;
    scene = sceneFromFile(imData, 'rgb', [], d);
    scene = sceneSet(scene, 'fov', fov);
    
    oi = oiCompute(scene, oi);
    
    sensor = sensorCompute(sensor, oi);
    tmp = sensorGet(sensor, 'photons');
    tmp = tmp((cc(1)-10):(cc(1)+10), (cc(2)-10):(cc(2)+10));
    voltsImage(:,ii) = tmp(:);
end

%% Unsupervised learning
%  Compute correlation matrix
corMatrix = corr(voltsImage');

% Multi-dimensional scaling
[Pos,E] = cmdscale(-log((corMatrix+1)/2));

% Kmeans
predLabels = reshape(kmeans(Pos, 3), [21 21]);

%% Visualize
trueLabels = sensorGet(sensor, 'cone type');
trueLabels = trueLabels(cc(1)-10:cc(1)+10, cc(2)-10:cc(2)+10);
figure;
subplot(1,2,1); imagesc(trueLabels);
subplot(1,2,2); imagesc(predLabels);
fprintf('Accuracy: %f\n', sum(abs(predLabels(:)-trueLabels(:)))/441);