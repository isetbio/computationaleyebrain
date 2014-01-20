%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
ppi = 500;          % Display property: points per inch
imgSz = [32 32];    % Row / columns of image in pixels

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
scene = cell(2, 1);
img = ones(imgSz)*0.5;            % Init to black
img(:, round(imgSz(2)/2)) = .99;  % Draw vertical straight line in the middle
scene{1} = sceneFromFile(img, 'rgb', [], display); % create scene

img(1:imgSz(1)/2, :) = circshift(img(1:imgSz(1)/2, :), [0 1]);
%img = circshift(img, [0 1]);
scene{2} = sceneFromFile(img, 'rgb', [], display);

% set scene fov
dist = 1.0;
fov = 2*atand(max(imgSz)*25.4 / ppi / 1000 /dist/2);

for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', fov);
    scene{ii} = sceneSet(scene{ii}, 'distance', dist);
end

%% Create Sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', 0.05);

%% Create Human Lens
%  Create a typical human lens
wave   = 380 : 4 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);
wvf    = wvfComputePSF(wvf);
oi     = wvf2oi(wvf,'shift invariant');

% Compute optical image
% Actually, we could wait to compute it in coneSamples
% But, we compute it here to do sanity check
vcAddAndSelectObject('scene', scene{1});
OIs{1} = oiCompute(scene{1}, oi);
vcAddAndSelectObject('scene', scene{2});
OIs{2} = oiCompute(scene{2}, oi);
% vcAddAndSelectObject('oi', OIs{1}); oiWindow;
% vcAddAndSelectObject('oi', OIs{2}); oiWindow;

%% Compute cone absorption (noise free)
sensor = sensorSet(sensor, 'exp time', 0.05);
sensor = sensorSetSizeToFOV(sensor, fov, scene{1}, OIs{1});
sensor = sensorComputeNoiseFree(sensor, OIs{1});
p1 = double(sensorGet(sensor, 'photons')); 

sensor = sensorSetSizeToFOV(sensor, fov, scene{2}, OIs{2});
sensor = sensorComputeNoiseFree(sensor, OIs{2});
p2 = double(sensorGet(sensor, 'photons'));

coneType = sensorGet(sensor, 'cone type');
coneL1 = zeros(size(p1, 2), 1); coneL2 = coneL1;
for curCol = 1 : size(p1, 2)
    coneL1(curCol) = mean(p1(coneType(:,curCol) == 2,curCol));
    coneL2(curCol) = mean(p2(coneType(:,curCol) == 2,curCol));
end

% plot
vcNewGraphWin; plot(coneL1); hold on; plot(coneL2, 'r');
title('Cone absorption for a horizontal line in two scenes');

%% Add some Gaussian blur (for eye-movement)
G = fspecial('gaussian', [16 1], 8); % eye movement is around 8 cones wide
coneLG1 = imfilter(coneL1, G, 'same'); 
coneLG2 = imfilter(coneL2, G, 'same');

% plot
vcNewGraphWin; plot(coneLG1); hold on; plot(coneLG2, 'r');
title('Cone absorption for a horizontal line in two scenes with eye move');

%% Compute error rate
%  This is the error rate from the ideal observer model
%  It uses only the information from one horizontal line of L cones
%  However, the accuracy can still be extremely high
errRate = 1/2 * exp(-1/4 * sum((coneLG1-coneLG2).^2 ./ ...
                sqrt(coneLG1 + coneLG2)));

%% Generate noise samples 
%  Init some parameters
nFrames = 1000; % number of samples to be generated

%  Randomly set eye-movement
sensor = sensorSet(sensor, 'exp time', 0.01);

sensor = ctInitEyeMovements(sensor, scene{1}, OIs{1}, 5*nFrames);
sensor = coneAbsorptions(sensor, OIs{1});
pSamples1 = double(sensorGet(sensor, 'photons'));


sensor = ctInitEyeMovements(sensor, scene{2}, OIs{2}, 5*nFrames);
sensor = coneAbsorptions(sensor, OIs{2});
pSamples2 = double(sensorGet(sensor, 'photons'));

%  Fit Gaussian
gMu1 = mean(pSamples1, 3); gMu2 = mean(pSamples2, 3);

%  Compute error rate


%% Do it by SVM
% Classification
svmOpts = '-s 0 -q';

nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
indx = randperm(5*nFrames);
refPhotons   = RGB2XWFormat(pSamples1);
szN = size(refPhotons, 1);
refPhotons = sum(reshape(refPhotons(:,indx), [szN, nFrames, 5]),3)';
matchPhotons = RGB2XWFormat(pSamples2);
matchPhotons = sum(reshape(matchPhotons(:,indx), [szN, nFrames, 5]),3)';

acc = svmClassifyAcc(cat(1,refPhotons, matchPhotons), ...
    labels, nFolds, 'svm', svmOpts);