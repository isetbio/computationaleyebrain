%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
if notDefined('ppi'), ppi = 500; end            % points per inch
if notDefined('imgFov'), imgFov = [6 6]/60; end % visual angle (degree)
if notDefined('nFrames'), nFrames = 1000; end   % Number of samples

vDist  = 1.0;                                   % viewing distance (meter)
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/39.37/vDist);     % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
scene = cell(2, 1);
img = ones(imgSz)*0.5;            % Init to black
img(:, round(imgSz(2)/2)) = .99;  % Draw vertical straight line in middle
scene{1} = sceneFromFile(img, 'rgb', [], display); % create scene

img(1:round(imgSz(1)/2), :) = circshift(img(1:round(imgSz(1)/2), :),[0 1]);
%img = circshift(img, [0 1]);
scene{2} = sceneFromFile(img, 'rgb', [], display);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Show radiance image (scene)
% vcAddAndSelectObject('scene', scene{1});
% vcAddAndSelectObject('scene', scene{2}); sceneWindow;

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

% Show irradiance (optical image) 
%vcAddAndSelectObject('oi', OIs{1}); oiWindow;
%vcAddAndSelectObject('oi', OIs{2}); oiWindow;

%% Create Sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', 0.05);

%% Compute cone absorption (noise free)
%  Compute for scene 1
sensor = sensorSetSizeToFOV(sensor, imgFov, scene{1}, OIs{1});
sensor = sensorComputeNoiseFree(sensor, OIs{1});
p1 = double(sensorGet(sensor, 'photons')); 
% show cone absorptions
% vcAddAndSelectObject(sensor); sensorWindow;

% Compute for scene 2
sensor = sensorSetSizeToFOV(sensor, imgFov, scene{2}, OIs{2});
sensor = sensorComputeNoiseFree(sensor, OIs{2});
p2 = double(sensorGet(sensor, 'photons'));
% show cone absorptions
% vcAddAndSelectObject(sensor); sensorWindow;

coneType = sensorGet(sensor, 'cone type');
coneL1 = zeros(size(p1, 2), 1); coneL2 = coneL1;
for curCol = 1 : size(p1, 2)
    coneL1(curCol) = mean(p1(coneType(:,curCol) == 2,curCol));
    coneL2(curCol) = mean(p2(coneType(:,curCol) == 2,curCol));
end
% Plot
% vcNewGraphWin; plot(coneL1); hold on; plot(coneL2, 'r');
% title('Cone absorption for a horizontal line in two scenes');

%% Add some Gaussian blur (for eye-movement)
% G = fspecial('gaussian', [16 1], 8); % eye movement is around 8 cones wide
% coneLG1 = imfilter(coneL1, G, 'same'); 
% coneLG2 = imfilter(coneL2, G, 'same');

% plot
% vcNewGraphWin; plot(coneLG1); hold on; plot(coneLG2, 'r');
% title('Cone absorption for a horizontal line in two scenes with eye move');

%% Compute error rate
%  This is the error rate from the ideal observer model
%  It uses only the information from one horizontal line of L cones
%  However, the accuracy can still be extremely high
% errRate = 1/2 * exp(-1/4 * sum((coneLG1-coneLG2).^2 ./ ...
%                 sqrt(coneLG1 + coneLG2)));

%% Generate noise samples
%  Set exposure time to 10 ms
emDuration = 0.01;
sensor = sensorSet(sensor, 'exp time', emDuration);

%  Init some parameters
%  Simga is etimated from Jon's data fitted by Brownian motion
params.center   = [0,0];
params.Sigma    = 1e-4 *[0.3280 0.0035; 0.0035 0.4873]*emDuration*1000;
params.nSamples = round(0.05 / emDuration)*nFrames;
params.fov      = sensorGet(sensor,'fov',scene{1},oi);

% Set up the eye movement properties
sensor = emInit('fixation',sensor,params);

% Compute the cone absopritons
sensor = coneAbsorptions(sensor, OIs{1});

% Store the photon samples
pSamples1 = double(sensorGet(sensor, 'photons'));

% Compute cone absorptions for the second stimulus and store photon
% absorptions
sensor = coneAbsorptions(sensor, OIs{2});
pSamples2 = double(sensorGet(sensor, 'photons'));

%  Fit Gaussian mean for ideal observer calculation.
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


%%