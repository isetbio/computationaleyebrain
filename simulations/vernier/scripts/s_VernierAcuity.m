%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
ieInit;

ppi = 1000;            % points per inch
imgFov = [.1 .1];      % image field of view
nFrames = 5000;        % Number of samples

vDist  = 1.0;                                   % viewing distance (meter)
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);   % number of pixels in image
imgFov = atand(max(imgSz)/ppi/39.37/vDist);     % Actual fov

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
%  create vernier scene
scene = cell(2, 1);

params.display = display;
params.sceneSz = imgSz;
params.offset  = 0;
params.barWidth = 1;
params.bgColor = 0.5;
scene{1} = sceneCreate('vernier', 'display', params);
params.offset = 1;
scene{2} = sceneCreate('vernier', 'display', params);

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
oi = oiCreate('wvf human');

% Compute optical image
% Actually, we could wait to compute it in coneSamples
% But, we compute it here to do sanity check
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

% Show irradiance (optical image) 
%vcAddAndSelectObject('oi', OIs{1}); oiWindow;
%vcAddAndSelectObject('oi', OIs{2}); oiWindow;

%% Create Sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', 0.05);

% Now sensor is created with color filer of normalized cone fundamentals.
% We want to replace it with the effective cone absorptions. Also, we will
% take lens transmittance, macular pigment transmittance into account.
%cone = coneCreate('human');
%effAbsorption = coneGet(cone, 'effective absorptance');
%sensor = sensorSet(sensor, 'filter spectra', effAbsorption);

%% Compute cone absorption (noise free)
%  Compute for scene 1
%sensor = sensorSetSizeToFOV(sensor, imgFov, scene{1}, OIs{1});
%sensor = sensorComputeNoiseFree(sensor, OIs{1});
%p1 = double(sensorGet(sensor, 'photons')); 
% % show cone absorptions
% % vcAddAndSelectObject(sensor); sensorWindow;
% 
% % Compute for scene 2
%sensor = sensorSetSizeToFOV(sensor, imgFov, scene{2}, OIs{2});
%sensor = sensorComputeNoiseFree(sensor, OIs{2});
%p2 = double(sensorGet(sensor, 'photons'));
% % show cone absorptions
% % vcAddAndSelectObject(sensor); sensorWindow;
% 
% coneType = sensorGet(sensor, 'cone type');
% coneL1 = zeros(size(p1, 2), 1); coneL2 = coneL1;
% for curCol = 1 : size(p1, 2)
%     coneL1(curCol) = mean(p1(coneType(:,curCol) == 2,curCol));
%     coneL2(curCol) = mean(p2(coneType(:,curCol) == 2,curCol));
% end
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
%  Set exposure time to 1 ms
expTime = sensorGet(sensor, 'exp time');
emDuration = 0.001;
sensor = sensorSet(sensor, 'exp time', emDuration);
sensor = sensorSetSizeToFOV(sensor, imgFov, scene{1}, OIs{1});

% Generate eyemovement
p.nSamples = nFrames + 50;
sensor = eyemoveInit(sensor, p);

% Compute the cone absopritons
sensor = coneAbsorptions(sensor, OIs{1}, 2);

% Store the photon samples
pSamples1 = sensorGet(sensor, 'photons');

% Compute cone absorptions for the second stimulus and store photon
% absorptions
sensor = coneAbsorptions(sensor, OIs{2}, 2);
pSamples2 = sensorGet(sensor, 'photons');


%% Do it by SVM
% Classification
nFolds = 10;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
refPhotons   = RGB2XWFormat(pSamples1);

% Add 50 samples to 1 to form the 50ms integration time
% We do use moving sum here
emPerExposure = 50;
refPhotons = cumsum(refPhotons, 2);
refPhotons = refPhotons(:, emPerExposure + 1:end) - ...
             refPhotons(:, 1:end-emPerExposure);
       
szN = size(refPhotons, 1);
matchPhotons = RGB2XWFormat(pSamples2);
matchPhotons = cumsum(matchPhotons, 2);
matchPhotons = matchPhotons(:, emPerExposure + 1:end) - ...
             matchPhotons(:, 1:end-emPerExposure);

acc = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
    labels, nFolds, 'linear');

fprintf('SVM acc:%f\n',acc(1));

%% Do it by KMean
% idx = kmeans(cat(1,refPhotons', matchPhotons'),2);
% idx = idx * 2 - 3;
% fprintf('kmeans accuracy: %f\n',sum(idx == labels) / length(labels));