%% s_rgcSpikesTS
%
%
%  (HJ) April, 2014

%% Init
totTime = 50; % 1000 ms
expTime = 0.001; % 1 ms
wave = 400:10:700; % wavelengths samples
pupilDiameterMm = 3; % pupil diameter
fov = 1; % scene field of view
nFrames = 2; % This number is too small, but it's how much we can have
imageSz = 256;
display = displayCreate('LCD-Apple');

%% Load real human data
dataDir = fullfile(frontendRootPath, 'data', 'cone-connection');
data = load(fullfile(dataDir, 'rgcConeConnections.mat'));
indx = (data.rgcType == 3);
data.rgcType = data.rgcType(indx);
data.coneWeights = data.coneWeights(:, indx);

%% Load cone temperal inpulse response
coneIR = load(fullfile(frontendRootPath, 'data', 'coneIR.mat'));
coneIR.tIRF = coneIR.tIRF * 50;

%% Create first scene
%  This scene is just used to estimate sensor size
I = zeros(imageSz);
scene = sceneFromFile(I, 'rgb', [], display);
scene = sceneSet(scene, 'h fov', fov);

%% Create human oi
% generate human optics structure
wvf    = wvfCreate('wave',wave);
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

oi = oiCompute(scene, oi);

%% Create sensor for human cones
density = sum(data.coneType-2)/length(data.coneType);
pparams.rgbDensities = [0 1-density density 0];
sensor = sensorCreate('human', [], pparams);
sensor = sensorSet(sensor, 'exp time', expTime);
sensor = sensorSetSizeToFOV(sensor, fov, scene, oi);

% Get cone positions
xy = sensorGet(sensor, 'xy');
coneType = sensorGet(sensor, 'cone type');
coneType = coneType(:);
indxL = (data.coneType == 2); xyL = data.conePos(indxL, :);
indxM = (data.coneType == 3); xyM = data.conePos(indxM, :);

% Set eye movement parameters
params.center   = [0,0];
params.Sigma    = 1e-4 *[0.3280 0.0; 0.0 0.0];

params.nSamples = nFrames;
params.fov      = sensorGet(sensor,'fov',scene, oi);

%% Scene -> oi -> Cone Absorptions -> rgc volts -> spikes
volts = zeros(length(data.coneType), 1);
rgcVolts = zeros(length(data.rgcType), totTime);
curV = zeros(length(data.rgcType), nFrames);
spikes = zeros([size(rgcVolts) nFrames]);

for ts = 1: totTime
    fprintf('Computing %d\n',ts);
    % Renew scene
    % I(:, 1:round(ts/totTime*imageSz)) = 1;
    I(:) = 1;
    scene = sceneFromFile(I, 'rgb', [], display);
    scene = sceneSet(scene, 'h fov', fov);
    
    % Renew eye movement
    sensor = emInit('fixation gaussian', sensor, params);
    
    % Compute cone absorptions
    oi = oiCompute(scene, oi);
    sensor = coneAbsorptions(sensor, oi, 0);
    
    % Get cone absorptions and interpolate
    v = double(sensorGet(sensor, 'volts'));
    v = reshape(v, [size(xy, 1) size(v, 3)]);
    
    for ii = 1 : size(v, 2)
        % interpolate for each cone type
        % in EJ's data, we only have L and M cones
        indx = (coneType == 2);
        F = scatteredInterpolant(xy(indx, 1), xy(indx, 2), ...
                    v(indx,ii), 'linear', 'nearest');
        volts(indxL) = F(xyL(:,1), xyL(:,2));
        indx = (coneType == 3);
        F = scatteredInterpolant(xy(indx, 1), xy(indx, 2), ...
                    v(indx,ii), 'linear', 'nearest');
        volts(indxM) = F(xyM(:,1), xyM(:,2));
                
        % Compute volts after rgc summing up
        % Logically, we should first compute all time sequence data,
        % convolve it with cone inpulse response and then do the spatial
        % summing up. Considering associative rule for matrix
        % multiplication, we could first do spatial summing and then
        % temporal convolution
        rgcVolts(:, ts) = data.coneWeights' * volts; % spatial
        
        startT = max(1, ts - length(coneIR.tIRF) + 1);
        curV(:, ii) = curV(:, ii) + rgcVolts(:, ts:-1:startT) * ...
                            coneIR.tIRF(1:ts-startT+1,1);
    end
    % We need to figure out how to set parameter for sigmoidal function
    % alpha = -log(1/0.95-1)/(quantile(curV(:), 0.9)-mean(curV(:)));
    % alpha = min(alpha, 100);
    threshold = sigmf(curV, [10 3]);
    accIndx = rand(length(data.rgcType), nFrames) < threshold;
    curV(accIndx) = 0;
    spikes(:,ts,:) = accIndx;
end

% spikes = sparse(spikes);

%% Visualize