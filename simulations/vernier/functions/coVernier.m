function [acc, w] = coVernier(varargin)
% compute discrimination probability for vernier acuity
%
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) ISETBIO TEAM, 2016

% Parse inputs
p = parseInputs(varargin{:});

d = p.Results.display;            % display structure
imgFov = p.Results.imgFov;        % image field of view
sensorFov = p.Results.spatialInt; % spatial integration
nFrames = p.Results.nFrames;      % number of frames in training
vDist = p.Results.vDist;          % viewing distance
expTime = p.Results.expTime;      % exposure time
emFlag = p.Results.emFlag;        % eye movement indicators
sensor = p.Results.sensor;

ppi = displayGet(d, 'ppi');       % display resolution
imgSz  = round(tand(imgFov)*vDist*39.37*ppi);  % number of pixels in image
imgFov = atand(imgSz/ppi/39.37/vDist);         % Actual fov

% verneir scene parameters
params.display = d;
params.sceneSz = [imgSz imgSz];
params.barWidth = 1;
params.barColor = p.Results.barColor;
params.bgColor = 0.5;

% Create Scene
scene = cell(2, 1);
params.offset = 0; scene{1} = sceneCreate('vernier', 'display', params);
params.offset = 1; scene{2} = sceneCreate('vernier', 'display', params);

% set scene fov
for ii = 1 : 2
    scene{ii} = sceneSet(scene{ii}, 'h fov', imgFov);
    scene{ii} = sceneSet(scene{ii}, 'distance', vDist);
end

% Create Human Lens
oi = oiCreate('wvf human');

% Compute optical image
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);

% Create Sensor
sensor = sensorSet(sensor, 'exp time', expTime);
sensor = sensorSetSizeToFOV(sensor, sensorFov(1), scene{1}, OIs{1});
sz = sensorGet(sensor, 'size');

if all(emFlag == 0)  % no eye movement
    sensor = sensorSet(sensor, 'sensor positions', zeros(nFrames, 2));
    sensor = coneAbsorptions(sensor, OIs{1});
    pSamples1 = sensorGet(sensor, 'photons');
    
    sensor = coneAbsorptions(sensor, OIs{2});
    pSamples2 = sensorGet(sensor, 'photons');
else
    % Set exposure time to 1 ms
    emDuration = 0.001;
    emPerExposure = expTime / emDuration;
    sensor = sensorSet(sensor, 'exp time', emDuration);
    
    % Generate eyemovement
    params.emFlag = emFlag;
    params.nSamples = nFrames * emPerExposure;
    sensor = eyemoveInit(sensor, params);
    
    % Compute the cone absopritons
    sensor = coneAbsorptions(sensor, OIs{1});
    
    % Store the photon samples and add photons in one exposure time
    pSamples1 = sensorGet(sensor, 'photons');
    pSamples1 = sum(reshape(pSamples1, [sz nFrames emPerExposure]), 4);
    
    % Compute cone absorptions for the second stimulus and store photon
    % absorptions
    sensor = coneAbsorptions(sensor, OIs{2});
    pSamples2 = sensorGet(sensor, 'photons');
    pSamples2 = sum(reshape(pSamples2, [sz nFrames emPerExposure]), 4);
end

% prepare data for SVM linear classification
nFolds = 5;
labels = [ones(nFrames,1); -1*ones(nFrames,1)];
data = cat(1, RGB2XWFormat(pSamples1)', RGB2XWFormat(pSamples2)');

% Normalize data
data = bsxfun(@rdivide, bsxfun(@minus, data, mean(data)), std(data));

% Choose parameter C in svm and compute cross-validation error
opts = sprintf('-s 2 -q -C -v %d', nFolds);
res = train(labels, sparse(data), opts);
acc = res(2);

% Get weights of svm linear classifier
if nargout > 1
    opts = sprintf('-s 2 -q -c %e', res(1));
    svmStruct = train(labels, sparse(data), opts);
    w = reshape(svmStruct.w, sz);
end

end

function p = parseInputs(varargin)
% Parse input parameters in name-value pairs
%
% Support parameters:
%   'display'    - display structure
%   'imgFov'     - image field of view (degrees)
%   'spatialInt' - spatial integration size (degrees)
%   'vDist'      - viewing distance (meters)
%   'expTime'    - exposure time of human cone
%   'barColor'   - color of the bar, can be scalar or 3 element vector
%
% HJ, ISETBIO TEAM, 2016

p = inputParser;

p.addParameter('display', displayCreate('LCD-Apple', 'dpi', 400));
p.addParameter('imgFov', 0.5, @isnumeric);
p.addParameter('spatialInt', 0.2, @isnumeric);
p.addParameter('expTime', 0.05, @isnumeric);
p.addParameter('nFrames', 5000, @isnumeric);
p.addParameter('vDist', 1.0, @isnumeric);
p.addParameter('barColor', 0.99, @isnumeric);
p.addParameter('emFlag', [1 1 1], @isnumeric);
p.addParameter('sensor', sensorCreate('human'));

p.parse(varargin{:});

end