%% s_VernierAcuity
%    Testing if people can see the difference between two cases:
%      1) A straight line
%      2) Two straight lines with 1 pixel apart
%    Vernier Acuity shows the positional acuity is around 6 sec of arc
%
% (HJ) Jan, 2014

%% Init Parameters
ppi = 100;          % Display property: points per inch
vDist = 1:0.1:2.0;  % Viewing distance in meters
nFrames = 100;      % Number of samples for each case
imgSz = [12 12];    % Row / columns of image in pixels

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', ppi);

%% Create Scene
scene = cell(2, 1);
img = zeros(imgSz);            % Init to black
img(:, round(imgSz(2)/2)) = .99; % Draw vertical straight line in the middle
scene{1} = sceneFromFile(img, 'rgb', [], display); % create scene

%img(1:imgSz(1)/2, :) = circshift(img(1:imgSz(1)/2, :), [0 1]);
img = circshift(img, [0 2]);
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
OIs{1} = oiCompute(scene{1}, oi);
OIs{2} = oiCompute(scene{2}, oi);
% vcAddAndSelectObject('oi', OIs{1}); oiWindow;
% vcAddAndSelectObject('oi', OIs{2}); oiWindow;

%% Loop over viewing distance
absorptions = cell(2, 1);
acc = zeros(length(vDist), 1); % classification accuracy
err = zeros(length(vDist), 1); % classification std error
svmOpts = '-s 0 -q';

for indxDist = 1: length(vDist)
    % Compute cone absorptions
    dist = vDist(indxDist);
    for ii = 1 : 2
        % Update scene information
        scene{ii} = sceneSet(scene{ii}, 'distance', dist);
        ncols = sceneGet(scene{ii}, 'ncols');
    end
    
    sensor = sensorSetSizeToFOV(sensor, fov, scene{1}, oi);
    sensor = sensorSet(sensor, 'exp time', 0.01); % 10 ms for each
    sensor = ctInitEyeMovements(sensor, scene{1}, oi, 5*nFrames);
    
    for ii = 1 : 2
        sensor = coneSamples(scene{ii}, nFrames, sensor, oi);
        absorptions{ii} = double(sensorGet(sensor, 'photons'));
    end
    
    % Do classification
    nFolds = 10;
    labels = [ones(nFrames,1); -1*ones(nFrames,1)];
    indx = randperm(5*nFrames);
    refPhotons   = RGB2XWFormat(absorptions{1});
    szN = size(refPhotons, 1);
    refPhotons = sum(reshape(refPhotons(:,indx), [szN, nFrames, 5]),3)';
    matchPhotons = RGB2XWFormat(absorptions{2});
    matchPhotons = sum(reshape(matchPhotons(:,indx), [szN, nFrames, 5]),3)';
    accuracy = svmClassifyAcc(cat(1,refPhotons, matchPhotons), ...
        labels, nFolds, 'svm', svmOpts);
    err(indxDist) = accuracy(2);
    acc(indxDist) = accuracy(1);
end

% sensor = sensorComputeNoiseFree(sensor, OIs{2});
% sensor = sensorSet(sensor, 'exp time', 0.05);
% sensor = sensorSetSizeToFOV(sensor, fov, scene{1}, OIs{1});
% sensor = sensorComputeNoiseFree(sensor, OIs{1});
% sensor = sensorCompute(sensor, OIs{2});
% p2 = sensorGet(sensor, 'photons');
% for curCol = 1 : 35
% tmpL2(curCol) = mean(p2(find(coneType(:,curCol) == 2),curCol));
% end;
% plot(tmpL); hold on; plot(tmpL2, 'r');