%% s_displayLetters
%
%
%
% (HJ) Dec, 2013

%% Init Parameters
if notDefined('letter'),     letter = ['c','o']; end
if notDefined('fontFamily'), fontFamily = 'georgia'; end
if notDefined('fontSize'),   fontSize   = 9; end
if notDefined('dpi'),        dpi = 100; end
if notDefined('vDist'),      vDist = 1:0.1:2.0; end
if notDefined('nFrames'),    nFrames = 500; end

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', dpi);

%% Create Scene
scene = cell(2, 1);
for ii = 1 : 2
    fName     = sprintf('letter_%s_%d.png', letter(ii), fontSize);
    scene{ii} = sceneFromFile(fName, 'rgb', [], display);
end

%% Create Sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', 0.05);

%% Create Human Lens
%  Create a typical human lens
%  Notice here the optical image haven't been computed yet
%  OI will be computed in together with cone samples
wave   = 380 : 4 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);
wvf    = wvfComputePSF(wvf);
oi     = wvf2oi(wvf,'shift invariant');

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
        
        % compute field of view
        fov = 2*atan(ncols * fontSize * 0.3514 / 1000 /dist/2);
        scene{ii} = sceneSet(scene{ii}, 'h fov', fov);
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
