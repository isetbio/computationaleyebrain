%% s_displayLetters
%
%
%
% (HJ) Dec, 2013

%% Init Parameters
if notDefined('letter'),     letter = 'g'; end
if notDefined('fontFamily'), fontFamily = 'georgia'; end
if notDefined('fontSize'),   fontSize   = 12; end
if notDefined('dpi'),        dpi = 100; end
if notDefined('vDist'),      vDist = 0.7:0.1:1.6; end
if notDefined('nFrames'),    nFrames = 4000; end
if notDefined('filter'),     filter = [2 3]; end % filter a = 0.2 or 0.3

%% Create virtual display
display = displayCreate('LCD-Apple');
display = displaySet(display, 'dpi', dpi);

%% Create Scene
scene = cell(2, 1);
for ii = 1 : 2
    fName     = sprintf('letter_%s_%d_Filter_%d.png',...
                        letter, fontSize, filter(ii));
    scene{ii} = sceneFromFile(fName, 'rgb', [], display);
end

%% Create Sensor
if exist('humanConeDensities', 'var')
    pparams.humanConeDensities = humanConeDensities;
else
    pparams = [];
end

sensor = sensorCreate('human', [], pparams);
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
wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

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
    
    vcAddAndSelectObject('scene', scene{1});
    OIs{1} = oiCompute(scene{1}, oi);
    vcAddAndSelectObject('scene', scene{2});
    OIs{2} = oiCompute(scene{2}, oi);
    
    sensor = sensorSetSizeToFOV(sensor, fov, scene{1}, OIs{1});
    expTime = sensorGet(sensor, 'exp time');
    emDuration = 0.001;
    sensor = sensorSet(sensor, 'exp time', emDuration);
    
    params.center   = [0,0];
    params.Sigma    = 1e-4 *[0.3280 0.0035; 0.0035 0.4873]*emDuration*1000;
    emPerExposure = round(expTime / emDuration);
    params.nSamples = nFrames + emPerExposure;
    params.fov      = sensorGet(sensor,'fov',scene{1},OIs{1});
    
    % Set up the eye movement properties
    sensor = emInit('fixation gaussian', sensor, params);
    
    for ii = 1 : 2
        sensor = coneAbsorptions(sensor, OIs{ii}, 2);
        absorptions{ii} = double(sensorGet(sensor, 'photons'));
    end
    
    % Do classification
    nFolds = 10;
    labels = [ones(nFrames,1); -1*ones(nFrames,1)];
    refPhotons   = RGB2XWFormat(absorptions{1});
    refPhotons = cumsum(refPhotons, 2);
    refPhotons = refPhotons(:, emPerExposure + 1:end) - ...
        refPhotons(:, 1:end-emPerExposure);
    
    matchPhotons = RGB2XWFormat(absorptions{2});
    matchPhotons = cumsum(matchPhotons, 2);
    matchPhotons = matchPhotons(:, emPerExposure + 1:end) - ...
        matchPhotons(:, 1:end-emPerExposure);
    
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'svm', svmOpts);
    err(indxDist) = accuracy(2);
    acc(indxDist) = accuracy(1);
end
