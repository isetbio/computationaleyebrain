%% s_displayLetters
%
%
%
% (HJ) Dec, 2013

%% Init Parameters
if notDefined('letter'),     letter = 'v'; end
if notDefined('fontFamily'), fontFamily = 'Georgia'; end
if notDefined('fontSize'),   fontSize   = 13; end
if notDefined('vDist'),      vDist = 0.4:0.1:1; end
if notDefined('nFrames'),    nFrames = 1000; end

%% Create virtual display
fName = fullfile(ctRootPath,'ctData','Display Models','Dell LCD Stripe Pixels.mat');
vDisp  = ctDisplayLoad(fName);
vDisp.sStimulus.oSample = 6;

%% Create Scene
scene = cell(2, 1);
OIs   = cell(2, 1);
stim = stimCreate;
stim = stimSet(stim,'letter',letter);
stim = stimSet(stim,'font Family',fontFamily);
stim = stimSet(stim,'font Size',fontSize);
dpi = 300;
stim = stimSet(stim,'font dpi',dpi);

%% Create Sensor
emDuration = 0.01;
expTime    = 0.05;
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'exp time', expTime);


%% Create Human Lens
%  Create a typical human lens
%  Notice here the optical image haven't been computed yet
%  OI will be computed in together with cone samples
wave   = 380 : 4 : 780;
wvf    = wvfCreate('wave',wave);
pupilDiameterMm = 3;
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

if exist('defocus', 'var')
    wvf    = wvfSet(wvf, 'calc observer focus correction', defocus);
end
wvf    = wvfComputePSF(wvf, false);
oi     = wvf2oi(wvf,'shift invariant', false);

%% Loop over viewing distance
absorptions = cell(2, 1);
acc = zeros(length(vDist), 1); % classification accuracy
err = zeros(length(vDist), 1); % classification std error
svmOpts = '-s 0 -q';
oSamples = vDisplayGet(vDisp, 'oSample');

for indxDist = 1: length(vDist)
    fprintf('viewing distance: %.2f\n', vDist(indxDist));
    
    dist = vDist(indxDist);
    % Enable adjusting the filters here, rather than in ctCreateFontScene
    scene{1} = ctCreateFontScene(vDisp, dist, stim, true);  
    scene{2} = ctCreateFontScene(vDisp, dist, stim, false);
    
    % Compute cone absorptions
    fov = atand(tand(0.3371/2) * 0.6 / dist) * 2;
    for ii = 1 : 2
        % compute field of view
        % fov = 2*atan(ncols /oSamples * fontSize * 0.3514 / 1000 /dist/2);
        scene{ii} = sceneSet(scene{ii}, 'h fov', fov);
        vcAddAndSelectObject('scene', scene{ii});
        OIs{ii} = oiCompute(scene{ii}, oi);
    end
    
    sensor = sensorSetSizeToFOV(sensor, fov, scene{1}, OIs{1});

    sensor = sensorSet(sensor, 'exp time', emDuration);
    
    params.center   = [0,0];
    params.Sigma    = 1e-4 *[0.3280 0.0035; 0.0035 0.4873]*emDuration*4000;
    emPerExposure = round(expTime / emDuration);
    params.nSamples = nFrames * emPerExposure;
    params.fov      = fov;
    
    % Set up the eye movement properties
    sensor = emInit('fixation gaussian', sensor, params);
    
    for ii = 1 : 2
        sensor = coneAbsorptions(sensor, OIs{ii}, 2);
        absorptions{ii} = double(sensorGet(sensor, 'photons'));
    end
    
    absorptions{2} = absorptions{2} * mean(absorptions{1}(:)) / mean(absorptions{2}(:));
    % Do classification
    nFolds = 10;
    labels = [ones(nFrames,1); -1*ones(nFrames,1)];
    refPhotons   = RGB2XWFormat(absorptions{1});
    %refPhotons = cumsum(refPhotons, 2);
    %refPhotons = refPhotons(:, emPerExposure + 1:end) - ...
    %    refPhotons(:, 1:end-emPerExposure);
    sz = size(refPhotons);
    refPhotons = sum(reshape(refPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);
    
    matchPhotons = RGB2XWFormat(absorptions{2});
%     matchPhotons = cumsum(matchPhotons, 2);
%     matchPhotons = matchPhotons(:, emPerExposure + 1:end) - ...
%         matchPhotons(:, 1:end-emPerExposure);
    matchPhotons = sum(reshape(matchPhotons, [sz(1) sz(2)/emPerExposure emPerExposure]), 3);
    
    accuracy = svmClassifyAcc(cat(1,refPhotons', matchPhotons'), ...
        labels, nFolds, 'linear', svmOpts);
    err(indxDist) = accuracy(2);
    acc(indxDist) = accuracy(1);
    fprintf('acc: %.2f\n', acc(indxDist));
end
