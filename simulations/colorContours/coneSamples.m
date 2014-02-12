function [sensor, oi, xy, coneType] = coneSamples(scene, nSamples, ...
                                                        sensor, oi, cbType)
%% function coneSamples(scene, nSamples, [sensor], [oi], [cbType])
%    Compute cone samples for scene
%
%  Inputs:
%    scene        - ISET scene structure 
%    nSamples     - number of samples to be generated 
%    sensor       - eye sensor to be used, if empty, this program will 
%                   create a standard one
%    oi           - Optical image data, if empty, 
%    cbType       - colorblind type
%                   0 - Normal; 1 - Protan; 2 - Deuteran; 3 - Tritan
%
%  Outputs:
%    sensor       - sensor with cone absorption data computed
%    oi           - optical image
%    xy           - cone xy positions, see sensorCreateMosaic
%    coneType     - cone types of K, L, M, S, see sensorCreateMosaic
%
%  Example:
%    sensor = ctConeSamples(scene, 10)
%
%  See also:
%    coneAbsorptions, sensorCreateMosaic
%
%  (HJ) Sep, 2013 

%% Check Inputs
%  Check number of inputs
if nargin < 1, error('Scene structure required'); end
if nargin < 2, nSamples = 1000; end
if nargin < 3, sensor = []; end
if nargin < 4, oi = []; end
if nargin < 5, cbType = 0; end

%  Parse cbType
switch cbType
    case 0
        coneDensity = [.1 .6 .2 .1];
    case 1
        coneDensity = [.1 .0 .8 .1];
    case 2
        coneDensity = [.1 .8 .0 .1];
    case 3
        coneDensity = [.1 .7 .2 .0];
    otherwise
        error('Unknown cone type');
end

%% Create OI if not given
% The data were collected by Thibos and are described in the
% wvfLoadThibosVirtualEyes function and reference therein.
if isempty(oi)
    wave = 400:10:700; wave = wave(:);
    pupilMM = 3;
    
    % Load human wvf
    zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
    wvfP = wvfCreate('wave', wave, ...
                     'zcoeffs', zCoefs, ...
                     'name',sprintf('human-%d',pupilMM));
    wvfP = wvfComputePSF(wvfP);
    
    oi = wvf2oi(wvfP,'human');
    oi = oiSet(oi,'name','Human WVF 3mm');
end

%% Sensor Setup
%  Create a sensor if not given
if isempty(sensor)
    sensor = sensorCreate('human'); 
    sensor = sensorSet(sensor,'exp time',0.05);
    [sensor,xy,coneType] = sensorCreateConeMosaic(sensor, [], coneDensity);
    sensor = sensorSetSizeToFOV(sensor, sceneGet(scene,'hfov'),scene,oi);
end

%  sensor for Thibos calculation
oi     = oiCompute(oi,scene);
sensor = sensorCompute(sensor,oi);

%% Compute sensor image samples
%  If needed, eye-movement should be set up here
%  Now, we set it with no eye-movement
% sensor = sensorSet(sensor,'movement positions', [0 0]);
% sensor =  sensorSet(sensor,'frames per position', nSamples);
% vcAddAndSelectObject('scene',scene);
% vcAddAndSelectObject('oi', oi);

%  Compute photon cone absorptions
sensor = coneAbsorptions(sensor, oi);

end