%% v_ptbIsetIsomerizations
%
% Confirm that when we start with the same irradiance function we get
% the same number of estimated isomerizations
%
%
% DHB/BW Copyright ISETBIO Team, 2013


%%
s_initISET

%% Create an radiance image
scene = sceneCreate('uniform ee');    % Equal energy
scene = sceneSet(scene,'fov',2);      % Two deg field of view
wave = sceneGet(scene,'wave');

sz = oiGet(scene,'size');
rect = [sz(2)/2,sz(1)/2,5,5];
roiLocs = ieRoi2Locs(rect);
radiance = vcGetROIData(scene,roiLocs,'energy');
radiance = mean(radiance);

vcNewGraphWin; 
plot(wave,radiance); 
xlabel('Wavelength'); ylabel('Radiance (energy)')
set(gca,'ylim',[min(radiance(:))/2 max(radiance(:))*1.5]); grid on

%% Compute the irradiance in ISETBIO
oi    = oiCreate('human');
oi    = oiCompute(oi,scene);
vcAddAndSelectObject(oi); oiWindow;

%%  ISETBIO sensor absorptions
sensor = sensorCreate('human');
sensor = sensorSet(sensor,'exp time',0.050);

sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor); sensorWindow('scale',1);

%% Extract ISETBIO irradiance
sz   = oiGet(oi,'size');
rect = [sz(2)/2,sz(1)/2,2,2];
roiLocs    = ieRoi2Locs(rect);
irradiance = vcGetROIData(oi,roiLocs,'energy');
irradiance = mean(irradiance);
vcNewGraphWin; plot(wave,irradiance); grid on

%% PTB Calculation

optics = oiGet(oi,'optics');
pupilDiameterMm    = opticsGet(optics,'pupil diameter','mm');
focalLengthMm      = opticsGet(optics,'focal length','mm');
integrationTimeSec = sensorGet(sensor,'exp time','sec');

% The v_ptbIsetIrradiance has a much closer value.  Figure out the problem.
% [isoPerCone,pupilDiamMm,photoreceptors,irradianceWattsPerM2]
[~, ~, ptbPhotoreceptors, ptbIrradiance] = ...
    ptbConeIsomerizationsFromRadiance(radiance(:), wave(:),...
    pupilDiameterMm, focalLengthMm, integrationTimeSec,0);

%% Compare the irradiances 

% They differ - this time more than I would like.
vcNewGraphWin; plot(wave,ptbIrradiance(:),'ro',wave,irradiance(:),'ko');
set(gca,'ylim',[0 max(ptbIrradiance(:))*1.5],'xlim',[0 max(ptbIrradiance(:))*1.5]);
legend('PTB','ISETBIO')

%% Then go on to check the cone isomerizations

%% End
