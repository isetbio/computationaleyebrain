%% v_ptbIsetIsomerizations
%
% Confirm that when we start with the same irradiance function we get
% the same number of estimated isomerizations in PTB and ISET.
%
%
% DHB/BW/HJ Copyright ISETBIO Team, 2013


%%
s_initISET

%% Create a radiance image
scene = sceneCreate('uniform ee');    % Equal energy
scene = sceneSet(scene,'name','Equal energy uniform field');
scene = sceneSet(scene,'fov',20);     % Big field required
wave  = sceneGet(scene,'wave');

sz = sceneGet(scene,'size');
rect = [sz(2)/2,sz(1)/2,5,5];
roiLocs = ieRoi2Locs(rect);
radianceData = plotScene(scene,'radiance energy roi',roiLocs);

title(sprintf(sceneGet(scene,'name')));

radiance = mean(radianceData.energy(:));

%% Compute the irradiance in ISETBIO
oi     = oiCreate('human');
optics = oiGet(oi,'optics');

% Turn off relative illumination (cos4th)
optics = opticsSet(optics,'off axis method','skip');

% Turn off the OTF
optics = opticsSet(optics,'otf method','skip otf');
oi     = oiSet(oi,'optics',optics);
oi     = oiCompute(oi,scene);
vcAddAndSelectObject(oi); oiWindow;

%% Extract ISETBIO irradiance from middle of data

sz         = oiGet(oi,'size');
rect       = [sz(2)/2,sz(1)/2,2,2];
roiLocs    = ieRoi2Locs(rect);
irData = plotOI(oi,'irradiance energy roi',roiLocs);
title(sprintf(oiGet(oi,'name')));
irradiance = mean(irData.y(:));


%% PTB Calculation

optics = oiGet(oi,'optics');
pupilDiameterMm    = opticsGet(optics,'pupil diameter','mm');
focalLengthMm      = opticsGet(optics,'focal length','mm');
integrationTimeSec = 0.05;
m = opticsGet(optics,'magnification',sceneGet(scene,'distance'));

% The v_ptbIsetIrradiance and this differ by about 1 part in 100.  Figure
% out the problem some day.
% [isoPerCone,pupilDiamMm,photoreceptors,irradianceWattsPerM2]
macularPigmentOffset = 0;
[isoPerCone, ~, ptbPhotoreceptors, ptbIrradiance] = ...
    ptbConeIsomerizationsFromRadiance(radiance(:), wave(:),...
    pupilDiameterMm, focalLengthMm, integrationTimeSec,macularPigmentOffset);

% ptb effective absorbtance
ptbCones = ptbPhotoreceptors.isomerizationAbsorptance'; % Appropriate for quanta

%% Compare the irradiances 

% They differ - this time more than I would like.
vcNewGraphWin; plot(wave,ptbIrradiance(:),'ro',wave,irradiance(:),'ko');
set(gca,'ylim',[0 max(ptbIrradiance(:))*1.5]);
legend('PTB','ISETBIO')
xlabel('Wave (nm)'); ylabel('Irradiance (q/s/nm/m^2)')
title('Without magnification correction');

% This one accounts for the magnification difference, so they should be
% really close.  The magnification difference results from how Peter
% Catrysse implemented the radiance to irradiance calculation, with a small
% formula difference that accounts for the conversion of visual angle to
% surface of the retina.
%
% But still, they differ by about 1/100.  This may be because of some
% spatial interpolation and blurring that we haven't accounted for in the
% oiCompute/optics.
%
vcNewGraphWin; 
plot(wave,ptbIrradiance(:)/(1+abs(m))^2,'ro',wave,irradiance(:),'ko');
set(gca,'ylim',[0 max(ptbIrradiance(:))*1.5]);
xlabel('Wave (nm)'); ylabel('Irradiance (q/s/nm/m^2)')
legend('PTB','ISETBIO')
title('Magnification corrected comparison');

%%  ISETBIO sensor absorptions
%
sensor    = sensorCreate('human');
isetCones = sensorGet(sensor,'spectral qe');
isetCones = isetCones(:,2:4);

%% Compare PTB sensor spectral responses with ISETBIO
vcNewGraphWin; plot(wave, isetCones);
hold on; plot(wave, ptbCones, '--');

vcNewGraphWin; plot(wave,ptbCones);
plot(ptbCones(:),isetCones(:),'o');
hold on; plot([0 1],[0 1], '--');
xlabel('PTB cones');
ylabel('ISET cones');


%% End
