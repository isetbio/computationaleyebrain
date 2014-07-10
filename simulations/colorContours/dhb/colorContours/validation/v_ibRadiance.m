%% An ISETBIO imperfection - radiance to irradiance across wavelength
%
% Note that a small amount of error (one-tenth of one percent) creeps into
% the calculation when we transform from radiance to irradiance.  The
% perfectly uniform scene energy as a function of wavelength has a very
% slight energy difference in the irradiance.  This must have to do with
% the numerical precision of the otf calculation in opticsOTF, the main
% loop, around line 104.
%
% Can someone get me some more numerical precision?
%
% I don't think it has any practical significance.  But, I do hate it.
%
% BW Copyright ISETBIO Team, 2013

%%
s_initISET

%% Create an radiance image
scene = sceneCreate('uniform ee');    % Equal energy
scene = sceneSet(scene,'fov',3);      % Two deg field of view
wave = sceneGet(scene,'wave');

sz = sceneGet(scene,'size');
rect = [sz(2)/2,sz(1)/2,5,5];
roiLocs = ieRoi2Locs(rect);
radiance = vcGetROIData(scene,roiLocs,'energy');
radiance = mean(radiance);
err = (max(radiance(:)) - min(radiance(:)))/mean(radiance(:));
fprintf('Percent error is %.3f\n',err*100)

vcNewGraphWin; 
plot(wave,radiance); 
set(gca,'ylim',[min(radiance(:))/2 max(radiance(:))*1.5]); grid on

%% Compute the irradiance in ISETBIO
oi    = oiCreate;
oi    = oiCompute(oi,scene);

%% Note that the irradiance is not completely equal energy.
sz      = oiGet(oi,'size');
rect    = [sz(2)/2,sz(1)/2,5,5];
roiLocs = ieRoi2Locs(rect);
irradiance = vcGetROIData(oi,roiLocs,'energy');
irradiance = mean(irradiance);
err = (max(irradiance(:)) - min(irradiance(:)))/mean(irradiance(:));
fprintf('Percent error is %.3f\n',err*100)

vcNewGraphWin; 
plot(wave,irradiance); 
set(gca,'ylim',[min(irradiance(:))/2 max(irradiance(:))*1.5]); grid on

%%
vcAddAndSelectObject(oi); oiWindow;

%%