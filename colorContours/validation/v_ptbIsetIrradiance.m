%% v_PTBISETBIOIrradiance
%
%
%  Plot comparison of iset and ptb irradiance, optionally
%
% PTB, conversion is pupilArea/(eyeLength^2). (pi /(1 + 4*fN^2*(1+abs(m))^2)
%
% [**] This plot is currently broken because backOiD is no longer
% computed at this level.
%if (runtimeParams.DO_SIM_PLOTS)

%% ISETBIO formulation

scene = sceneFromFile('backFile.png','rgb',[],'LCD-Apple');
vcAddAndSelectObject(scene); sceneWindow;

%% Create standard human polychromatic PSF using wavefront tools
% Make it an iset OI thingy.
pupilDiameterMm = 3;   
wave = 380:4:780;
wvf = wvfCreate('wave',wave);
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf = wvfSet(wvf,'zcoeffs',sample_mean);
wvf = wvfComputePSF(wvf);
oiD = wvf2oi(wvf,'shift invariant');
optics = oiGet(oiD,'optics');

oiD = oiCompute(oiD,scene);
vcAddAndSelectObject(oiD); oiWindow;

% The difference is because of the magnification difference between ISET
% and PTB.  See note in doOneSimulation.  Document it here better.
vcNewGraphWin; plot(uData.y,ptbIrradiance,'.')
identityLine;
grid on


%% PTB computes the ptbIrradiance w/m2
d = displayCreate('LCD-Apple');
backRGB = [0.5, 0.5, 0.5]';                    % Define background for experment in monitor RGB
backSpd = displayGet(d,'spd')*backRGB;
wave    = displayGet(d,'wave');
focalLengthMm = opticsGet(optics,'focal length','mm');
integrationTimeSec = 0.05;
[~,~,ptbPhotoreceptors,ptbIrradiance] = ptbConeIsomerizationsFromSpectra(backSpd,wave,...
    pupilDiameterMm,focalLengthMm,integrationTimeSec,0);
vcNewGraphWin; plot(wave,ptbIrradiance); grid on

%%
if (0)
    % Get background irradiance out of the optical image.
    %
    % [**] This currently works be using an ROI that was selected
    % by hand an then stored in a .mat file.  May want to
    % make this more programmatic.  But, we get this just
    temp = load('roiLocs');
    backUdata = plotOI(backOiD,'irradiance energy roi',temp.roiLocs);
    isetIrradianceWattsPerM2 = backUdata.y';
    
    % Make a new plot of PTB and iset irradiance.  Not quite
    % sure why we don't just add this to the window that
    % comes up in the call to PlotOI above.
    figure; hold on
    plot(staticComputedValues.wavelengthsNm,isetIrradianceWattsPerM2,'r');
    plot(staticComputedValues.wavelengthsNm,ptbAdjustedIrradianceWattsPerM2,'k');
    theRatio = isetIrradianceWattsPerM2 ./ ptbAdjustedIrradianceWattsPerM2;
end