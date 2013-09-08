%% v_PTBISETBIOIrradiance
%
% Compare iset and ptb irradiance calculations
%
%
% PTB, conversion is pupilArea/(eyeLength^2). (pi /(1 + 4*fN^2*(1+abs(m))^2)
%
% [**] This plot is currently broken because backOiD is no longer
% computed at this level.
%if (runtimeParams.DO_SIM_PLOTS)

%% 
s_initISET

monitorName = 'LCD-Apple.mat';
wave   = 380:4:780;
bLevel = 0.5;    % Linear value of display primary

%% Create human optics with standard wave front aberrations

wvf    = wvfCreate('wave',wave);

pupilDiameterMm = 3;   
sample_mean = wvfLoadThibosVirtualEyes(pupilDiameterMm);
wvf    = wvfSet(wvf,'zcoeffs',sample_mean);

wvf    = wvfComputePSF(wvf);
oiB    = wvf2oi(wvf,'shift invariant');

%% ISETBIO formulation

% Write an image file with a gray background.  We want it set to one half
% the max primary intensity.  To do this, we need to read the DAC for the
% display calculate the right level.

d       = displayCreate(monitorName);
gTable  = displayGet(d,'gamma table');  % plot(gTable)
igTable = ieLUTInvert(gTable);        % Maps linear values to DAC
bDAC    = ieLUTLinear([0.5,0.5,0.5],igTable)/size(gTable,1);

bImage = ones(128,128,3);
for ii=1:3
    bImage(:,:,ii) = bImage(:,:,ii)*bDAC(ii);
end

bFile = fullfile(isetbioRootPath,'tmp','bFile.png');
imwrite(bImage,bFile);

% Build the scene from the image file
bScene = sceneFromFile(bFile,'rgb',[],monitorName,wave);
vcAddAndSelectObject(bScene); sceneWindow;

% % Create a scene with a background of [0.5,0.5,0.5]
% scene = sceneFromFile('backFile.png','rgb',[],'LCD-Apple');
% vcAddAndSelectObject(scene); sceneWindow;

%% ISETBIO path for creating an irradiance

oiB    = oiCompute(oiB,bScene);
vcAddAndSelectObject(oiB); oiWindow;

%%
sz = oiGet(oiB,'size');
rect = [sz(2)/2,sz(1)/2,5,5];
roiLocs = ieRoi2Locs(rect);

ibIrradiance = vcGetROIData(oiB,roiLocs,'energy');

ibIrradiance = mean(ibIrradiance,1);
vcNewGraphWin; plot(wave,ibIrradiance); grid on; title('ISET irradiance')

%% PTB computes the ptbIrradiance w/m2
d       = displayCreate(monitorName);
wave    = displayGet(d,'wave');

backRGB = [0.5, 0.5, 0.5]';                    % Define background for experment in monitor RGB
backSpd = displayGet(d,'spd')*backRGB;

optics = oiGet(oiB,'optics');
focalLengthMm   = opticsGet(optics,'focal length','mm');
pupilDiameterMm = opticsGet(optics,'pupil diameter','mm');
integrationTimeSec = 0.05;  % Irrelevant for irradiance

[~,~,ptbPhotoreceptors,ptbIrradiance] = ...
    ptbConeIsomerizationsFromSpectra(backSpd,wave,...
    pupilDiameterMm,focalLengthMm,integrationTimeSec,0);

vcNewGraphWin; plot(wave,ptbIrradiance); grid on; title('PTB irradiance')

%% Compare
% The difference is because of the magnification difference between ISET
% and PTB.  See note in doOneSimulation.  Document it here better.
vcNewGraphWin; plot(ibIrradiance,ptbIrradiance,'.')
identityLine;
grid on

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