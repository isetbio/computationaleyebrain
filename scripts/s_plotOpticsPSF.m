%% s_plotOpticsPSF
%
%  Plot human optics PSF as a function of pupil diameter, defocus and
%  astigmatism
%
%  HJ, VISTA TEAM, 2016

%% Init
ieInit;

%% Plot PSF as a function of pupil diameter
pupil_list = [3 4.5 6 7.5];  % pupil diameter in mm
wave = 400:10:700;           % wavelength samples in nm

for pupilMM = pupil_list
    zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
    wvf = wvfCreate('wave', wave, 'zcoeffs', zCoefs);
    wvf = wvfSet(wvf,'calc pupil size',pupilMM);
    wvf = wvfComputePSF(wvf);
    
    % convert wavefront object to oi structure
    oi = wvf2oi(wvf);
    
    % plot
    oiPlot(oi, 'psf', [], 500);
end

%% Plot PSF as a function of defocus and astigmatism
pupilMM = 3;
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
astig = 0.5;

for defocus = -0.5:0.5:0.5
    % adjust coefficients
    % a simpler way to set defocus is wvfSet(wvf, 'defocusdiopters', val)
    Z = zCoefs;
    Z(5) = Z(5) + wvfDefocusDioptersToMicrons(defocus, pupilMM);
    Z(6) = Z(6) + wvfDefocusDioptersToMicrons(astig, pupilMM);
    
    wvf = wvfCreate('wave', wave, 'zcoeffs', Z);
    wvf = wvfSet(wvf, 'measured pupil size', pupilMM);
    wvf = wvfSet(wvf, 'calc pupil size',pupilMM);
    wvf = wvfComputePSF(wvf);
    
    oi = wvf2oi(wvf);
    oiPlot(oi, 'psf 550');
end