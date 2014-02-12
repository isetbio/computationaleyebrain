% decompStockLens.m

% Load this from the PsychToolBox
load('den_lens_ssf');

wavelength  = 390:830;
lenspigment = den_lens_ssf;

% cropped the curve from 460 nm to 730 and extrapolated linearly to get
% dens_lens1, which can vary with aging or other factors. 
cropWave = 460:730;
ind460nm = find(wavelength == min(cropWave));
ind730nm = find(wavelength == max(cropWave));
dens_lens1_partial = lenspigment(ind460nm:ind730nm);
dens_lens1 = interp1(cropWave, dens_lens1_partial, wavelength,'linear','extrap')';

% subtructed lens1 function from Stockman and interpolated linearly again
% to get lens2 function, which is fixed.

dens_lens2 = lenspigment - dens_lens1;

% cropped the curve from 390 nm to 730 nm because all of calculations have
% been done between the range.

comments = 'Lens density as a function of wavelength from stockman 1999. It was decomposed into two factors based on Pokorny paper (1987).';
dens_lens1 = dens_lens1(1:ind730nm);
dens_lens2 = dens_lens2(1:ind730nm);
S_lens_s2f = [390 1 341];
wave       = [390:730]';

%
fName = fullfile(cmPublicRootPath,'calibration','decompStockLens.mat');
save(fName, 'comments', 'dens_lens1', 'dens_lens2', 'S_lens_s2f', 'wave')
