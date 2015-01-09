%% s_bradleyJOV2014
%
%    This script replicates Bradley et al. paper in JOV 2014
%
%       Retina-V1 model of detectability across the visual field
%
%    This scirpt compares the components, performance and results between
%    their implementation and ISETBIO implementation
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit;

%% Compare human optics
%  create human optics in ISETBIO
oi = oiCreate('human'); % standard human optics

%  get otf function of human optics
otf = oiGet(oi, 'optics otf', 540); % otf at 540 nm

%  get frequency support in cycles / degree
fx = oiGet(oi, 'optics otf fx');   % frequency support in cycles / mm
fx = fx * tand(1) * oiGet(oi, 'focal length') * 1000; % units conversion
fx = fx(fx >= 0); % use non-negative frequencies only

%  plot mtf function
vcNewGraphWin; grid on; hold on;
plot(fx, otf(1, 1:size(otf, 2)/2), '--r');

% plot MTF in Bradley et al. paper
% MTF function is given by function:
%   MTF(f) = 0.78 exp(-0.172f) + 0.22 exp(-0.037f)
MTF = @(f) 0.78 * exp(-0.172*f) + 0.22 * exp(-0.037*f);
plot(fx, MTF(fx), '-b');

% add label, title and lengend
xlabel('Spatial Frequency (cycles/deg)'); ylabel('MTF');
title('Comparison of Optics'); legend('ISETBIO Optics', 'Bradley Paper');

% The MTF function used in the paper is not wavelength dependent and it
% seems like it's close to MTF in ISETBIO for wavelength around 540 nm

%% Light adaptation


%% 