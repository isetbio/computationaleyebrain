%% s_wavelengthAcuity
%   Compute sensitivity for different color
%   The computation is done in cone absorption level, openency and rgc
%   features are not included
%
%  (HJ) March, 2014

refWave = 480:10:580;
jndWave = zeros(length(refWave), 1);

for ii = 1 : length(refWave)
    fprintf('refWave:%d\t\t', refWave(ii));
    jndWave(ii) = wavelengthAcuity(refWave(ii));
    fprintf('jndWave:%.1f\t dist:%.1f\n', jndWave(ii), jndWave(ii)-refWave(ii));
end

% Plot
