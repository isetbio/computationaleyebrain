% wavelength samples
wave = 400:700;

% Compute equal energy light position
energy = ones(length(wave), 1) * 0.001;
XYZ_ee = ieXYZFromEnergy(energy', wave);
LMS_ee = xyz2lms(XYZ_ee);

% Compute position for invariant wavelength
% Here, we do for deutan and protan: 475 and 575
invWave = [475 575];
LMS_inv = zeros(length(invWave), 3);
for ii = 1 : length(invWave);
    energy = zeros(length(wave), 1);
    energy(wave == invWave(ii)) = 0.15;
    XYZ_inv = ieXYZFromEnergy(energy', wave);
    LMS_inv(ii,:) = xyz2lms(XYZ_inv);
end

% Plot points and half-planes
vcNewGraphWin; hold on;
lc = [0.9 0.6 0.2];
plot3([0 LMS_ee(1) LMS_inv(:, 1)'], ...
      [0 LMS_ee(2) LMS_inv(:, 2)'], ...
      [0 LMS_ee(3) LMS_inv(:, 3)'], ...
      'o', 'MarkerSize', 10, 'LineWidth', 2);
plot3([0 LMS_ee(1) LMS_inv(1, 1), 0], ...
      [0 LMS_ee(2) LMS_inv(1, 2), 0], ...
      [0 LMS_ee(3) LMS_inv(1, 3), 0], ...
      '-', 'LineWidth', 2, 'Color', lc);
plot3([0 LMS_ee(1) LMS_inv(2, 1), 0], ...
      [0 LMS_ee(2) LMS_inv(2, 2), 0], ...
      [0 LMS_ee(3) LMS_inv(2, 3), 0], ...
      '-', 'LineWidth', 2, 'Color', lc);
view(-15, 25); xlabel('L'); ylabel('M'); zlabel('S'); grid on;

% Plot confusion line for protanope
M = 40; S = 20;
plot3([0 120], [M M], [S S], '--', 'LineWidth', 2);
LMS = lms2lmsDichromat(reshape([0 M S], [1 1 3]), 'protan', [], LMS_ee);
LMS = LMS(:);
plot3(LMS(1),M,S, 'x', 'MarkerSize', 10, 'LineWidth', 2);