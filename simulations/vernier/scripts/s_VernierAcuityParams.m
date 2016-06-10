%% s_VernierAcuityParams.m
%
%    Test the effects of parameters on Vernier discriminability
%
% HJ, ISETBIO TEAm, 2016

%% Number of samples
%  Init parameters
d = displayCreate('LCD-Apple');
d = displaySet(d, 'dpi', 400);

sampleList = [3000 4200 5500 7500 10000];
acc = zeros(length(sampleList), 1);

% compute discrimination accuracy for using different number of samples
cprintf('*blue', 'Effects of number of samples\n');
for ii = 1 : length(sampleList)
    acc(ii) = coVernier('display', d, 'nFrames', sampleList(ii));
    fprintf('nSamples: %d\t Accuracy:%.2f\n', sampleList(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(sampleList, acc, 'LineWidth', 2);

% save results
res.nSamples.sampleList = sampleList;
res.nSamples.acc = acc;

%% Spatial integration size
d = displayCreate('LCD-Apple');
d = displaySet(d, 'dpi', 400);

intSz = [0.13 0.15 0.17 0.19 0.21]; % sensor field of view in degree
acc = zeros(length(intSz), 1);

cprintf('*blue', 'Effects of spatial integration size\n');
for ii = 1 : length(intSz)
    acc(ii) = coVernier('display', d, 'spatialInt', intSz(ii));
    fprintf('sensor fov: %.2f\t Accuracy:%.2f\n', intSz(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(intSz, acc, 'LineWidth', 2);

% save results
res.sensorFov.intSz = intSz;
res.sensorFov.acc = acc;

%% Exposure time
d = displayCreate('LCD-Apple');
d = displaySet(d, 'dpi', 400);

expTime = [0.03 0.035 0.04 0.045 0.055]; % cone integration time
acc = zeros(length(expTime), 1);

cprintf('*blue', 'Effects of spatial integration size\n');
for ii = 1 : length(expTime)
    acc(ii) = coVernier('display', d, 'expTime', expTime(ii));
    fprintf('exposure time: %.3f\t Accuracy:%.2f\n', expTime(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(expTime, acc, 'LineWidth', 2);

% save results
res.expTime.expTime = expTime;
res.expTime.acc = acc;