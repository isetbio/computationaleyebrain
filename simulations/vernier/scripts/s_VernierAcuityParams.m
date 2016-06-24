%% s_VernierAcuityParams.m
%
%    Test the effects of parameters on Vernier discriminability
%
% HJ, ISETBIO TEAm, 2016

%% Viewing distance
vDistList = [0.7:0.1:1.2 1.5 2];
acc = zeros(length(sampleList), 1);

% compute discrimination accuracy for different viewing distance
cprintf('*blue', 'Effects of viewing distance\n');
for ii = 1 : length(vDistList)
    acc(ii) = coVernier('vDist', vDistList(ii));
    fprintf('vDist: %.1f\t Accuracy:%.2f\n', vDistList(ii), acc(ii));
end

% plot
vcNewGraphWin; plot(vDistList, acc, 'LineWidth', 2);
xlabel('Viewing distance (m)'); ylabel('Discrimination Probability');
grid on;

% save results
res.vDist.vDistList = vDistList;
res.vDist.acc = acc;
    
%% Number of samples
%  Init parameters
sampleList = [500 1000 2000 3000 4200 5500 7500 10000];
acc = zeros(length(sampleList), 1);

% compute discrimination accuracy for using different number of samples
cprintf('*blue', 'Effects of number of samples\n');
for ii = 1 : length(sampleList)
    acc(ii) = coVernier('nFrames', sampleList(ii));
    fprintf('nSamples: %d\t Accuracy:%.2f\n', sampleList(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(sampleList, acc, 'LineWidth', 2);
xlabel('Number of samples'); ylabel('Discrimination Probability');
grid on;

% save results
res.nSamples.sampleList = sampleList;
res.nSamples.acc = acc;

%% Spatial integration size
intSz = [0.13 0.15 0.17 0.19 0.21]; % sensor field of view in degree
acc = zeros(length(intSz), 1);

cprintf('*blue', 'Effects of spatial integration size\n');
for ii = 1 : length(intSz)
    acc(ii) = coVernier('spatialInt', intSz(ii));
    fprintf('sensor fov: %.2f\t Accuracy:%.2f\n', intSz(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(intSz, acc, 'LineWidth', 2);
xlabel('Cone mosaic size (deg)'); ylabel('Discrimination Probability');
grid on;

% save results
res.sensorFov.intSz = intSz;
res.sensorFov.acc = acc;

%% Exposure time
expTime = [0.03 0.035 0.04 0.045 0.055]; % cone integration time
acc = zeros(length(expTime), 1);

cprintf('*blue', 'Effects of exposure time\n');
for ii = 1 : length(expTime)
    acc(ii) = coVernier('expTime', expTime(ii));
    fprintf('exposure time: %.3f\t Accuracy:%.2f\n', expTime(ii), acc(ii));
end

% plot the curve
vcNewGraphWin; plot(expTime, acc, 'LineWidth', 2);
xlabel('Cone integration time (s)'); ylabel('Discrimination Probability');
grid on;

% save results
res.expTime.expTime = expTime;
res.expTime.acc = acc;

%% Bar vs background contrast
barColor = [0.55 0.6 0.7 0.8 0.9];
acc = zeros(length(barColor), 1);

cprintf('*blue', 'Effects of bar contrast\n');
for ii = 1 : length(barColor)
    acc(ii) = coVernier('barColor', barColor(ii));
    fprintf('barColor: %.2f\t Accuracy:%.2f\n', barColor(ii), acc(ii));
end

% plot
vcNewGraphWin; plot(barColor, acc, 'LineWidth', 2);
xlabel('Bar color'); ylabel('Discrimination Probability');
grid on;

% save results
res.contrast.barColor = barColor;
res.contrast.acc = acc;

%% Cone isolated color bar
d = displayCreate('LCD-Apple', 'dpi', 400, 'gamma', 'linear');

% get cone isolated color direction
bgLMS = [0.5 0.5 0.5] * displayGet(d, 'rgb2lms');
acc = zeros(3, 5);
dirName = ['L', 'M', 'S'];
vDistList = 0.5:0.1:1;

% create human sensor
sensor = sensorCreate('human');

for ii = 1 : length(dirName)
    stimuliLMS = bgLMS; stimuliLMS(ii) = stimuliLMS(ii) * 1.15;
    barColor = stimuliLMS * displayGet(d, 'lms2rgb');
    for jj = 1 : length(vDistList)
        acc(ii, jj) = coVernier('barColor', barColor, ...
            'nFrames', 3000, 'emFlag', [1 1 1], 'vDist', vDistList(jj), ...
            'sensor', sensor, 'display', d);
        fprintf('Color direction: %s\t vDist:%.2f\t Accuracy:%.2f\n', ...
            dirName(ii), vDistList(jj), acc(ii, jj));
        % vcNewGraphWin; imagesc(w);
    end
end

% save results
res.color.dirName = dirName;
res.color.vDistList = vDistList;
res.color.acc = acc;