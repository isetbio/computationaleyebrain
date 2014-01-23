%% Test and Plot for different ppi of display
%   test svm performance for Vernier acuity under different ppi and patch
%   size
ppiRange = 300:50:1200;
nSamples = 3000;
fovRange = [0.2 : 0.1 : 1.2; 0.2:0.1:1.2]';
accPPI   = zeros(length(fovRange), length(ppiRange));
errPPI   = accPPI;

for iFov = 1 : length(fovRange)
    for iPPI = 1 : length(ppiRange)
        imgFov = fovRange(iFov, :);
        ppi    = ppiRange(iPPI);
        s_VernierAcuity;
        accPPI(iFov, iPPI) = acc(1);
        errPPI(iFov, iPPI) = acc(2);
    end
end

% save result
save expResults.mat ppiRange nSamples fovRange accPPI errPPI

vcNewGraphWin; hold on;
for iFov = 1 : 6
    plot(ppiRange, accPPI(iFov, :));
end