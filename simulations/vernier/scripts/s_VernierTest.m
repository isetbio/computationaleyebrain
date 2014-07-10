%% Test and Plot for different ppi of display
%   test svm performance for Vernier acuity under different ppi and patch
%   size
ppiRange = 300:100:800;
nSamples = 3000;
fovRange = [0.02 : 0.02 : 0.1; 0.02:0.02:0.1]'; % 3 minutes to 12 minutes
accPPI   = zeros(length(fovRange), length(ppiRange));
errPPI   = accPPI;

for iFov = 1 : length(fovRange)
    for iPPI = 1 : length(ppiRange)
        fprintf('Exp for FOV: %.2f\t PPI:%d\n', ...
            fovRange(iFov,1), ppiRange(iPPI));
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
for iFov = 1 : length(fovRange)
    plot(ppiRange, accPPI(iFov, :));
end