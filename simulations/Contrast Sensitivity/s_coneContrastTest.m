%% s_coneContrastTest
%
%

cRange = [0 0.2 0.3 0.4 0.42 0.44 0.45 0.46 0.47 0.48];
resultA = zeros(7, length(cRange));
resultE = zeros(7, length(cRange));
for ppc = 7:-1:1
    fprintf('ppc:%d\n', ppc);
    for colorIndx = 1 : length(cRange)
        testColor = cRange(colorIndx);
        fprintf('testColor: %f\n', testColor);
        for iter = 1 : 2
            s_contrastSensitivity;
            resultA(ppc, colorIndx) = resultA(ppc, colorIndx) + acc/2;
            resultE(ppc, colorIndx) = resultE(ppc, colorIndx) + err/2;
        end
        fprintf('ppc:%d\t color:%f\t acc:%f\n', ppc, testColor, resultA(ppc, colorIndx));
    end
end

save resultALM2.mat resultA resultE