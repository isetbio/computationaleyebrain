%% s_coneContrastTest
%
%

totContrast{1} = [.1 .05 .04 .03 .02 .01 .005 0.004 0.003 .002];
totContrast{2} = [.4 .2 .1 .05 .04 .03 .02 .01 .005 0.004];
totContrast{3} = [.8 .6 .4 .2 .1 .05 .04 .03 .02 .01 .005];
totContrast{4} = [1 .8 .6 .4 .2 .1 .05 .02 .01];
totContrast{5} = [1 .8 .6 .4 .2 .1 .05];
totContrast{6} = [1 .8 .6 .4 .2];
totContrast{7} = [1 .8 .6 .4];
totContrast{8} = [1 .8 .6 .4];
resultA = cell(8, 1);
resultE = cell(8, 1);
fqIndx = 1;
for frequency = [5 10 12 15 17 20 30 40]
    fprintf('frequency:%d\n', frequency);
    for colorIndx = 1 : length(totContrast{fqIndx})
        contrast = totContrast{fqIndx}(colorIndx);
        fprintf('testColor: %f\n', contrast);
        s_contrastSensitivity;
        resultA{fqIndx}(colorIndx) = acc;
        resultE{fqIndx}(colorIndx) = err;
        fprintf('freq:%d\t contrast:%f\t acc:%f\n', frequency, contrast, acc);
    end
    fqIndx = fqIndx + 1;
end

save resultALM2.mat resultA resultE