%% s_contrastSensitivity_Proclus
%    This scripts help initialize and allocate jobs on proclus for
%    computing CSF by computation observer
%
%
%  (HJ) March, 2014

frequency = [5 10 12 15 17 20 30 40];
totContrast{1} = [.1 .05 .04 .03 .02 .01 .005 0.004 0.003 .002];
totContrast{2} = [.4 .2 .1 .05 .04 .03 .02 .01 .005 0.004];
totContrast{3} = [.8 .6 .4 .2 .1 .05 .04 .03 .02 .01 .005];
totContrast{4} = [1 .8 .6 .4 .2 .1 .05 .02 .01];
totContrast{5} = [1 .8 .6 .4 .2 .1 .05];
totContrast{6} = [1 .8 .6 .4 .2];
totContrast{7} = [1 .8 .6 .4];
totContrast{8} = [1 .8 .6 .4];


cmd = 'params.testContrast = totContrast{jobindex};';
cmd = [cmd '[jndContrast, acc, err, tContrast] =' ...
            'coContrastSensitivity(frequency(jobindex), params);'];
cmd = [cmd 'save(sprintf(''~/CSF%d.mat'', jobindex));'];
sgerun2(cmd,'coCSF',1, 1:length(frequency));
