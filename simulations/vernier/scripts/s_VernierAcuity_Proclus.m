%% s_VernierAcuity_Proclus
%    This script helps run vernier acuity in parallel on
%    Proclus
%
%  (HJ) March, 2014

barWidth = [1 2 3 5 8 10 15 20];
cmd = 'params.barWidth = barWidth(jobindex);';
cmd = [cmd '[jndDist, acc, err, tDist] =' ...
            'coVernierAcuity(params);'];
cmd = [cmd 'save(sprintf(''~/vernier%d.mat'', jobindex));'];
sgerun2(cmd,'vernierAcuity',1, 1:length(barWidth));
