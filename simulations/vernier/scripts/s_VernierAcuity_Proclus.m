%% s_VernierAcuity_Proclus
%    This script helps run vernier acuity in parallel on
%    Proclus
%
%  (HJ) March, 2014

barWidth = [15 18 20 22 25 30 32 35 38 40 45];
cmd = 'params.barWidth = barWidth(jobindex);';
cmd = [cmd '[jndDist, acc, err, tDist] =' ...
            'coVernierAcuity(params);'];
cmd = [cmd 'save(sprintf(''~/vernier%d.mat'', jobindex));'];
sgerun2(cmd,'vernierAcuity',1, 1:length(barWidth));
