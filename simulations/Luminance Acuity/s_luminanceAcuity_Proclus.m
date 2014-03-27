%% s_luminanceAcuity_Proclus
%    This script help run coLuminaneAcuity on stanford proclus cluster
%
%
%  (HJ) March, 2014

intensity = [10 15 20 30 40 50 60 70 80 90 100];
totIntensity = cell(length(intensity), 1);
for ii = 1 : length(intensity)
    totIntensity{ii} = intensity(ii)*[1.002 1.005 1.008 1.01 1.02 1.04 1.08 1.1 1.2];
end


cmd = 'params.tIntensity = totIntensity{jobindex};';
cmd = [cmd '[jndIntensity, acc, err, tIntensity] =' ...
            'coLuminanceAcuity(intensity(jobindex), params);'];
cmd = [cmd 'save(sprintf(''~/luminance%d.mat'', jobindex));'];
sgerun2(cmd,'coLuminanceAcuity',1, 1:length(frequency));