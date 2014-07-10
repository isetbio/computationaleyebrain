%% s_waveDiscrimination_Proclus
%    This script helps run wavelength discrimination in parallel on
%    Proclus
%
%  (HJ) March, 2014

refWave = 400:10:700;
refWave = [refWave' refWave'];

cmd = [cmd '[jndWave, acc, err, tWave] =' ...
            'wavelengthDiscrimination(refWave(jobindex,:));'];
cmd = [cmd 'save(sprintf(''~/waveDisc%d.mat'', jobindex));'];
sgerun2(cmd,'waveDiscrimination',1, 1:length(refWave));
