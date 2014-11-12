%% s_waveDiscrimination_Proclus
%    This script helps run wavelength discrimination in parallel on
%    Proclus
%
%    If proclus is not available, it will run locally instead.
%
%  (HJ) ISETBIO TEAM, 2014

refWave = 400:5:600;

try % try using proclus to accelerate computation
    cmd = ['[jndWave, expData] =' ...
      'coWaveDiscrimination([refWave(jobindex) refWave(jobindex)+0.01]);'];
    cmd = [cmd 'save(sprintf(''~/waveDisc%d.mat'', jobindex));'];
    sgerun2(cmd,'waveDiscrimination',1, 1:length(refWave));
catch
    jndWave = zeros(length(refWave), 1);
    for ii = 1 : length(refWave)
        fprintf('refWave:%d\t\t', refWave(ii));
        [jndWave(ii), expData] =  ...
            coWaveDiscrimination([refWave(ii) refWave(ii)+0.01]);
        fprintf('jndWave:%.1f\t dist:%.1f\n', ...
                 jndWave(ii), jndWave(ii)-refWave(ii));
    end
end