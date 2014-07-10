function [lmsrmContrast lmsrmBack stimRYGCBM backxy] = ct_MakeIsoStims(display, sensors, stimLMSRm, backRYGCBM, wls, backgroundSPD)
% [lmsrmContrast lmsrmBack stimRYGCBM backxy] = ct_MakeIsoStims(display, sensors, stimLMSRm, backRYGCBM, wls, scalingfactor)
% 
% <Input>
%   display ... ISET display structure containing 341 x n display primaries
%   sensors ... LMS fundamentals + alpha (rods and/or melanopsin) is 341 x n
%   stimLMSRm . stimLMSRm.dir (a stimulus direction) is n x 1 
%   backRYGCBM. .dir is m x 1 (all of them are 1)
%               .scale is 0.5
%   wls     ... wavelength is 1 x 341 (390:730)
%   scalingfactor is 1.
%
%
% <Output>
%   lmsrmContrast ... a contrast of sensor is n x 1
%   lmsrmBack     ... a sensor response to Background is n x 1
%   stimRYGCBM    ...dir stimulus direction of primaries is m x 1
%                   .maxScale 
%                   .scale
%   backxy        ... chormaticity of background is 1 x 2
%
% see also script_fittingprocedureSGE.m
%
% C) Vista Lab, HH 2012 

%% set params
if ieNotDefined('wls'),             wls           = cm_getDefaultWls;  end
if ieNotDefined('backgroundSPD'),   backgroundSPD = zeros(size(wls))';  end
% if ieNotDefined('scalingfactor'),   scalingfactor = 1;       end % scalingfactor = [1 0.5 0.25 0.1];

% transpose for calculation
% scalingfactor = scalingfactor(:)';
this.spectra = displayGet(display,'spd');

% spectrum and chromaticity of backgound 
backSpectrum = this.spectra * backRYGCBM.dir + backgroundSPD;
backXYZ      = ieXYZFromEnergy(backSpectrum', wls);
backxy       = chromaticity(backXYZ);
scalingfactor= 1; 
%%  Calculate the maximum scale factor for the stimulus 

[stimLMSRm stimRYGCBM] = ...
    findMaxConeRodMelanopsinScale(this, stimLMSRm, backRYGCBM, sensors, backgroundSPD);

stimLMSRm.scale = stimLMSRm.maxScale * scalingfactor;

%% Calculate the values of the primaries

stimRYGCBM = conesrm2RYGCBM(this, stimLMSRm, backRYGCBM, sensors, backgroundSPD);

%% Compute the vector of cone contrasts

[lmsrmContrast lmsrmBack] = ...
    RYGCBM2conesrmContrast(this, stimRYGCBM, backRYGCBM, sensors, backgroundSPD);

return
