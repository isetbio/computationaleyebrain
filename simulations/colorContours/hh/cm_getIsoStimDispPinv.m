function [IsoAmp1P MaxCont] = cm_getIsoStimDispPinv(display, sensors)
% Color matching get cone isolating stimulus display with pseudoinverse.
%
%   [IsoAmp1P MaxCont] = cm_getIsoStimDispPinv(display, sensors)
%
% Input:
%   display:  A display structure
%   sensors:  Matrix with columns equal to sensors through inert pigments
%
% Hiroshi (c) Stanford VISTASOFT Team, 2012

%%
backRYGCBM.dir   = ones(displayGet(display,'n primaries'),1);
backRYGCBM.scale = 0.5;

nLCDs     = displayGet(display,'n primaries');
numSens   = size(sensors,2);
MaxCont   = zeros(1,numSens);
IsoAmp1P  = zeros(nLCDs,numSens);

for ii = 1:numSens
    stimLMSRm.dir     = zeros(numSens,1);
    stimLMSRm.dir(ii) = 1;        
    [lmsrmContrast,~,stimRYGCBM] = ct_MakeIsoStims(display, sensors, stimLMSRm, backRYGCBM);
    MaxCont(ii)    = lmsrmContrast(ii);
    IsoAmp1P(:,ii) = (stimRYGCBM.dir / (lmsrmContrast(ii) * 100)); 
end

end