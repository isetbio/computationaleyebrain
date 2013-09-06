function [LMSQuanta, LMSEnergyFunc] = humanConeQE(inertP, melanopsinflag)
% [cFilters, LMSEnergyFunc] = ieReadHumanQE(inertP)
% 
% Input
%    inertP ... is created by odParams
%     (fields).lens
%             .macular
%             .LPOD
%             .MPOD
%             .SPOD
%             .melPOD
%             .wave
%
% Output
%   cFilters ... LMS cones and melanposin function with quantal unit
%
%   inertFilter .. if it's IR filter, should be 1 at all wavelength
%
%   LMSEnergyFunc ... LMS energy function
%
%
% See also: odParams.m
%
% (c) VISTA lab 2012 HH

%% The PTB way is right here.  We mimic it.

%% TODO
%  Use the structure as a photoreceptor structure for photoreceptorCreate
%  Separate out the sets/gets using the functions inside of
%    FillInPhotoreceptors for some of the gets
%  
%  whatCalc = 'LivingHumanFovea';
%  photoreceptors = DefaultPhotoreceptors(whatCalc);
%  photoreceptors.eyeLengthMM.source = '17';
%  photoreceptors.nomogram.S = [380,4,101];
%  photoreceptors = FillInPhotoreceptors(photoreceptors);
%
%% 

if ~exist('melanopsinflag','var') || isempty(melanopsinflag)
    melanopsinflag = false;
end

% get LMS cone and melanopsin resuponse function at energy unit
modelparams = [inertP.lens inertP.macular inertP.LPOD inertP.MPOD inertP.SPOD inertP.melPOD];

%
LMSEnergyFunc = cm_PigmentResponsefunction(inertP.visfield, modelparams, melanopsinflag);

% change an unit from energy to quanta for ISET calculation
q2e         = Quanta2Energy(inertP.wave, ones(length(inertP.wave), 1));
LMSQuanta   = diag(q2e) * LMSEnergyFunc;

% normalize
mx          = max(LMSQuanta);
LMSQuanta    = LMSQuanta*diag(1./mx);

end

function PigResfunc = cm_PigmentResposefunction(foveaflag, modelparams, melanopsinflag)
% PigResfunc = cm_PigmentResposefunction(foveaflag, modelparams, melanopsinflag)
%
% <Input>
%   foveaflag       ... fovea or not
%   modelparams     ... model parameters: 1) lens pigment density parameter
%                       2) macular pigment density 3)-6) photopigment
%                       optical density 
%
%   melanopsinflag  ... add melanopsin response function or not
%
% <Output>
%
%   PigResfunc      ... response functions for each photopigment
%
%
%
% (c) VISTA lab 2012 HH
%%
switch foveaflag
    case {1,'f','fovea','fov'}
        foveaflag = true;
    case {0,'p','peri','periphery'}
        foveaflag = false;
end

wls = cm_getDefaultWls;

lensfactor = modelparams(1);
macfactor  = modelparams(2);

% get LMS absorbance
absorbanceSpectra = cm_loadLMSabsorbance(foveaflag);

% need melanopsin?
if melanopsinflag
    
    melpeak           = 482;
    melabsorbance     = PhotopigmentNomogram(wls',melpeak,'StockmanSharpe');
    absorbanceSpectra = [absorbanceSpectra melabsorbance'];
    
    PODs = modelparams(3:6);
else
    PODs = modelparams(3:5);
end


%% transmittance of lens and macular
% lens transmittance function
lT = cm_LensTransmittance(lensfactor, wls,'stockman2');

% macular transmittance function
mT = macular(macfactor, wls);

% whole eye transmittance function
eT = lT .* mT.transmittance';

% assume no peak shift
Lambdashift = 0;

%%  pigment response function
PigResfunc = cm_variableLMSI_PODandLambda(absorbanceSpectra, PODs, Lambdashift, eT);

end
