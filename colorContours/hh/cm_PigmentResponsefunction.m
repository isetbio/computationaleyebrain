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