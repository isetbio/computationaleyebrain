function [sensorIsolateAmp1P AbsSpectra densities sensors] = cm_prepCalibrationMtx(fovealflag, PhotmeternameForAxis)
% [sensorIsolateAmp1P AbsSpectra densities] = ...
%           cm_prepCalibrationMtx(fovealflag, PhotmeternameForAxis)
%
% This function prepares pigment-isolated stimulus, LMS cones and
% estimated melanopsin absorbance and photopigment optical densities that
% were used for this calculation.
%
% <Input>
%   fovealflag           ... foveal data or not
%   PhotmeternameForAxis ... 'pr650' or 'pr715' (char)
%
% <Output>
%   sensorIsolateAmp1P   ... pigment-isolated stimulus (delta = 1%)
%   AbsSpectra           ... LMS cones and estimated melanopsin absorbance
%   densities            ... photopigment optical densities that were used for this calculation
%   sensors              ... sensor response functions
%
% see also cm_CalibrationRevise.m
%
% HH (c) Vista lab 2012. 
% 
%% prep
wls = cm_getDefaultWls;

if ~exist('fovealflag','var') || isempty(fovealflag)
    fovealflag = false;
    disp('calc peripheral isolated stimuli.')
end

if ~exist('PhotmeternameForAxis','var') || isempty(PhotmeternameForAxis)
    PhotmeternameForAxis = 'pr715';
    disp('SPDs were calibrated with PR-715.')
end

display.wavelengths = wls;
display.spectra     = cm_getledSPD(display.wavelengths, PhotmeternameForAxis);
[LMSperi, LMSfovea] = ct_loadStandardObserverData;

if fovealflag == false;
    LMS = LMSperi;
    modelparams = [1 0];
else
    LMS = LMSfovea;
    modelparams = [1 0.28];
end


%% LMS absorbance spectra
[absorbanceSpectra ,~, densities] = cm_loadLMSabsorbance(fovealflag);

%% estimate melanopsin absorbance and response function

% melanopsin parmeters based on Dacey 2005
melanopsinpeak = 482; axialphotopigmentdensity = 0.5; 

% melanopsin spectral absorbance 
melanopsinAbsorbance = PhotopigmentNomogram(wls',melanopsinpeak,'StockmanSharpe');
AbsSpectra = [absorbanceSpectra melanopsinAbsorbance'];

% melanopsin response function
melanopsinResponse = ...
    cm_SensorRespfromNomogram(melanopsinpeak, modelparams, axialphotopigmentdensity, wls);

% add melanopsin optical density in this calc
densities(end+1) = axialphotopigmentdensity;

%% calc pigment-isolated stimulus
% 3 cones and melanopsin absobance
sensors = [LMS melanopsinResponse];

% get stimulus we used as isolated-stimuli
backRYGCBM.dir = ones(size(display.spectra,2),1);
backRYGCBM.scale = .5;
numSens = size(sensors, 2);

% calc LED amplitude of pigment-isolated stimuli
for ii = 1:numSens
    
    stimLMSRm.dir     = zeros(numSens,1);
    stimLMSRm.dir(ii) = 1;
    
    [lmsrmContrast, ~, stimRYGCBM] = ct_MakeIsoStims(display, sensors, stimLMSRm, backRYGCBM);
    MaxContrasts(ii)     = lmsrmContrast(ii);
    sensorIsolateAmp(:,ii)   = stimRYGCBM.dir;
    sensorIsolateAmp1P(:,ii) = stimRYGCBM.dir / (lmsrmContrast(ii) * 100);
    
end


end
