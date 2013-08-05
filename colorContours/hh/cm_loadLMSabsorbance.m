function [absorbanceSpectra wls opticaldensities]= cm_loadLMSabsorbance(fovealflag, wls)
%  [absorbanceSpectra wls opticaldensities]= ...
%            cm_loadLMSabsorbance(stLMSname, receptorTypes, wls)
% 
% This function loads the absorbance spectra. They were aquired by deviding
% stockman LMS fucntios by standard lens and macular pigment transmittance.
%
%
% absorbanceSpectra = cm_loadLMSabsorbance(fovealflag, receptorTypes);
%
% see also cm_prepCalibrationMtx.m
%
% HH (c) Vista lab 2012. 
% 
%%
if ~exist('wls','var') || isempty(wls)
    wls = cm_getDefaultWls;
end

% for PhotopigmentAxialDensity.m (from PTB)
if fovealflag == false;
    stLMSname   = 'stockmanLMSwithoutMacLensPeri';
    receptorTypes = {'LCone','MCone','SCone'};
    
elseif fovealflag == true;
    stLMSname   = 'stockmanLMSwithoutMacLens';
    receptorTypes = {'FovealLCone','FovealMCone','FovealSCone'};
end

stLMSabsorbtance = load(stLMSname);
species = 'Human'; source = 'StockmanSharpe'; 

% get absorbanceSpectra (which required lens transmittance and axial pigment density to become LMS function)
for ii = 1:length(receptorTypes)
    
    % photopigment optical densities (from psychtoolbox)
    opticaldensities(ii) = PhotopigmentAxialDensity(receptorTypes{ii},species,source);
    
    % calc absorbance from absorbtance
    [absorbanceSpectra(:,ii), absorbanceSpectraWls] = AbsorbtanceToAbsorbance(stLMSabsorbtance.data(:,ii), wls, opticaldensities(ii));

    % for comfirmation, recombert to absorbtance
    % abtcheck(:,ii) = AbsorbanceToAbsorbtance(absorbanceSpectra(:,ii)', wls, opticaldensities(ii));
end
