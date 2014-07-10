function [absorbanceSpectra, absorbanceSpectraWls] =...
    AbsorbtanceToAbsorbance(absorbtanceSpectra, absorbtanceSpectraWls, axialOpticalDensities)
%
% [absorbanceSpectra, absorbanceSpectraWls] =...
%   AbsorbtanceToAbsorbance(absorbtanceSpectra, absorbtanceSpectraWls, axialOpticalDensities);
%
% This code is just a flipside of AbsorbanceToAbsorbtance in Psychotoolbox.
% These words were defined as below here;
% Absorbance indicates photopigment's absorption.
% Absorbtacne indicates cone's absorption.
%
% Exapmple
%
% % try to estimate L-cone absorbance from the data of stockman 2000.
%
% % LcDenseF is axial pigment density of foveal L-cone
% LcDenseF = PhotopigmentAxialDensity('FovealLCone','Human','StockmanSharpe');
%
% % load stockman's data
% st2 = load('stockman2000');   Lcent = st2.data(:,1); % L-cone should be first column.
% t = macular(0.28,st2.wavelength); % use 'standard' macular denstity
%
% % lens and macular pigment transmittance at fovea
% transF = t.transmittance' .* LensTransmittance(st2.wavelength);
% % LAtF is L-cone Absorbtance at fovea 
% LAtF = Lcent' ./ transF;
%
% % convert Absorbtance to Absorbance 
% LAbF = AbsorbtanceToAbsorbance(LAtF, st2.wavelength, LcDenseF);
% figure, plot(st2.wavelength,LAbF,'r','LineWidth',3)
% hold on, plot(st2.wavelength, st2.data(:,1),'rx');
% legend({'Absorbance','Absorbtacne'})
%
%
% Copyright HH, Vista Lab, 2010
%
%%

if ~exist('absorbtanceSpectra','var'); Help AbsorbtanceToAbsorbance; return; end
if ~exist('axialOpticalDensities','var'); disp('axialOpticalDensities is required.'); return; end
if ~exist('absorbtanceSpectraWls','var'); absorbtanceSpectraWls = []; end

%%
% Normalize absorbtanceSpectance, should have (1-10^-axialOpticalDensities) as maximam value
absorbtanceSpectra = absorbtanceSpectra ./ max(absorbtanceSpectra) * (1-10^-axialOpticalDensities);
absorbanceSpectra = (- 1 ./ axialOpticalDensities) .* log10(absorbtanceSpectra - 1);

absorbanceSpectra = real(absorbanceSpectra); % remove small value of imaginary number
absorbanceSpectraWls = absorbtanceSpectraWls; % just mimic the original file

return
