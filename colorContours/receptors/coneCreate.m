function cone = coneCreate(species, varargin)
% Create a cone structure that includes pre-retinal filter parameters
%
% The pre-retinal filters and related cone parameters are
%
%   macular pigment structure
%   lens structure
%   cone photopigment optical densities
%
% Absorbance spectra:   Normalized to a peak value of 1.
% Optical density (OD): Pigment specific
% Absorbtance spectra:  Proportion of quanta actually absorbed
%
%    absorbtanceSpectra = 1 - 10.^(-OD * absorbanceSpectra), where
%
% Example:
%    cone = coneCreate('human');
%    cone = coneCreate('human','sensor',sensorCreate);
%
% See also:
%   coneGet, coneSet, macularCreate, lensCreate
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check 
if notDefined('species'),     species = 'human'; end

species = ieParamFormat(species);

%% Create cone structure
switch species
    case {'human','humanfovea'}
        
        cone.species = 'human';
        cone.type    = 'cone';
        cone.wave    = (400:10:700)';
        
        % Default macular pigment has density of 0.28.
        cone.macular = macularCreate;
        
        % lens transmittance, default stockman's lens pigment densitiy
        cone.lens  = lensCreate;
        % lensTransmittance(...
          %  modelParams.lensTransF, modelParams.wave, ...
           % modelParams.lensTransSrc);
        
        % cone photopigment optical densities for L,M,S
        cone.opticalDensity = [0.5 0.5 0.4]';
        
        % peak spectrum position
        % Only stores for L,M,S and peak for K is ignored
        % cone.peakLambda = modelParams.peakLambda;
        
        % Normalized spectral absorption (absorbance) of the three cone
        % types.  We NEED TO GET THIS RIGHT, not load the stockmanQuanta.
        % We should load the Sharpe absorbance data like PTB.
        absorptance = ieReadSpectra('stockmanQuanta',cone.wave);
        % cone.absorbance = log10(1 - absorptance)/-cone.opticalDensity;
        cone.absorbance = log10(1 - absorptance)*diag( 1 ./ -cone.opticalDensity);


        % quantal efficiency
        cone.quantalEfficiency = [.3 .3 .3];
        
    otherwise
        error('Unknown type encountered');
end


%% Now over-write default parameters with what user sent in

% Check that we have an even number of arguments
n = length(varargin);
if isodd(n), error('Must have param,val pairs'); end

% Set the parameter,val pairs into the cone structure
for ii=1:2:(n-1)
    cone = coneSet(cone,varargin{ii},varargin{ii+1});
end

return

