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
% Absorptance spectra:  Proportion of quanta actually absorbed
%
%    absorptanceSpectra = 1 - 10.^(-OD * absorbanceSpectra), where
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
        cone.name    = 'default human cone';
        cone.wave    = (400:10:700)';
        
        % Default macular pigment has density of 0.28.
        cone.macular = macularCreate([],cone.wave);
        
        % lens transmittance, default stockman's lens pigment densitiy
        cone.lens  = lensCreate([],cone.wave);
        % lensTransmittance(...
          %  modelParams.lensTransF, modelParams.wave, ...
           % modelParams.lensTransSrc);
        
        % photopigment optical densities for L,M,S
        cone.opticalDensity = [0.5 0.5 0.4]';
        
        % peak spectrum position
        % Only stores for L,M,S and peak for K is ignored
        % cone.peakLambda = modelParams.peakLambda;
        
        
        % We should load the Sharpe absorbance data like PTB.
        %    absorptance = 10.^ieReadSpectra('coneAbsorbance',cone.wave);
        %
        % Note that the original data is from Stockman Sharp and it's in
        % log10 mode. Thus, we need to convert it back.
        %
        % The relationship between absorbance and absorptance is as below:
        % cone.absorbance = log10(1 - absorptance)*...
        %                               diag( 1 ./ -cone.opticalDensity);

        cone.absorbance = 10.^ieReadSpectra('coneAbsorbance',cone.wave);
        
        % Peak absorptance efficiency
        % Note that the actual peak efficiency is the product of this value
        % and the opticalDensity. The effective value should be around .3
        % for all three types of cones
        cone.peakEfficiency = [2 2 2]/3;

        cone.spatialDensity = [0 .6 .3 .1];  % No blanks, L,M,S
        
    otherwise
        error('Unknown type encountered');
end


%% Over-write default parameters with what user sent in
% Check that we have an even number of arguments
n = length(varargin);
if isodd(n), error('Must have param,val pairs'); end

% Set the parameter,val pairs into the cone structure
for ii=1:2:(n-1)
    cone = coneSet(cone,varargin{ii},varargin{ii+1});
end

return

