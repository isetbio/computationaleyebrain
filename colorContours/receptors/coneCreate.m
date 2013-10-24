function cone = coneCreate(species, modelParams)
% Create a cone structure that includes pre-retinal filter parameters
%
% See the PTB routines,
%     DefaultPhotoreceptors
%     ptbConeIsomerizationsFromSpectra
%     FillInPhotoreceptors
%   
%  And check HH's routines.  They need to be integrated into isetbio
%  structures/calls.
%
% I am leaning to coneCreate, rodCreate, iprgcCreate rather than
% photoreceptorCreate.  At least, I will start at the lower level and then
% maybe create the top level one when the separate ones are created.
%
% There are many types of cones, also.  I am thinking they all have the
% same slots, but the parameters differ depending on fovea, periphery, and
% so forth.  I am guessing we stick macular pigment into the parameter
% field, and optical density, and so forth, and that when we do a
% coneGet(cone,'spectral qe') we combine the whole thing.
%
%  Example:
%    cone = coneCreate('human');
%
%  See also:
%    coneGet, coneSet, sensorCreate
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check 
if notDefined('type'), species = 'human'; end
if notDefined('modelParams'), modelParams = odParams; end

species = ieParamFormat(species);

%% Check essential fields
%  Check visual field
if ~isfield(modelParams, 'visualfield')
    modelParams.visualfield = 'fovea'; 
end
%  Check wave length
if ~isfield(modelParams, 'wave'), modelParams.wave = (400:10:700)'; end

%  Check lens transmittance factor
if ~isfield(modelParams, 'lensTransF'), modelParams.lensTransF = 1; end

%  Check lens transmittance data source
if ~isfield(modelParams, 'lensTransSrc')
    modelParams.lensTransSrc = 'stockman2'; 
end

%  Check macular pigment density
if ~isfield(modelParams, 'macDensity')
    modelParams.macDensity = 0.3521; 
end

%  Check L,M,S,mel PODS
if ~isfield(modelParams, 'LPOD'), modelParams.LPOD = 0.5; end
if ~isfield(modelParams, 'MPOD'), modelParams.MPOD = 0.5; end
if ~isfield(modelParams, 'SPOD'), modelParams.SPOD = 0.4; end
if ~isfield(modelParams, 'melPOD'), modelParams.melPOD = 0.5; end

%  Check peak shift
if ~isfield(modelParams, 'peakLambda')
    % Set to default stockman and sharp S-, M- and L-cone pigment spectra
    % Stockman & Sharpe (2000), p. 1730 or http://www.cvrl.org/
    % Keep in L,M,S order
    modelParams.peakLambda = [558.9 530.3 420.7]';
end

%  Check cone density for [K,L,M,S]
if ~isfield(modelParams, 'coneDensity')
    modelParams.coneDensity  = [0.1 0.6 0.2 0.1];
end

%  Check quanta efficiency
%  The value comes from Rodieck RW, The First Steps in Seeing, page 472
if ~isfield(modelParams, 'quantalEfficiency')
    modelParams.quantalEfficiency = 0.667;
end

%  Check cone absorbtance
if ~isfield(modelParams, 'coneAbsorbtance')
    % Should generate with lambda shift, change this later
    % load T_log10coneabsorbance_ss.mat;
    % modelParams.absorbance = 10.^SplineCmf(S_log10coneabsorbance_ss,...
    %    T_log10coneabsorbance_ss, modelParams.wave,2)';
    modelParams.absorbance = StockmanSharpeNomogram(modelParams.wave, ...
                                modelParams.peakLambda)';
    % Add absorbance for K-cones 
    modelParams.absorbance = padarray(modelParams.absorbance,[0 1],'pre');
end

% Exposure time
if ~isfield(modelParams, 'expTime')
    modelParams.expTime = 0.05; % 50 ms
end

%% Create cone structure
switch species
    case 'human'
        cone.species = 'human';
        cone.type    = modelParams.visualfield;
        cone.wave    = modelParams.wave(:);
        % lens transmittance, default stockman's lens pigment densitiy
        cone.lensTrans = lensTransmittance(...
            modelParams.lensTransF, modelParams.wave, ...
            modelParams.lensTransSrc);
        cone.lensTrans = cone.lensTrans(:);
        
        % macular pigment
        cone.macular   = macular(modelParams.macDensity, ...
                                 modelParams.wave);
        
        % cone density
        cone.coneDensity = modelParams.coneDensity(:);
        
        % cone PODs
        cone.PODs = [modelParams.LPOD modelParams.MPOD ...
                     modelParams.SPOD modelParams.melPOD];
        cone.PODs = cone.PODs(:);

        % peak spectrum position
        % Only stores for L,M,S and peak for K is ignored
        cone.peakLambda = modelParams.peakLambda;
        
        % quantal efficiency
        cone.quantalEfficiency = modelParams.quantalEfficiency;
        
        % absorbance
        cone.absorbance = modelParams.absorbance;
        
        % create the underlying sensor structure
        sensor = sensorCreate('human');
        sensor = sensorSet(sensor, 'wave', modelParams.wave);
        sensor = sensorSet(sensor, 'exp time',modelParams.expTime);
        
        % set parameters to cone sensor
        sensor = sensorCreateConeMosaic(sensor, ...
                   sensorGet(sensor, 'size'), cone.coneDensity);
        sensor = sensorSet(sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
        
        cone.sensor = sensor;
        
        % set adaptation
        cone.adaptType = 2; % Adaptation by cone type, see coneAdaptation
        
    otherwise
        error('Unknown type encountered');
end


end
