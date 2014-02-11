function val = coneGet(cone, param, varargin)
%% function coneGet(cone, params, [varargin])
%    Get properties from cones
%
%  Inputs:
%    cone     - cone structure, created by coneCreate
%    param    - parameter name to get, spaces and case insensitive
%    varargin - possible units for some parameters
%
%  Outputs:
%    val      - value for parameter, if not found, return empty
%
%  Supported params:
%    {'name'}                        - user defined name of the cone
%    {'type'}                        - should be 'cone'
%    {'species', 'kind'}             - cone species, generally 'human'
%    {'density', 'cone density'}     - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'wave', 'wavelength'}          - wavelength of samples in cones
%    {'lens'}                        - underlying lens structure
%    {'lens transmittance'}          - lens transmittance
%    {'lens absorption'}             - lens absorbtance
%    {'macular', 'macular pigments'} - macular pigment structure
%    {'mac dens','macular density'}  - macular density
%    {'macular transmittance'}       - macular transmittance
%    {'macular absorption'}          - macular absorption
%    {'eye trans'}                   - totally transmittance for lens and
%                                      macular pigments
%    {'PODs','POD'}                  - pigment density vector for [L,M,S]
%    {'LPOD'}                        - L pigment density
%    {'MPOD'}                        - M pigment density
%    {'SPOD'}                        - S pigment density
%    {'peak efficiency'}             - quantal efficiency
%    {'absorbance'}                  - cone absorbance
%    {'absorptance'}                 - cone absorbtance without any
%                                      pre-receprocal transimitance
%    {'effetive absorptance'}        - cone absorbtance with lens and
%                                      macular pigment transmittance
%    {'quantal fundamentals'}        - quantal fundamentals of the cones
%
%
%  Example:
%    cone = coneCreate('human');
%    expTime = coneGet(cone, 'species');
%
%  See also:
%    coneSet, coneCreate
%
%  TODO:
%    For most parameters, we should accept a third parameter as wavelength
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('Cone structure required'); end
if notDefined('param'), error('Parameter name required'); end

%% Get property value
param = ieParamFormat(param);  % Lower case and remove spaces
switch param   
    case {'name'}
        val = cone.name;
    case {'type', 'visualfield'}
        val = cone.type; % Should always be 'cone'
    
    case {'wave', 'wavelength'}
        val = cone.wave;
    case {'species', 'kind'}
        % Animal species.
        % Currently supports only 'human'
        val = cone.species;
        
    case {'spatialdensity', 'conespatialdensity'}
        % Spatial density of the cone samples.
        % The vector should be in form [K, L, M, S]
        % Sum of the four elememts should be 1
        val = cone.spatialDensity;
 
    case {'lens'}
        % Lens structure, see lensCreate and lensGet
        val = cone.lens;
    case {'lenstrans', 'lenstransmittance'}
        % Lens transmittance
        val = lensGet(coneGet(cone,'lens'), ' transmittance');
    case {'lensabsorption','lensabsorbtance'}
        val = lensGet(coneGet(cone,'lens'), 'absorption');
        
    case {'macular', 'macularpigment'}
        % coneGet(cone,'macular pigment')
        val = cone.macular;

    case {'macdens','maculardensity'}
        val = macularGet(coneGet(cone,'macular'), 'density');
    case {'maculartransmittance','maculartrans', 'mactrans'}
        % coneGet(cone,'macular transmittance')
        val = macularGet(coneGet(cone,'macular'), 'transmittance');
    case {'macularabsorption'}
        % coneGet(cone,'macular absorption')
        val = macularGet(coneGet(cone,'macular'), 'absorption');
        
    case {'oculartransmittance','eyetrans','eyetransmittance'}
        % coneGet(cone,'ocular transmittance')
        % pre-retina eye transimittance, including lens and macular pigment
        val = coneGet(cone,'lens trans') .* ...
                coneGet(cone, 'macular trans');
            
    case {'pods','pod'}
        % pigment density inside each cone cell
        % 3-element vector, for [L,M,S] repectively
        val = cone.opticalDensity;
    case {'lpod'}
        val = cone.opticalDensity(1);
    case {'mpod'}
        val = cone.opticalDensity(2);
    case {'spod'}
        val = cone.opticalDensity(3);
        
        %     case {'peaklambda', 'lambdamax'}
        %         % Position with max cone absorption
        %         % This can be used to simulate color anormalous
        %         val = cone.peakLambda;
        
    case {'peakefficiency'}
        val = cone.peakEfficiency;
        
    case {'absorbance'}
        val = cone.absorbance;
        
    case {'conespectralabsorptance','absorptance'}
        % coneGet(cone,'cone spectral absorptance')
        %
        % This is the cone absorptance without the ocular media
        absorbance = coneGet(cone, 'absorbance');
        PODs = coneGet(cone, 'PODs');
        val = 1-10.^(-absorbance*diag(PODs));
        
    case {'effectivespectralabsorptance', 'effectiveabsorptance'}
        % coneGet(cone,'effective spectral absorptance')
        %
        % Combines cone photopigment, ocular transmittance and peak
        % efficiency
        absorptance = coneGet(cone, 'absorptance');
        eyeTrans    = coneGet(cone, 'ocular transmittance');
        peakEfficiency = coneGet(cone,'peak efficiency');
        
        val = (absorptance.*repmat(eyeTrans, [1 size(absorptance,2)]))* ...
            diag(peakEfficiency);
        
    case {'quantalfundamentals','photonfundamentals'}
        % coneGet(cone,'photon fundamentals')
        %
        % Cone absorptance scaled by peak efficiency
        % Note that this is for pure cones, without any oclus absorptance
        val = coneGet(cone, 'cone spectral absorptance');
        qe  = coneGet(cone, 'peak efficiency');
        if length(qe) == size(val,2)
            for ii = 1 : size(val, 2)
                val(:,ii) = val(:,ii) * qe(ii);
            end
        end
        val = val ./ repmat(max(val), [size(val, 1) 1]);
        
    case {'energyfundamentals'}
        % coneGet(cone, 'energy fundamentals')
        %
        % cone absorptance in energy
        % Note that this is for pure cones, without any oclus absorptance
        val = coneGet(cone, 'quantal fundamentals');
        % Get constants
        h = vcConstants('planck');
        c = vcConstants('speed of light');
        wave = coneGet(cone, 'wave');
        % Convert to energy. Note that the sensitivity is the inverse of
        % spetral absorption. Thus, we use conversion of energy to quanta
        % on quantal fundamentals
        val = val / (h*c) .* (1e-9 * repmat(wave, [1 size(val, 2)])); 
        
        % Normalize
        val = val ./ repmat(max(val), [size(val, 1) 1]);
        
    otherwise
        error('Unknown parameter encountered');
end

end
