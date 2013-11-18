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
%    {'species', 'kind'}             - cone species, generally 'human'
%    {'density', 'cone density'}     - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'sensor'}                      - underlying sensor structure
%    {'wave', 'wavelength'}          - wavelength of samples in cones
%    {'type', 'visual field'}        - visual field of eye, can be 'fovea'
%                                      or 'periphery'
%    {'lens trans'}                  - lens transmittance
%    {'lens absorption'}             - lens absorbtance
%    {'macular', 'macular pigments'} - macular pigment structure
%    {'mac dens','macular density'}  - macular density
%    {'macular trans'}               - macular transmittance
%    {'macular absorption'}          - macular absorption
%    {'eye trans'}                   - totally transmittance for lens and
%                                      macular pigments
%    {'PODs','POD'}                  - PODs vector for [L,M,S,mel]
%    {'LPOD'}                        - L POD density
%    {'MPOD'}                        - M POD density
%    {'SPOD'}                        - S POD density
%    {'melPOD'}                      - mel POD density
%    {'peak lambda', 'lambda max'}   - peak spectra position
%    {'qe', 'quantal eff'}           - quantal efficiency
%    {'absorbance'}                  - cone absorbance
%    {'absorbtance'}                 - cone absorbtance without any
%                                      pre-receprocal transimitance
%    {'effetive absorbtance'}        - cone absorbtance with lens and
%                                      macular pigment transmittance
%    {'quantal fundamentals'}        - quantal fundamentals of the cones
%    {'adapted volts'}               - volts image after cone adaptation,
%                                      which could be used for rgc
%                                      computation
%    {'adapt gain','gain map'}       - cone adaptation gain map
%    {'adaptation type','adapt type} - cone adaptation type, see
%                                      coneAdaptation for more details
%    {'adapt offset'}                - cone adaptation offset
%
%    MORE PARAMETERS ABOUT UNDERLYING SENSOR CAN BE FOUND IN sensorGet
%    FUNCTION
%
%  Example:
%    cone = coneCreate('human');
%    expTime = coneGet(cone, 'species');
%
%  See also:
%    coneSet, coneCreate, sensorGet
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
        val = cone.type;
        case {'wave', 'wavelength'}
        val = cone.wave;
    case {'species', 'kind'}
        % Which animal species.
        % Currently supports only 'human'
        val = cone.species;
        
    case {'density', 'conedensity'}
        % Spatial density of the cone samples.
        % The sum of this vector should be 1
        % The ...
        val = cone.coneDensity;
    case {'sensor'}
        % Not sure this should be here
        val = cone.sensor;
        

    case {'lenstrans', 'lenstransmittance'}
        val = cone.lensTrans;
    case {'lensabsorption'}
        val = 1 - cone.lensTrans;
        
    case {'macular', 'macularpigments'}
        val = cone.macular;
    case {'macdens','maculardensity'}
        val = cone.macular.density;
    case {'maculartrans', 'mactrans'}
        val = cone.macular.transmittance;
    case {'macularabsorption'}
        val = cone.macular.absorption;
    case {'eyetrans','eyetransmittance'}
        val = cone.lensTrans .* cone.macular.transmittance;
    case {'pods','pod'}
        val = cone.PODs;
    case {'lpod'}
        val = cone.PODS(1);
    case {'mpod'}
        val = cone.PODS(2);
    case {'spod'}
        val = cone.PODS(3);
    case {'melpod'}
        val = cone.PODS(4);
    case {'peaklambda', 'lambdamax'}
        val = cone.peakLambda;
    case {'qe', 'quantalefficiency', 'quantaleff'}
        val = cone.quantalEfficiency;
    case {'absorbance'}
        val = cone.absorbance;
    case {'absorbtance'}
        absorbance = coneGet(cone, 'absorbance');
        wave = coneGet(cone, 'wave');
        PODs = coneGet(cone, 'PODs');
        PODs = [PODs(4);PODs(1:3)]; % Change to K, L, M, S order
        val = AbsorbanceToAbsorbtance(absorbance', wave, PODs)';
    case {'effetive absorbtance', 'effabsorbtance'}
        absorbtance = coneGet(cone, 'absorbtance');
        eyeTrans = coneGet(cone, 'eye trans');
        val = absorbtance .* repmat(eyeTrans, [1 size(absorbtance, 2)]);
    case {'quantalfundamentals'}
        val = coneGet(cone, 'eff absorbtance');
        qe  = coneGet(cone, 'qe');
        if length(qe) == size(val,2)
            for i = 1 : size(val, 2)
                val(:,i) = val(:,i) * qe(i);
            end
        end
        val = val ./ repmat(max(val), size(val, 2));
    case {'adaptedvolts'}
        if isfield(cone, 'adaptVolts')
            val = cone.adaptVolts;
        else
            val = [];
            warning('Adaptation image not computed');
        end
    case {'adaptgain','gainmap'}
        if isfield(cone, 'adaptGain')
            val = cone.adaptGain;
        else
            val = [];
            warning('Adaptation image not computed');
        end
    case {'adaptationtype', 'adapttype'}
        val = cone.adaptType;
    case {'adaptoffset'}
        if isfield(cone, 'adaptOffset')
            val = cone.adaptOffset;
        else
            val = [];
            warning('Adaptation image not computed');
        end
        
    otherwise
        % Try to get parameter value from the underlying sensor
        if isempty(varargin)
            val = sensorGet(cone.sensor, param);
        else
            val = sensorGet(cone.sensor, param, varargin);
        end
end

end
