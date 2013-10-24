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
%    {'peak shift', 'lambda shift'}  - peak shift, varargin is used to
%                                      decide units, default 'nm'
%    {'qe', 'quantal eff'}           - quantal efficiency
%    {'absorbance'}                  - cone absorbance
%    {'absorbtance'}                 - cone absorbtance without any
%                                      pre-receprocal transimitance
%    {'effetive absorbtance'}        - cone absorbtance with lens and
%                                      macular pigment transmittance
%
%    MORE PARAMETERS ABOUT UNDERLYING SENSOR CAN BE FOUND IN sensorGet
%    FUNCTION
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
    case {'species', 'kind'}
        val = cone.species;
    case {'density', 'conedensity'}
        val = cone.coneDensity;
    case {'wave', 'wavelength'}
        val = cone.wave;
    case {'type', 'visualfield'}
        val = cone.type;
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
    case {'peakshift', 'lambdashift'}
        val = cone.peakShift;
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
    otherwise
        % Try to get parameter value from the underlying sensor
        val = sensorGet(cone.sensor, param, varargin);
end

end
