function cone = coneSet(cone, param, val, varargin)
% Set parameters in a cone structure.
%
%      cone = coneSet(cone, param, val, [varargin])
%
%  Inputs:
%    cone     - cone structure, created by coneCreate
%    param    - parameter name to set, spaces and case insensitive
%    val      - new value of parameter to be set to
%    varargin - possible units for some parameters
%
%  Outputs:
%    cone     - cone structure with parameter and its related information
%               set to new value
%
%  Supported parameters:
%    {'density', 'cone density'}     - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'fov', 'h fov'}                - sensor fov, scene and oi are
%                                      accepted in varargin
%    {'wave', 'wavelength'}          - wavelength of samples in cones
%    {'sensor'}                      - underlying ISET sensor structure
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
%    {'absorbance'}                  - cone absorbance, not recommended to
%                                      set this parameter directly unless
%                                      you know what you are doing
%    {'adaptation type','adapt type} - cone adaptation type, see
%                                      coneAdaptation for more details
%
%    MORE SUPPORTED PARAMETERS CAN BE FOUND IN FUNCTION sensorSet
%
%  Example:
%    cone = coneCreate('human');
%    cone = coneSet('density',[0.1 0.65 0.2 0.05]);
%    cone = coneSet('exp time', 0.02);
%
%  See also:
%    coneCreate, coneGet, sensorSet
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('cone structure required'); end
if notDefined('param'), error('parameter name required'); end
if notDefined('val'), error('new value of parameter required'); end

%% Set parameters
param = ieParamFormat(param);  % Lower case and remove spaces
switch param
    case {'density', 'conedensity'}
        cone.coneDensity = val;
        cone.sensor = sensorCreateConeMosaic(cone.sensor, ...
                   sensorGet(cone.sensor, 'size'), cone.coneDensity);
               
    case {'wave', 'wavelength'}
        cone.wave = val;
        cone.sensor  = sensorSet(cone.sensor, 'wave', val);
        cone.macular = macular(coneGet(cone, 'mac density'), val);
        
    case {'sensor'}
        cone.sensor = val;
        
    case {'fov', 'hfov'}
        if ~isempty(varargin), scene = varargin{1}; end
        if length(varargin) > 1, oi = varargin{2}; end
        cone.sensor = sensorSetSizeToFOV(cone.sensor, val, scene, oi);

    case {'lenstrans', 'lenstransmittance'}
        cone.lensTrans = val;
        cone.sensor = sensorSet(cone.sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
        
    case {'lensabsorption'}
        %
        cone = coneSet(cone, 'lens trans', 1 - val);
        cone.sensor = sensorSet(cone.sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
    
        % Macular pigment related sets
    case {'macular'}
        % cone = coneSet(cone,'macular',m);
        cone.macular = val;
        
    case {'maculardensity'}
        % cone = coneSet(cone,'macular density',val)
        % val is typically between 0 and 0.7, a range of macular pigment
        % densities.
        %
        m    = coneGet(cone,'macular');
        m    = macularSet(m,'density',val);
        cone = coneSet(cone,'macular',m);
        
    case {'pods','pod'}
        % Pigment optical densities for the cones and melanopsin 
        if (any(size(cone.PODs)~= size(val)))
            error('PODs size mismatch');
        end
        cone.PODs = val;
        
    case {'lpod'}
        if ~isscalar(val), error('val should be scaler'); end
        cone.PODS(1) = val;
    case {'mpod'}
        if ~isscalar(val), error('val should be scaler'); end
        cone.PODS(2) = val;
    case {'spod'}
        if ~isscalar(val), error('val should be scaler'); end
        cone.PODS(3) = val;
    case {'melpod'}
        if ~isscalar(val), error('val should be scaler'); end
        cone.PODS(4) = val;
        
    case {'peaklambda', 'lambdamax'}
        cone.peakShift = val;
        cone.absorbance = StockmanSharpeNomogram(cone.wave, val)';
        cone.absorbance = padarray(cone.absorbance,[0 1],'pre');
        
    case {'qe', 'quantalefficiency', 'quantaleff'}
        % Is this spectral?
        cone.quantalEfficiency = val;
        
    case {'absorbance'}
        disp('Overwrite absorbtance, please be sure what you are doing');
        cone.absorbance = val;
    
    case {'adaptationtype', 'adapttype'}
        if cone.adaptType ~= val
            cone = coneAdapt(cone, val);
        end
        
    otherwise % Try set the sensor parameter
        cone.sensor = sensorSet(cone.sensor,param,val,varargin);
end

end
