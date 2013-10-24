function cone = coneSet(cone, param, val, varargin)
%% function coneSet(cone, param, val, [varargin])
%    Set parameters to cone structure. Please refer to coneGet for more
%    detailed information about supported parameters
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
%    COMMENT TO BE ADDED HERE
%
%  Example:
%    cone = coneCreate('human');
%    cone = coneSet('exp time', 0.02);
%
%  See also:
%    coneCreate, coneGet
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
                   sensorGet(sensor, 'size'), cone.coneDensity);
               
    case {'wave', 'wavelength'}
        cone.wave = val;
        cone.sensor  = sensorSet(cone.sensor, 'wave', val);
        cone.macular = macular(coneGet(cone, 'mac density'), val);

    case {'lenstrans', 'lenstransmittance'}
        cone.lensTrans = val;
        cone.sensor = sensorSet(cone.sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
        
    case {'lensabsorption'}
        cone = coneSet(cone, 'lens trans', 1 - val);
        cone.sensor = sensorSet(cone.sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
        
    case {'macular', 'macularpigments'}
        cone.macular = val;
        cone.sensor = sensorSet(cone.sensor, 'color filters', ...
            coneGet(cone, 'eff absorbtance'));
        
    case {'macdens','maculardensity'}
        cone.macular = macular(val, coneGet(cone, 'wave'));

    case {'maculartrans', 'mactrans'}
        if (any(size(cone.macular.transmittance)~= size(val)))
            error('Transmittance size mismatch');
        end
        cone.macular.transmittance = val;
        cone.macular.absorption = 1 - val;
        
    case {'macularabsorption'}
        if (any(size(cone.macular.absorption)~= size(val)))
            error('Absorption size mismatch');
        end
        cone.macular.absorption = val;
        cone.macular.transmittance = 1 - val;

    case {'pods','pod'}
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
    case {'peakshift', 'lambdashift'}
        % Should re-generate a lot of things, think about it later
        cone.peakShift = val;
        
    case {'qe', 'quantalefficiency', 'quantaleff'}
        cone.quantalEfficiency = val;
        
    case {'absorbance'}
        disp('Overwrite absorbtance, please be sure what you are doing');
        cone.absorbance = val;
        
    otherwise % Try set the sensor parameter
        cone.sensor = sensorSet(cone.sensor,param,val,varargin);
end

end
