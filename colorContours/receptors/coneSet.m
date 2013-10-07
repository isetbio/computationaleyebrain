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
    case {'position', 'pos', 'conepos'}
        assert(coneGet(cone, 'size') == size(val), 'pos size mismatch');
        cone.conePos = val;
    case {'density', 'conedensity'}
        cone.density = val(:);
        % Adjust cone pattern to new density
        % Shall I keep the rSeed the same as before?
        [sensor, xy, coneType] = sensorCreateConeMosaic(cone.sensor, ...
                          sensorGet(cone.sensor, 'size'), cone.density);
        cone.sensor   = sensor;
        cone.conePos  = xy;
        cone.coneType = coneType;
        
    case {'conetype', 'pattern'}
        % I shall think about how to do this one
        
    % Should add more parameter about second site noise here
    otherwise % Try set the sensor parameter
        cone.sensor = sensorSet(cone.sensor,param,val,varargin);
end

end
