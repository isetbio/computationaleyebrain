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
%    val      - value for parameter
%
%  Supported params:
%    {'position', 'pos', 'cone pos'} - x,y position of each cone
%    {'density', 'cone density'}     - density of each cone type, for
%                                      human, it should be [K,L,M,S]
%    {'cone type', 'pattern'}        - type of each cone, should have same
%                                      size as sensor size of the cone
%    {'size', 'sz'}                  - cone sensor size
%    {'exp time', 'exposure time'}   - eye integration time
%    {'fov', 'h fov'}                - horizontal field of view
%    OTHER SUPPORTED PARAMETERS CAN BE FOUND IN 'sensorGet' FUNCTION
%
%  Example:
%    cone = coneCreate('human');
%    expTime = coneGet(cone, 'exp time');
%
%  See also:
%    coneSet, coneCreate
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('Cone structure required'); end
if notDefined('params'), error('Parameter name required'); end

%% Get property value
param = ieParamFormat(param);  % Lower case and remove spaces
switch param
    case {'position', 'pos', 'conepos'}
        val = cone.conePos;
    case {'density', 'conedensity'}
        val = cone.density;
    case {'conetype', 'pattern'}
        val = cone.coneType;
    otherwise
        val = sensorGet(cone.sensor, param, varargin);
end

end
