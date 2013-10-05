function cone = coneCreate(type, varargin)
%% function coneCreate
%  Default cone structure
%
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

%% Check 
if notDefined('type'), type = 'human'; end
if nargin < 2, density = [.1 .6 .2 .1]; else density = varargin{1}; end

%% Create cone structure
switch type
    case 'human'
        cone.species = 'human';
        sensor = sensorCreate('human');
        sensor = sensorSet(sensor, 'exp time',0.05); % Init to 50 ms
        sensor = sensorSet(sensor, 'fov', 2);        % Init to 2 degree
        [sensor, xy, coneType] = sensorCreateConeMosaic(sensor, ...
                          sensorGet(sensor, 'size'), density);
        cone.sensor   = sensor;
        cone.conePos  = xy;
        cone.coneType = coneType;
        cone.density  = density; 

    otherwise
        error('Unknown type encountered');
end


end
