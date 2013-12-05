function cone = coneImport(filename, type, fillData, varargin)
% Import cone structure from plain text data
%   cone = coneImport(filename, [type], [fillData], [varargin])
%
% This function is used to import cone structure from some readable plain
% text file, or directly from a string
%
% Inputs:
%   filename  - name of file to be read from. If filename is a double, not
%               string, we treat it as fid
%   type      - now only 'json' or 'string' format is supported, when type
%               is 'string', the filename contains string of plain text
%               Default: 'json' 
%   fillData  - bool, indicating whether or not to fill data field with
%               default value after importing the cone structure. If not,
%               missing data field will be set to empty
%   varargin  - key-value pairs for overwriting field in cone structure. If
%               not set here, user can also set it afterwards with coneSet
%               function
%
% Outputs:
%   cone      - cone structure, see coneCreate for more information
%
% Example:
%   cone = coneExport('defaultCone.json', 'json', true);
%
% See also:
%   coneExport, coneCreate, coneSet
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('filename'), error('filename / plain text required'); end
if notDefined('type'), type = 'json'; end
if notDefined('fillData'), fillData = true; end

if mod(length(varargin),2)~=0, error('varargin should be in pairs'); end

%% Init
%  open file for reading
if ~isempty(filename) && ~strcmp(type, 'string')
    if ischar(filename)
        fid = fopen(filename, 'r');
    else
        fid = filename;
    end
    assert(fid > 0, 'file cannot be opened');
end

%% Generate cone structure from plain text
switch type
    case 'json'
        cone = json.load(fscanf(fid, '%c', inf));
        cone = cone.data; % Get the data field and schema is ignored
    case 'string'
        cone = json.load(filename);
        cone = cone.data;
    otherwise
        error('Unknown output format type encountered');
end

%% Fill in missing field
%  fill default value
if fillData
    % build a default cone. we will get value from here
    defaultCone = coneCreate;
    defaultM = coneGet(defaultCone, 'macular');
    defaultLens = coneGet(defaultCone, 'lens');
    
    m = coneGet(cone, 'macular');
    lens = coneGet(cone, 'lens');
    
    % set fields to macular
    % Since coneSet assume we have a complete cone structure, we could not
    % use it directly here. Sad. HJ
    if ~isfield(m, 'unit density')
        m.unitDensity = macularGet(defaultM, 'unit density');
    end
    if ~isfield(m, 'wave')
        m.wave = macularGet(defaultM, 'wave');
    end
    
        
    % set fields to lens
    if ~isfield(lens, 'unit density')
        lens.unitDensity = macularGet(defaultLens, 'unit density');
    end
    if ~isfield(lens, 'wave')
        lens.wave = macularGet(defaultLens, 'wave');
    end
    
    % set fields to cone
    if ~isfield(cone, 'absorbance')
        cone.absorbance = coneGet(defaultCone, 'absorbance');
    end
    if ~isfield(cone, 'wave')
        cone.wave = coneGet(defaultCone, 'wave');
    end
    cone = coneSet(cone, 'macular', m);
    cone = coneSet(cone, 'lens', lens);
    
else
    % Just set to empty if missing
    % Since coneSet assume we have a complete cone structure, we could not
    % use it directly here. Sad. HJ
    if ~isfield(cone.macular, 'unit density')
        cone.macular.unitDensity = [];
    end
    if ~isfield(cone.macular, 'wave')
        cone.macular.wave = [];
    end
    
        
    % set fields to lens
    if ~isfield(cone.lens, 'unit density')
        cone.lens.unitDensity = [];
    end
    if ~isfield(cone.lens, 'wave')
        cone.lens.wave = [];
    end
    
    % set fields to cone
    if ~isfield(cone, 'absorbance')
        cone.absorbance = [];
    end
    if ~isfield(cone, 'wave')
        cone.wave = [];
    end
end

%  now we should have a valid cone structure for coneGet and coneSet
%  parse varargin and overwrite fields
for ii = 1 : 2 : length(varargin)
    cone = coneSet(cone, varargin{ii}, varargin{ii+1}); 
end