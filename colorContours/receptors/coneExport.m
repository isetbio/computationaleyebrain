function txtStr = coneExport(cone, filename, type, clearData, varargin)
% Export cone structure to plain text data format
%   txtStr = coneExport(cone, [filename], [type], [clearData], [varargin])
%
% This function is used to export cone structure to some readable plain
% text file. The plain text file would be used in parameter editing and it
% could be loaded back by calling coneImport function
%
% Inputs:
%   cone      - cone structure, created by coneCreate
%   filename  - name of file to be written into, if empty, the function
%               will not export to files. If filename is a double, not
%               string, we treat it as fid
%   type      - now only 'json' format is supported
%   clearData - bool, indicating whether or not to clear data field bofore
%               exporting the cone structure
%   varargin  - key-value pairs for overwriting comment and display formats
%               of fields (this could help improving the appearance of the
%               json based webpage)
%
% Outputs:
%   txtStr    - string, containing the cone parameters in plain text
%
% Example:
%   cone = coneCreate;
%   jsonStr = coneExport(cone, 'defaultCone.json', 'json');
%
% See also:
%   coneCreate, coneImport, coneClearData
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if notDefined('cone'), error('cone structure required'); end
if notDefined('filename'), filename = []; end
if notDefined('type'), type = 'json'; end
if notDefined('clearData'), clearData = true; end
if mod(length(varargin),2)~=0, error('varargin should be in pairs'); end

assert(strcmp(cone.type, 'cone'), 'first input should be cone structure');

%% Init
%  open file for writing
if ~isempty(filename)
    if ischar(filename)
        fid = fopen(filename, 'w');
    else
        fid = filename;
    end
    assert(fid > 0, 'file cannot be opened');
end

%  clear data fields if needed
if clearData, cone = coneClearData(cone); end

%% Export data field
switch type
    case 'json' % output as json format
        % Init json
        json.startup;
        txtStr = json.dump(cone);
    otherwise
        error('Unknown output format type encountered');
end

%% Generate json scheme
schema = [];
schemaStr = '';

% Will add something for schemaStr soon. HJ
% This is only useful when loading json to web-browser, and it will only
% affect the appearance there

for ii = 1 : 2 : length(varargin)
    schema.(varargin{ii}) = varargin{ii+1};
end

if ~isempty(schema), schemaStr = json.dump(schema); end

%% Write to file
if ~isempty(filename)
    fprintf(fid, '{"data":%s}', txtStr);
    if ~isempty(schemaStr)
        fprintf(fid, '{"schema":%s}', schemaStr);
    end
    fclose(fid);
end

end