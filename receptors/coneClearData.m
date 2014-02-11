function cone = coneClearData(cone, varargin)
% Clear data field and only retain parameters in cone structure
%
%      cone = coneClearData(cone, [varargin])
%
% This function is mainly used to clear data field in cone structure.
% Generally, it will be used before export cone parameters to text files,
% e.g. json files. Clearing data before export can help make the json files
% smaller and more readable
%
% Inputs:
%   cone      - cone structure, created by coneCreate
%   varargin  - key-value pairs to set special field
%
% Outputs:
%   cone      - cone with empty data field
%   
% Example:
%   cone = coneCreate;
%   cone = coneClearData;
%
% See also:
%   coneCreate, coneExport
%
% HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check inputs
if nargin < 1, error('cone structure required'); end
if mod(length(varargin),2)~=0, error('varargin should be in pairs'); end
assert(strcmp(cone.type, 'cone'), 'first input should be cone structure');

%% Clear data field
%  remove data field from cone
if isfield(cone, 'wave'), cone = rmfield(cone, 'wave'); end
if isfield(cone, 'absorbance'), cone = rmfield(cone, 'absorbance'); end

%  remove data field from macular
cone.macular = rmfield(cone.macular, 'wave');
cone.macular = rmfield(cone.macular, 'unitDensity');

%  remove data field from lens
cone.lens = rmfield(cone.lens, 'wave');
cone.lens = rmfield(cone.lens, 'unitDensity');


%% Overwrite field with varargin
for ii = 1 : 2 : length(varargin)
    cone.(varargin{ii}) = varargin{ii + 1};
end

end