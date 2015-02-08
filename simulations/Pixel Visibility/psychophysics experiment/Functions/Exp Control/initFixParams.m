function d = initFixParams(d, fov, varargin)
%% function display = initFixParams(display, fov, [varargin])
%    Initialize display fixation point parameters
%
%  Inputs:
%    d        - ISETBIO compatible display structure
%    fov      - Fixation point size in degrees
%    varargin - Name-value pair for fixation parameters
%             - Acceptable name: 'Position','fixType', 'fixColor'
%
%  Outputs:
%    d        - display structure with fixation point information set
%
%  Example:
%    d = initDisplay('LCD-Gechic.mat');
%    d = initFixParams(d, 1);
%
%  See also:
%    initDisplay
%
%  (HJ) ISETBIO TEAM, 2015

%% Check inputs
if nargin < 1, error('Display structure required'); end
if nargin < 2, error('Size of fixation point required (in degree)'); end
if mod(length(varargin),2)~=0, error('Parameters should be in pairs'); end

%% Init fixation by default value
d.fixType       = 'dot';
d.fixSizePixels = 2;        % number of pixels
d.fixColorRgb   = [0 0 0]'; % black
              
ecc = angle2pix(d, fov);
d.fixStim = round([0 -1 1] * ecc + d.numPixels(1)/2); 
d.fixPosY = round(d.numPixels(2)/2);
d.fixPosX = round(d.numPixels(1)/2);

% If using bits++ in color++ mode, we should further divide this by 2
% d.fixPosX = d.fixPosX / 2;

%% Set Parameter by varargin
for i = 1 : 2 : length(varargin)
    switch lower(varargin{i})
        case 'position'
            pos = varargin{i+1};
            d.fixPosY = round(pos(2)/2);
            d.fixPosX = round(pos(1)/2);
        case 'fixtype'
            d.fixType = varargin{i+1};
        case 'fixcolor'
            d.fixColorRgb = varargin{i+1};
        otherwise
            warning(['Unknown parameter ' varargin{i} ', ignored']);
    end
end

%% END