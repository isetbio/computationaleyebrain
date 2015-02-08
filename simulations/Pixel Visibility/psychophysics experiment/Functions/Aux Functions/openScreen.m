function d = openScreen(d, varargin)
%% function display = openScreen(display,[varargin])
%    Open PTB window under settings in display
%
%  Inputs:
%    d        - ISETBIO compatible display structure, see displayCreate
%    varargin - name-value pair for settings
%             - supports: hideCursorFlag, bitDepth, numBuffers, drawFix
%
%  Output:
%    d        - display structure with window pointer and rect set
%
%  General Process:
%    1. Parse input parameters
%    2. If needed, test 10 bit color support
%    3. Init screen parameters to PTB: resolution, framerate, gamma, etc.
%    4. Open a PTB screen with a background color (0.5 gray by default)
%    5. Draw fixation point & Hide cursor
%  
%  See also:
%    closeScreen, drawFixation
%
% (HJ) ISETBIO TEAM, 2015

%% Check inputs
%  Check number of inputs
if nargin < 1, error('Display structure required'); end
if mod(length(varargin),2) ~= 0
    error('Parameter should be in name-value pairs');
end

%  Check fields in display structure
if ~isfield(d,'screenNumber'), d.screenNumber = 0; end
if ~isfield(d,'frameRate'),    d.frameRate = 60;   end
if ~isfield(d,'resolution')
    res = Screen('Resolution', d.screenNumber);
    d.resolution = [res.width res.height]; 
end

% check gamma table
if isempty(displayGet(d,'gamma table'))
	error('Gamma table not found in display structure');
end

% check background color
if ~isfield(d, 'backColorRgb')
    d.backColorRgb = [0.5 0.5 0.5]';
end

%% Parse varargin
bitDepth       = 8;    % color bit depth
hideCursorFlag = true; % whether to hide mouse in experiment
numBuffers     = 2;    % number of buffers
drawFix        = true; % whether to draw fixation point

for i = 1 : 2 : length(varargin)
    switch lower(varargin{i})
        case 'hidecursorflag'
            hideCursorFlag = varargin{i+1};
        case 'bitdepth'
            bitDepth = varargin{i+1};
        case 'numbuffers'
            numBuffers = varargin{i+1};
        case 'drawfix'
            drawFix  = varargin{i+1};
    end
end

% check bitDepths, only support 8 bit and 10 bit
assert(bitDepth == 8 || bitDepth == 10,'Error: Unknown bit depth');

%% Test 10 bit support
if bitDepth == 10
    try
        AdditiveBlendingForLinearSuperpositionTutorial('Native10Bit');
    catch
        disp('10 bit not supported on this machine, use 8 bit instead');
        bitDepth = 8;
    end
end

%% Init screen parameters
% Skip flickering warning
Screen('Preference', 'SkipSyncTests', 0);

% Set the resolution
try
    % Try to set spatial resolution, then spatial and temporal
    Screen('Resolution', d.screenNumber, d.resolution(1), d.resolution(2));
    Screen('Resolution', d.screenNumber, ...
        d.resolution(1), d.resolution(2), d.frameRate);
catch ME
    warning(ME.identifier, ME.message)
end

WaitSecs(2);

% Save current gamma table
d.oldGamma=Screen('ReadNormalizedGammaTable', d.screenNumber);
try
    %Screen('LoadNormalizedGammaTable', display.screenNumber,display.gamma);
    % Load linearlized gamma table
    invGamma = displayGet(d, 'inverse gamma table');
    invGamma = invGamma / max(invGamma(:));
    Screen('LoadNormalizedGammaTable', d.screenNumber, invGamma);
    %LoadIdentityClut(display.windowPtr);
catch ME
    warning(ME.identifier, ME.message)
    % 10 bit gamma table not supported, reduce it to 8 bit
    pGamma = d.gamma(round(linspace(1,size(d.gamma,1),256)),:);
    Screen('LoadNormalizedGammaTable', d.screenNumber,pGamma);
end

%% Open PTB Screen
if bitDepth == 8 % Open screen for 8 bit
    [d.windowPtr,d.rect] = Screen('OpenWindow', ...
        d.screenNumber, d.backColorRgb, [],[], numBuffers);
else % Open screen for 10 bit
    PsychImaging('PrepareConfiguration');
	PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    PsychImaging('AddTask', 'General', 'EnableNative10BitFrameBuffer');
    [d.windowPtr,d.rect] = PsychImaging('OpenWindow', ...
        d.screenNumber, d.backColorRgb, [],[], numBuffers);
end

%% Handle cursor and fixation
%  Draw fixation point and hide cursor if necessary
if drawFix, drawFixation(d); end
if(hideCursorFlag), HideCursor;  end

% Flip and show
Screen('Flip', d.windowPtr);
%LoadIdentityClut(display.windowPtr);

end
