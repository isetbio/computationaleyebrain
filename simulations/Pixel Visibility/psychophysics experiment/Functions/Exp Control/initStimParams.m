function stimParams = initStimParams(varargin)
%% function stimParams = initStimParams([varargin])
%    initialize stimulus parameters
%
%  Input:
%    varargin   - name value pair, could be anything you want to be saved
%                 in stimParams
%  Output:
%    stimParams - stimulus structure
%
%  Example:
%    stimParams = initStimParams('Gsig',10);
%  
% (HJ) ISETBIO TEAM, 2015

%% Check inputs
if mod(length(varargin),2)~=0, 
    error('Inputs should be in name-value pairs'); 
end

%% Initialize default values
%  Set size of stimulus in degrees
stimParams.visualSize = [2 0.5]; %(deg)

% gapSize
stimParams.gapSize = 1/5;

% Spacial Blur
stimParams.Gsig = 0;

% duration of stimulus presentation
stimParams.duration = 1; % seconds

% init isi
stimParams.isi = 0.01;

%% Parse varargin and set to stimParams
for i = 1 : 2 : length(varargin)
    stimParams.(varargin{i}) = varargin{i+1};
end

%% END