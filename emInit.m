function [sensor, params] = emInit(emType, sensor, params)
% Init eye movement parameters in the sensor structure
%
%   [sensor, params] = emInit(emType, sensor, params)
%
% emType:  Eye movement type (fixation, ....)
% sensor:  The sensor
% params:  Depends on type
%   fixation:  sdx, sdy, center, nSamples, randSeed (optional).  The sd
%   units are in deg of visual angle.
%
% General Process:
%   1. Check eyeMoveType and set random seed
%   2. Generate moving postion x and y
%   3. Generate frames per position according to distribution and nSamples
%   4. Generate linear eye movement for testing
%
%
% Output Parameter:
%   sensor       - sensor with eye movement related parameters set
%
% Example:
%    sensor = emInit(sensor, scene, oi, 100, 1)
%
% (HJ) Copyright PDCSOFT TEAM 2013

%% Check inputs and Init
if ieNotDefined('emType'), error('eye movement type required'); end
if ieNotDefined('sensor'), error('sensor required.'); end
if ieNotDefined('params'), error('parameters required.'); end

% Initialize random number generation
if isfield(params,'randSeed'), rng(randSeed);
else
    % Figure out how to return the random number condition that will
    % produce the same result.  Stick it in params.randSeed.
    rng('shuffle');
end

% Each case builds the (x,y) and count variables for every position
emType = ieParamFormat(emType);
switch emType
    case {'fixation'}
        
        % Algorithm description to be placed here.  Probably should change
        % the algorithm for deltaX, deltaY steps and thus continuity of the
        % path.
        
        % Init parameters for Gaussian RV.  Should check for existence.
        center   = params.center;
        sigmaX   = params.sigmaX;
        sigmaY   = params.sigmaY;
        nSamples = params.nSamples;
        fov      = params.fov;
        
        % Compute sensor fov and size
        sz  = sensorGet(sensor, 'size');
        
        %% Initialize eye movements
        
        % This should be a case, we might call it 'fixating'.
        %
        % The eye wanders around the center. The positions are random in a
        % disk around the center. the distances here are in deg of visual
        % angle.
        x    = randn(nSamples,1)*sigmaX + center(1);
        y    = randn(nSamples,1)*sigmaY + center(2);
        
        % For efficiency, we round the calculations
        % that are centered less than 1 detector's width.
        xPos = round((x/fov)*sz(1)) * fov/sz(1);
        yPos = round((y/fov)*sz(2)) * fov/sz(2);
        
        % Group the same positions
        [~,~,ic]  = unique([xPos yPos],'rows');
        
        % Compute frame per position
        f    = hist(ic,unique(ic));      % frames per position
        f(1) = f(1) + nSamples - sum(f); % make sure sum(f) == nSamples
        
    otherwise
        error('Unknown emType %s\n',emType);
end

% Set sensor movement positions.
sensor = sensorSet(sensor,'movement positions',[xPos yPos]);
sensor = sensorSet(sensor,'frames per position',f);


end