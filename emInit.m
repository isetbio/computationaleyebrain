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
if ~isfield(params, 'center'), params.center = [0 0]; end
if ~isfield(params, 'Sigma')
    warning('Covariance matrix for eye movement missing. Use default(.01)');
    params.Sigma = 0.01 * eye(2);
end
if ~isfield(params, 'nSamples'), params.nSamples = 1000; end
if ~isfield(params, 'fov'), error('Field of view in params required'); end

if isfield(params,'randSeed') % Set up random seed
    rng(randSeed);
else % Store the newly generated random seed
    params.randSeed = rng('shuffle');
end

%% Initialize eye movements
% Each case builds the (x,y) and count variables for every position
emType = ieParamFormat(emType);
switch emType
    case {'fixation'}
        % The eye wanders around the center. The positions are random in a
        % disk around the center. the distances here are in deg of visual
        % angle. 
        % Eye movement is simulated by Brownian motion with increment of
        % Gaussian N(0, Sigma). The initial value is given by params.center
        
        % Compute sensor fov and size
        sz  = sensorGet(sensor, 'size');
        
        % Generate gaussian move
        pos  = ieMvnrnd(params.center, params.Sigma, params.nSamples);
        
        % Adding up to form brownian motion
        % It's a little tricky here. If using brownian motion, we could
        % goes to a very fall distance with high probability. So, we should
        % make it bounce back to center if it gets too large
        %pos = cumsum(pos - repmat(params.center, [params.nSamples,1]))+...
        %       repmat(params.center, [params.nSamples,1]);
        
        % For efficiency, we round the calculations
        % that are centered less than 1 detector's width.
        pos  = round((pos/params.fov).*repmat(sz, [params.nSamples,1])).*...
            params.fov ./ repmat(sz, [params.nSamples, 1]);
        
        % Group the same positions
        [pos,~,ic]  = unique(pos,'rows');
        
        % Compute frame per position
        f    = hist(ic,unique(ic));      % frames per position
        f(1) = f(1) + params.nSamples - sum(f); % make sure sum(f)=nSamples
        
    otherwise
        error('Unknown emType %s\n',emType);
end

% Set sensor movement positions.
sensor = sensorSet(sensor,'movement positions', pos);
sensor = sensorSet(sensor,'frames per position',f);

end