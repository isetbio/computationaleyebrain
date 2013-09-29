function simParams = constructSimulationParameters(tParams,sParams)
%% function simParams = constructSimulationParameters(tParams, sParams)
%
%  Deal out the parameters structures to create a struct array of
%  parameters, one entry of which will be used in each (potentially
%  parallel) simulation run.
%
%  See also:
%    setParameters
%
%  History:
%    8/24/13  dhb  Pulled this out.
%    9/09/13  hj   re-structure and clean up

%% Check inputs
if nargin < 1, error('theParams is required'); end
if nargin < 2, error('static params structure required'); end

%% Construct grid for dynamic parameters
cdAngles = linspace(0,sParams.dirAngleMax,sParams.nColorDirections+1)';
cdAngles = cdAngles(1:end-1);

os  = 1:length(tParams.OBSERVER_STATES);           % observer states
ct  = 1:length(tParams.DO_TAFC_CLASSIFIER_STATES); % classifier states
mp  = 1:length(tParams.macularPigmentDensityAdjustments);
cdi = 1:sParams.nColorDirections;                  % color direction index
tl  = 0.005 + 0.002 * (1 : sParams.nTestLevels);

[os, ct, mp, cdi, tl] = ndgrid(os, ct, mp, cdi, tl);   % create grid
os = os(:); ct = ct(:); tl = tl(:);
mp = mp(:); cdi = cdi(:); % convert to column array


%% Create simulation structure
simParams = struct(...
    'OBSERVER_STATE', tParams.OBSERVER_STATES(os), ...
    'DO_TAFC_CLASSIFIER', ...
        num2cell(tParams.DO_TAFC_CLASSIFIER_STATES(ct)),...
    'macularPigmentDensityAdjust', ...
        num2cell(tParams.macularPigmentDensityAdjustments(mp)), ...
    'cdAngle',         num2cell(cdAngles(cdi)), ...
    'nTestLevels',     num2cell(tl), ...
    'noiseType',       tParams.noiseType, ...
    'surroundType',    tParams.surroundType, ...
    'surroundSize',    tParams.surroundSize, ...
    'surroundWeight',  tParams.surroundWeight, ...
    'integrationArea', tParams.integrationArea, ...
    'ssFanoFactor',    tParams.secondSiteFanoFactor, ...
    'fractionUse', 0.005 ... % Kluge to select subregion of total cones
);


end