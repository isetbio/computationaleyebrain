function cone = coneCompute(cone, oi)
%% function coneCompute(cone, oi, [nSamples])
%    Compute cone samples for scene / oi
%
%  Inputs:
%    cone       - cone structure, see cone create for more information
%    oi         - Optical image data
%
%  Outputs:
%    cone       - cone structure, with volts image set in the underlying
%                 sensor structure, photon absorptions / volts images can
%                 be fetched with coneGet(cone, 'photons')
%
%  Notes:
%    If need to compute for multiple images, please use command
%      coneSet(cone, 'frames per position', nSamples)
%    Of course, you can use coneSet(cone, 'movement position', [x y]) to
%    handle eye-movements
%
%
%  Example:
%    cone   = coneCreate;
%    sensor = ctConeSamples(cone, oi, 10)
%
%  See also:
%    coneAbsorptions, coneCreate, coneGet
%
%  HJ/BW/DHB (c) ISETBIO Team, 2013

%% Check Inputs
%  Check number of inputs
if notDefined('cone'), error('Cone structure required'); end
if notDefined('oi')
    oi = vcGetSelectedObject('oi');
    if isempty(oi), error('Optical image required'); end
end

%%  Sensor compute oi
sensor = coneGet(cone, 'sensor');
sensor = sensorCompute(sensor, oi);
vcAddAndSelectObject('oi', oi);

%% Compute sensor image samples
%  Compute photon cone absorptions
%  We might need to update coneAbsorption to accept a third parameter to
%  turn off the waitbar. This is only important when we do parallel
%  computation
sensor = coneAbsorptions(sensor, oi);
cone   = coneSet(cone, 'sensor', sensor);

%% Cone adaptation
%  Compute cone adaptation using isetbio routine coneAdaptation
%  I hope the interface of coneAdaptation could be simplified
cone = coneAdaptation(cone);

end